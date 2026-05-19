% =========================================================================
% matrix_formulator.m
% Resilient Microgrid EMS – LP Matrix Formulator
%
% Purpose : Construct all LP matrices (f, A, b, Aeq, beq, lb, ub) for a
%           single 24-hour block (N_steps = 48 half-hour intervals) given
%           the forecast profiles and system parameters.
%
% Syntax  :
%   [f, A, b, Aeq, beq, lb, ub] = matrix_formulator( ...
%       P_pv, P_load_CL, P_load_NCL, params, SOC_start, SOC_end_target)
%
% Inputs  :
%   P_pv          – [N_steps × 1]  Forecasted PV power           [kW]
%   P_load_CL     – [N_steps × 1]  Critical load profile          [kW]
%   P_load_NCL    – [N_steps × 1]  Non-critical load profile      [kW]
%   params        – struct with all fields from config_parameters.m
%   SOC_start     – scalar  Initial SOC at t=0 of this day block   [-]
%   SOC_end_target– scalar  Target terminal SOC (NaN → unconstrained)
%
% Outputs :
%   f    – [N_tot × 1]       LP cost (objective) vector
%   A    – [n_ineq × N_tot]  Inequality constraint matrix  (A*x ≤ b)
%   b    – [n_ineq × 1]      Inequality RHS vector
%   Aeq  – [n_eq   × N_tot]  Equality constraint matrix   (Aeq*x = beq)
%   beq  – [n_eq   × 1]      Equality RHS vector
%   lb   – [N_tot × 1]       Lower bounds on decision variables
%   ub   – [N_tot × 1]       Upper bounds on decision variables
%
% Decision-variable layout (N_tot = N_steps × NV = 48 × 5 = 240):
%   For each time step t = 1 … N_steps the five variables are stacked:
%
%     global index = (t-1)*NV + local_offset
%
%   local_offset  Variable            Symbol          Unit
%   ──────────────────────────────────────────────────────
%        1        P_grid_buy(t)       P_gb(t)         kW
%        2        P_grid_sell(t)      P_gs(t)         kW
%        3        P_BESS_ch(t)        P_ch(t)         kW
%        4        P_BESS_dch(t)       P_dch(t)        kW
%        5        P_load_NCL_shed(t)  P_shed(t)       kW
%
% Power balance at time t (enforced as LP equality row):
%   P_gb(t) + P_pv(t) + P_dch(t)
%       = P_load_CL(t) + P_load_NCL(t) - P_shed(t) + P_ch(t) + P_gs(t)
%
%   Rearranged (moving all variables to LHS):
%   P_gb(t) - P_gs(t) - P_ch(t) + P_dch(t) + P_shed(t)
%       = P_load_CL(t) + P_load_NCL(t) - P_pv(t)
%         ────────────── net_demand(t) ──────────────
%
% SOC dynamics (Wh-balance, enforced as LP equality row for t = 1…N_steps):
%   E_bess * SOC(t) = E_bess * SOC(t-1)
%                   + eta_ch  * P_ch(t)  * Δt
%                   - (1/eta_dch) * P_dch(t) * Δt
%
%   SOC(0) ≡ SOC_start  (initial condition substituted directly)
%
%   Substituting and defining:
%     α = eta_ch  * Δt / E_bess
%     β = (1/eta_dch) * Δt / E_bess
%
%   The cumulative SOC after n steps:
%     SOC(n) = SOC_start + α * Σ_{t=1..n} P_ch(t) − β * Σ_{t=1..n} P_dch(t)
%
%   For the LP we write SOC(t) explicitly at every step using running sum
%   constraints rather than recursive substitution – this keeps Aeq sparse
%   and block-diagonal in structure.
%
% =========================================================================

function [f, A, b, Aeq, beq, lb, ub] = matrix_formulator( ...
    P_pv, P_load_CL, P_load_NCL, params, SOC_start, SOC_end_target)

% ── Unpack parameters ─────────────────────────────────────────────────────
N   = numel(P_pv);         % horizon length may vary for MPC
NV  = params.NV;            % 5 variables per interval
Dt  = params.Time_step;     % 0.5 h

E   = params.E_bess_kWh;    % BESS energy capacity [kWh]
eta_c   = params.eta_ch;    % Charging efficiency
eta_d   = params.eta_dch;   % Discharging efficiency
SOC_lo  = params.SOC_min;   % Lower SOC bound
SOC_hi  = params.SOC_max;   % Upper SOC bound

Pgb_max  = params.P_grid_max;   % Max grid import  [kW]
Pgs_max  = params.P_gsell_max;  % Max grid export  [kW]
Pb_max   = params.P_bess_max;   % Max BESS power   [kW]

C_buy    = params.C_buy;        % [Yen/kWh], scalar or horizon-length vector
C_sell   = params.C_sell;       % [Yen/kWh], scalar or horizon-length vector
C_p      = params.C_p;          % [Yen/kWh] NCL shedding penalty

if isscalar(C_buy)
    C_buy = repmat(C_buy, N, 1);
end
if isscalar(C_sell)
    C_sell = repmat(C_sell, N, 1);
end

% Local index offsets (within a single time step's variable block)
iB  = params.IDX_BUY;   % 1 → P_grid_buy
iS  = params.IDX_SELL;  % 2 → P_grid_sell
iC  = params.IDX_CH;    % 3 → P_BESS_ch
iD  = params.IDX_DCH;   % 4 → P_BESS_dch
iSh = params.IDX_SHED;  % 5 → P_load_NCL_shed

N_tot = N * NV;   % Total number of decision variables = 240

% Helper: global index of variable v at step t
%   gidx(t, v) = (t-1)*NV + v
gidx = @(t, v) (t-1)*NV + v;

% =========================================================================
% 1.  OBJECTIVE VECTOR  f  (N_tot × 1)
% =========================================================================
%
%   min  Σ_t [ C_buy*Δt * P_gb(t) - C_sell*Δt * P_gs(t) + C_p*Δt * P_shed(t) ]
%
%   f(gidx(t, iB))  =  C_buy  * Δt   (penalise import)
%   f(gidx(t, iS))  = -C_sell * Δt   (reward export)
%   f(gidx(t, iC))  =  0              (charging has no direct cost)
%   f(gidx(t, iD))  =  0              (discharging has no direct cost)
%   f(gidx(t, iSh)) =  C_p   * Δt   (penalise load shedding)

f = zeros(N_tot, 1);
for t = 1:N
    f(gidx(t, iB))  =  C_buy(t)  * Dt;   % grid purchase cost per step
    f(gidx(t, iS))  = -C_sell(t) * Dt;   % grid sell revenue per step (negative = gain)
    f(gidx(t, iC))  =  0;                % BESS charge: cost captured via grid term
    f(gidx(t, iD))  =  0;                % BESS discharge: cost captured via grid term
    f(gidx(t, iSh)) =  C_p   * Dt;      % NCL shedding penalty per step
end

% =========================================================================
% 2.  LOWER AND UPPER BOUNDS  lb, ub  (N_tot × 1 each)
% =========================================================================
%
%   All decision variables are non-negative (rectified power quantities).
%
%   Variable        lb      ub
%   ─────────────────────────────────────────────────────
%   P_grid_buy(t)    0      P_grid_max
%   P_grid_sell(t)   0      P_gsell_max
%   P_BESS_ch(t)     0      P_bess_max
%   P_BESS_dch(t)    0      P_bess_max
%   P_NCL_shed(t)    0      P_load_NCL(t)   ← cannot shed more than exists

lb = zeros(N_tot, 1);
ub = zeros(N_tot, 1);

for t = 1:N
    ub(gidx(t, iB))  = Pgb_max;
    ub(gidx(t, iS))  = Pgs_max;
    ub(gidx(t, iC))  = Pb_max;
    ub(gidx(t, iD))  = Pb_max;
    ub(gidx(t, iSh)) = P_load_NCL(t);  % cannot shed more NCL than scheduled
end

% =========================================================================
% 3.  EQUALITY CONSTRAINTS  Aeq * x = beq
% =========================================================================
%
%   Two sets of equality rows:
%     (a) Power balance at each step          →  N rows
%     (b) SOC dynamics at each step           →  N rows
%         (optional) Terminal SOC constraint  →  1 extra row if requested
%
%   Total equality rows: 2*N  (+ 1 if SOC_end_target is finite)

n_soc_terminal = double(~isnan(SOC_end_target));  % 1 or 0
n_eq = 2*N + n_soc_terminal;

Aeq  = zeros(n_eq, N_tot);
beq  = zeros(n_eq, 1);

% ── (a) Power Balance Rows  (rows 1 … N) ─────────────────────────────────
%
%   Physical equation at step t:
%     P_gb(t) + P_pv(t) + P_dch(t)               ← supply side
%         = P_CL(t) + P_NCL(t) - P_shed(t)        ← net demand
%           + P_ch(t) + P_gs(t)                    ← storage charge + export
%
%   LP row (all variables on LHS, constants on RHS):
%     +1 * P_gb(t)
%     -1 * P_gs(t)
%     -1 * P_ch(t)
%     +1 * P_dch(t)
%     +1 * P_shed(t)
%     = P_CL(t) + P_NCL(t) - P_pv(t)   ≡ net_demand(t)
%
for t = 1:N
    row = t;   % row index for power balance at step t
    Aeq(row, gidx(t, iB))  =  1;   % P_gb  contributes to supply
    Aeq(row, gidx(t, iS))  = -1;   % P_gs  diverts supply to grid
    Aeq(row, gidx(t, iC))  = -1;   % P_ch  diverts supply to battery
    Aeq(row, gidx(t, iD))  =  1;   % P_dch adds to supply
    Aeq(row, gidx(t, iSh)) =  1;   % P_shed reduces net demand (positive = demand reduced)

    % RHS: net demand = load - PV generation
    %   (positive → more demand than PV; negative → PV surplus)
    beq(row) = P_load_CL(t) + P_load_NCL(t) - P_pv(t);
end

% ── (b) SOC Dynamics Rows  (rows N+1 … 2*N) ─────────────────────────────
%
%   Wh-balance for step t (half-hour increment):
%     E * SOC(t) = E * SOC(t-1)
%                + eta_c  * P_ch(t)  * Δt
%                - (1/eta_d) * P_dch(t) * Δt
%
%   Substituting SOC(t-1) recursively down to SOC(0) = SOC_start:
%
%     E * SOC(t) = E * SOC_start
%                + eta_c  * Δt * Σ_{k=1..t} P_ch(k)
%                - (1/eta_d) * Δt * Σ_{k=1..t} P_dch(k)
%
%   LP row (move variable sums to LHS, constant E*SOC_start to RHS):
%     eta_c/E  * Δt * Σ_{k=1..t} P_ch(k)
%   - 1/(E*eta_d) * Δt * Σ_{k=1..t} P_dch(k)
%   = SOC(t) - SOC_start
%
%   To avoid introducing auxiliary SOC variables we write:
%   For each t, the row enforces that the running energy change
%   equals (SOC(t) - SOC_start) * E.
%
%   Since SOC(t) must stay in [SOC_lo, SOC_hi], these bounds are
%   translated to equivalent energy bounds and placed in the
%   inequality section (see Section 4 below).
%
%   Here we only encode the SOC-update equality for each t.
%
alpha = eta_c  * Dt / E;   % Energy credited per kW of charging
beta  = Dt / (eta_d * E);  % Energy debited  per kW of discharging

for t = 1:N
    row = N + t;   % row index for SOC dynamics at step t

    % Running sum: fill non-zero entries for k = 1 … t
    for k = 1:t
        Aeq(row, gidx(k, iC))  =  alpha;   %  + eta_c*Δt/E per kW charged
        Aeq(row, gidx(k, iD))  = -beta;    %  - Δt/(eta_d*E) per kW discharged
    end

    % RHS: desired cumulative SOC increment from start
    %   We treat SOC(t) as a FREE quantity here – it is bounded through
    %   the inequality constraints in Section 4.  The equality rows only
    %   define the SOC trajectory; the inequality rows bound it.
    %
    %   To avoid a circular dependency we leave beq blank (= 0) for the
    %   dynamics rows and instead reformulate:
    %
    %     0 ≤ SOC_start + (row accumulation) ≤ SOC_hi
    %
    %   This means the equality rows here serve ONLY as the inequality
    %   seeds below.  We therefore REMOVE these rows from Aeq and instead
    %   encode the SOC bounds as inequalities.
    %
    %   ── Implementation choice ──────────────────────────────────────────
    %   Rather than sparse auxiliary variables, we use the cleaner
    %   approach of putting SOC bounds in the A/b system (Section 4).
    %   The Aeq rows N+1…2N are CLEARED here and repopulated only if a
    %   terminal SOC equality is requested.
    %   ───────────────────────────────────────────────────────────────────
    Aeq(row, :) = 0;   % cleared; SOC bounds → inequality section
    beq(row)    = 0;
end

% ── (c) Terminal SOC Equality  (row 2*N + 1, optional) ───────────────────
%
%   Forces SOC at end of day to equal SOC_end_target:
%     SOC_start + alpha * Σ_all P_ch(t) - beta * Σ_all P_dch(t) = SOC_target
%
%   LP row:
%     alpha * Σ_{t=1..N} P_ch(t) - beta * Σ_{t=1..N} P_dch(t)
%     = SOC_end_target - SOC_start

if n_soc_terminal
    row = 2*N + 1;
    for t = 1:N
        Aeq(row, gidx(t, iC)) =  alpha;
        Aeq(row, gidx(t, iD)) = -beta;
    end
    beq(row) = SOC_end_target - SOC_start;
end

% Remove the all-zero SOC dynamics rows we cleared above
%   (keeps Aeq as compact as possible for linprog)
keep_eq = any(Aeq ~= 0, 2) | (beq ~= 0);
Aeq = Aeq(keep_eq, :);
beq = beq(keep_eq);

% =========================================================================
% 4.  INEQUALITY CONSTRAINTS  A * x ≤ b
% =========================================================================
%
%   (a) SOC lower bound at each step t:
%         SOC_start + alpha*Σ_{k≤t} P_ch - beta*Σ_{k≤t} P_dch  ≥  SOC_lo
%       → -alpha*Σ P_ch + beta*Σ P_dch  ≤  SOC_start - SOC_lo
%
%   (b) SOC upper bound at each step t:
%         SOC_start + alpha*Σ_{k≤t} P_ch - beta*Σ_{k≤t} P_dch  ≤  SOC_hi
%       → +alpha*Σ P_ch - beta*Σ P_dch  ≤  SOC_hi - SOC_start
%
%   (c) Mutual exclusion of simultaneous charge/discharge is NOT enforced
%       via LP (requires binary variables).  Instead, the cost function
%       naturally prevents simultaneous charge/discharge when C_buy > 0.
%
%   Total inequality rows: 2*N

n_ineq = 2 * N;
A = zeros(n_ineq, N_tot);
b = zeros(n_ineq, 1);

for t = 1:N
    % ── SOC lower-bound row  (rows 1 … N) ──────────────────────────────
    row_lo = t;
    for k = 1:t
        A(row_lo, gidx(k, iC)) = -alpha;   % -eta_c*Δt/E per kW charge
        A(row_lo, gidx(k, iD)) =  beta;    % +Δt/(eta_d*E) per kW discharge
    end
    b(row_lo) = SOC_start - SOC_lo;        % SOC_start - SOC_min

    % ── SOC upper-bound row  (rows N+1 … 2*N) ──────────────────────────
    row_hi = N + t;
    for k = 1:t
        A(row_hi, gidx(k, iC)) =  alpha;   % +eta_c*Δt/E per kW charge
        A(row_hi, gidx(k, iD)) = -beta;    % -Δt/(eta_d*E) per kW discharge
    end
    b(row_hi) = SOC_hi - SOC_start;        % SOC_max - SOC_start
end

% =========================================================================
% 5.  DIMENSION CHECKS (diagnostic guard)
% =========================================================================
assert(length(f)  == N_tot, 'f dimension mismatch');
assert(size(A,2)  == N_tot, 'A column dimension mismatch');
assert(size(Aeq,2) == N_tot || isempty(Aeq), 'Aeq column dimension mismatch');
assert(length(lb) == N_tot, 'lb dimension mismatch');
assert(length(ub) == N_tot, 'ub dimension mismatch');

fprintf('  [matrix_formulator] Matrices built: N_tot=%d, n_ineq=%d, n_eq=%d\n', ...
    N_tot, size(A,1), size(Aeq,1));

end  % function matrix_formulator
