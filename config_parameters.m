% =========================================================================
% config_parameters.m
% Resilient Microgrid Energy Management System – System Constants
%
% Purpose : Define and return all physical, economic, and operational
%           parameters used throughout the LP-based EMS optimisation.
%           Values are taken directly from the reference paper baseline.
%
% Usage   : Call this script (or run it) before matrix_formulator.m or
%           ems_main.m so that all variables exist in the caller workspace.
%
% Reference: "Resilient Microgrid EMS via Two-Stage LP Optimisation"
%            (see paper Table I and Section III-A for parameter sources)
% =========================================================================

% ── Temporal discretisation ──────────────────────────────────────────────
Time_step   = 0.5;          % [h]  Length of each decision interval (30 min)
N_steps     = 48;           % [-]  Number of intervals per 24-hour day
t_vec       = (1:N_steps).' * Time_step; % [h]  Cumulative hour vector (0.5 … 24)

% ── PV Generator ─────────────────────────────────────────────────────────
P_pv_max    = 100;          % [kW] Rated DC capacity of the PV array
eta_pv      = 1.0;          % [-]  Inverter + coupling efficiency (lumped)
%  Note: The actual per-interval PV profile P_pv(t) is generated in
%        ems_main.m as a synthetic curtailable forecast; P_pv_max sets
%        the upper bound for every interval.

% --- Wind Turbine Generator ─────────────────────────────────
rho        = 1.225;   % kg/m³
R_wtg      = 2.5;     % rotor radius (m) — tune to get 2 kW rated
omega_wtg  = 4.0;     % rad/s
v_cut_in   = 3.0;     % m/s
v_cut_out  = 25.0;    % m/s
P_wtg_rated = 2.0;    % kW

% ── Battery Energy Storage System (BESS) ─────────────────────────────────
BESS_Ah     = 82;           % [Ah]  Rated capacity at C/10 rate
V_dc        = 400;          % [V]   DC bus voltage (nominal)
E_bess_kWh  = (BESS_Ah * V_dc) / 1000; % [kWh] = 82 × 400 / 1000 = 32.8 kWh

SOC_min     = 0.10;         % [-]  Minimum allowable state-of-charge
SOC_max     = 0.90;         % [-]  Maximum allowable state-of-charge
SOC_init    = 1.00;         % [-]  SOC at the start of the planning horizon
%  Usable energy window:
E_usable    = (SOC_max - SOC_min) * E_bess_kWh; % [kWh] = 26.24 kWh

eta_ch      = 0.95;         % [-]  BESS charging efficiency  (one-way)
eta_dch     = 0.95;         % [-]  BESS discharging efficiency (one-way)

P_bess_max  = 10.0;         % [kW] Maximum charge / discharge power
%  The symmetric power limit is applied separately to P_BESS_ch and
%  P_BESS_dch via the lb/ub vectors in matrix_formulator.m.

% ── Grid Connection ───────────────────────────────────────────────────────
P_grid_max  = 10.0;         % [kW] Maximum import power from utility
P_gsell_max = 10.0;         % [kW] Maximum export power to utility

% ── Load Profiles ─────────────────────────────────────────────────────────
%  Critical load (CL)   : must be served at all times – no shedding allowed.
%  Non-Critical load (NCL): may be partially shed subject to penalty C_p.
%  Both profiles are generated synthetically in ems_main.m.
P_cl_max    = 5.0;          % [kW] Peak critical load (reference magnitude)
P_ncl_max   = 3.0;          % [kW] Peak non-critical load (reference magnitude)

% ── Economic Parameters ──────────────────────────────────────────────────
C_p         = 5;            % [Yen/kWh] Penalty cost for NCL shedding
%  Grid tariff vectors (flat rates used here; replace with TOU if needed):
C_buy       = 20;           % [Yen/kWh] Unit cost of purchasing from grid
C_sell      = 8;            % [Yen/kWh] Unit revenue from selling to grid
%  Combined cost vector per unit energy exchanged in one Time_step:
%    cost_buy  = C_buy  * Time_step   [Yen per kW of P_grid_buy  per step]
%    cost_sell = C_sell * Time_step   [Yen per kW of P_grid_sell per step]
%    cost_shed = C_p    * Time_step   [Yen per kW of P_shed      per step]

% ── Two-Stage SOC Target (Paper Baseline) ────────────────────────────────
SOC_target_D1_end = 0.45;   % [-]  Target SOC at end of Day 1 (= start Day 2)
%  This value comes from the Day-2 optimisation (Step 1 in ems_main.m).
%  The paper's baseline sets it to 0.45; ems_main.m overwrites this
%  programmatically once the Day-2 LP is solved.

% ── Variable Index Mapping (per time step) ───────────────────────────────
%  Each time step t contributes NV = 5 decision variables arranged as:
%
%  Index within step  |  Variable         |  Unit  |  Physical meaning
%  ─────────────────────────────────────────────────────────────────────
%   1  (IDX_BUY)      |  P_grid_buy(t)    |  kW    |  Power imported from grid
%   2  (IDX_SELL)     |  P_grid_sell(t)   |  kW    |  Power exported to grid
%   3  (IDX_CH)       |  P_BESS_ch(t)     |  kW    |  BESS charging power
%   4  (IDX_DCH)      |  P_BESS_dch(t)    |  kW    |  BESS discharging power
%   5  (IDX_SHED)     |  P_load_NCL_shed(t)|  kW   |  Non-critical load shed
%
%  Power balance (enforced as equality at every t):
%    P_grid_buy(t) + P_pv(t) + P_BESS_dch(t)
%        = P_load_CL(t) + P_load_NCL(t) - P_load_NCL_shed(t)
%          + P_BESS_ch(t) + P_grid_sell(t)
%
NV          = 5;            % Number of decision variables per time step
IDX_BUY     = 1;            % Column offset → P_grid_buy
IDX_SELL    = 2;            % Column offset → P_grid_sell
IDX_CH      = 3;            % Column offset → P_BESS_ch
IDX_DCH     = 4;            % Column offset → P_BESS_dch
IDX_SHED    = 5;            % Column offset → P_load_NCL_shed

% ── MPC Parameters ───────────────────────────────────────────────────────
Horizon     = 48;           % [steps] Lookahead horizon for MPC (48 = 24 h)

% ── PV Physics Model Parameters ──────────────────────────────────────────
PV_Area         = 40;            % [m^2] Total area of the PV array
PV_Efficiency   = 0.1875;        % [-] Standard efficiency of the PV panels
Temp_Ref        = 25;            % [°C] Standard Test Condition (STC) reference temperature
Beta_Temp       = -0.004;        % [1/°C] Temperature coefficient of power
NOCT            = 45;            % [°C] Nominal Operating Cell Temperature

% ── Display Summary ──────────────────────────────────────────────────────
fprintf('=== Microgrid EMS – Configuration Loaded ===\n');
fprintf('  BESS capacity  : %.2f kWh  (usable: %.2f kWh)\n', E_bess_kWh, E_usable);
fprintf('  SOC window     : [%.2f, %.2f]  |  SOC_init = %.2f\n', SOC_min, SOC_max, SOC_init);
fprintf('  Time step      : %.1f h  (%d steps/day)\n', Time_step, N_steps);
fprintf('  C_buy / C_sell : %d / %d Yen/kWh  |  C_penalty = %d Yen/kWh\n', C_buy, C_sell, C_p);
fprintf('=============================================\n\n');
