% =========================================================================
% ems_main.m
% Resilient Microgrid Energy Management System – Main Optimisation Script
%
% Purpose : Implements a model predictive control (MPC) scheduler for
%           microgrid energy management using repeated short-horizon LPs.
%           Solved via MATLAB's built-in linprog (Dual-Simplex algorithm).
%
% Workflow:
%   1. Load 96-step forecast and tariff data.
%   2. Use a rolling 48-step lookahead horizon to make decisions at each step.
%   3. Apply only the first control action, then update SOC and repeat.
%   4. Plot the resulting 48-hour MPC schedule.
%
% Dependencies : config_parameters.m, matrix_formulator.m
%
% Author  : [Your Name]
% Date    : 2025
% =========================================================================

clear; clc; close all;

% ── 0. Load system parameters ─────────────────────────────────────────────
run('config_parameters.m');   % populates all constants into workspace

% Package parameters into a struct for clean function passing
params = struct( ...
    'N_steps',      N_steps, ...
    'NV',           NV, ...
    'Time_step',    Time_step, ...
    'E_bess_kWh',   E_bess_kWh, ...
    'eta_ch',       eta_ch, ...
    'eta_dch',      eta_dch, ...
    'SOC_min',      SOC_min, ...
    'SOC_max',      SOC_max, ...
    'P_grid_max',   P_grid_max, ...
    'P_gsell_max',  P_gsell_max, ...
    'P_bess_max',   P_bess_max, ...
    'C_buy',        C_buy, ...
    'C_sell',       C_sell, ...
    'C_p',          C_p, ...
    'IDX_BUY',      IDX_BUY, ...
    'IDX_SELL',     IDX_SELL, ...
    'IDX_CH',       IDX_CH, ...
    'IDX_DCH',      IDX_DCH, ...
    'IDX_SHED',     IDX_SHED, ...
    'Horizon',      Horizon  ...
);

N   = N_steps;   % shorthand: 48 steps/day
NV  = params.NV; % shorthand: 5 variables/step
Dt  = Time_step; % shorthand: 0.5 h

% ── 1. Load Forecast Profiles from CSV ──────────────────────────────────────
%
%   Read GHI, Temp, Load, and Pricing data from microgrid_data.csv
%   Calculate PV power using the physics-informed model.

try
    data = readtable('microgrid_data.csv');
catch
    error('Could not find or read "microgrid_data.csv". Please ensure the file exists in the current directory with the correct 96-row structure.');
end

if height(data) ~= 96
    error('Dataset must contain exactly 96 rows (48 hours at 30-min intervals). Found %d rows.', height(data));
end

% Extract Day 1 Data (Rows 1 to 48)
GHI_D1        = data.GHI_W_m2(1:48);
Temp_D1       = data.Temp_C(1:48);
Load_D1       = data.Load_kW(1:48);
Wind_speed_D1  = data.wind_speed(1:48);
Wind_dir_D1    = data.wind_dir(1:48);

% Extract Day 2 Data (Rows 49 to 96)
GHI_D2        = data.GHI_W_m2(49:96);
Temp_D2       = data.Temp_C(49:96);
Load_D2       = data.Load_kW(49:96);
wind_speed_D2  = data.wind_speed(49:96);
Wind_dir_D2    = data.wind_dir(49:96);


% Extract time-varying electricity tariffs for MPC
C_buy_vec     = data.Buy_Price;
C_sell_vec    = data.Sell_Price;

% Calculate PV Power using the Physics-Informed Model
% Note: Since config_parameters is a script, its variables are already loaded 
% into the workspace. We package the relevant PV constants into a struct.
pv_config = struct('PV_Area', PV_Area, 'PV_Efficiency', PV_Efficiency, ...
                   'Temp_Ref', Temp_Ref, 'Beta_Temp', Beta_Temp, 'NOCT', NOCT);
PV_D1 = pv_physics_model(pv_config, GHI_D1, Temp_D1);
PV_D2 = pv_physics_model(pv_config, GHI_D2, Temp_D2);

% Calculate WTG power using the physics-informed Model
wtg_config = struct('rho',rho,'R',R_wtg,'omega',omega_wtg,...
    'P_rated',P_wtg_rated,'v_cut_in',v_cut_in,'v_cut_out',v_cut_out);

wind_speed = data.wind_speed;
wind_dir   = data.wind_dir;

WTG_D1 = wtg_physics_model(wtg_config, wind_speed(1:48), wind_dir(1:48));
WTG_D2 = wtg_physics_model(wtg_config, wind_speed(49:96), wind_dir(49:96));
WTG_2d = [WTG_D1; WTG_D2];

PV_2d = [PV_D1; PV_D2];

% Update combined generation passed to matrix_formulator:
P_gen_2d = PV_2d + WTG_2d;  % treat as total non-dispatchable generation

% Split total load into Critical (70%) and Non-Critical (30%)
% Adjust these ratios based on actual microgrid requirements.
P_CL_D1    = Load_D1 * 0.70;
P_NCL_D1   = Load_D1 * 0.30;

P_CL_D2    = Load_D2 * 0.70;
P_NCL_D2   = Load_D2 * 0.30;

% ── linprog options ───────────────────────────────────────────────────────
lp_opts = optimoptions('linprog', ...
    'Algorithm',    'dual-simplex', ...
    'Display',      'none',           ...  % suppress per-iteration output
    'OptimalityTolerance', 1e-8,      ...
    'ConstraintTolerance', 1e-8);

% =========================================================================
% MPC – Rolling Horizon Scheduling
% =========================================================================
%
%   Build a 96-step profile and solve a 48-step LP at every interval.
%   The first interval decision is applied in a receding horizon fashion.

PV_2d      = [PV_D1; PV_D2];
P_CL_2d    = [P_CL_D1; P_CL_D2];
P_NCL_2d   = [P_NCL_D1; P_NCL_D2];
C_buy_2d   = C_buy_vec;
C_sell_2d  = C_sell_vec;
T_total    = numel(PV_2d);
Horizon    = params.Horizon;          % configured in config_parameters.m

fprintf('─────────────────────────────────────────────\n');
fprintf('MPC : rolling horizon control with H = %d steps\n', Horizon);
fprintf('─────────────────────────────────────────────\n');

P_gb      = zeros(T_total, 1);
P_gs      = zeros(T_total, 1);
P_ch      = zeros(T_total, 1);
P_dch     = zeros(T_total, 1);
P_shed    = zeros(T_total, 1);
SOC_2d    = zeros(T_total, 1);

alpha = params.eta_ch * Dt / params.E_bess_kWh;
beta  = Dt / (params.eta_dch * params.E_bess_kWh);

SOC_cur = SOC_init;
realized_cost = 0;

for t = 1:T_total
    h = min(Horizon, T_total - t + 1);

    P_pv_h   = PV_2d(t : t+h-1);
    P_CL_h   = P_CL_2d(t : t+h-1);
    P_NCL_h  = P_NCL_2d(t : t+h-1);

    params.C_buy  = C_buy_2d(t : t+h-1);
    params.C_sell = C_sell_2d(t : t+h-1);

    [f_h, A_h, b_h, Aeq_h, beq_h, lb_h, ub_h] = matrix_formulator( ...
        P_pv_h, P_CL_h, P_NCL_h, params, SOC_cur, NaN);

    [x_h, cost_h, exitflag_h, output_h] = linprog(f_h, A_h, b_h, Aeq_h, beq_h, lb_h, ub_h, lp_opts);

    if exitflag_h ~= 1
        warning('MPC horizon at step %d did not converge (exitflag = %d).', t, exitflag_h);
        x_h = zeros(numel(f_h), 1);
    end

    X_h = reshape(x_h, NV, h).';

    P_gb(t)   = X_h(1, IDX_BUY);
    P_gs(t)   = X_h(1, IDX_SELL);
    P_ch(t)   = X_h(1, IDX_CH);
    P_dch(t)  = X_h(1, IDX_DCH);
    P_shed(t) = X_h(1, IDX_SHED);

    SOC_cur = SOC_cur + alpha * P_ch(t) - beta * P_dch(t);
    SOC_2d(t) = SOC_cur;

    realized_cost = realized_cost + (C_buy_2d(t) * P_gb(t) ...
                       - C_sell_2d(t) * P_gs(t) ...
                       + params.C_p * P_shed(t)) * Dt;
end

P_net_grid_2d = P_gb - P_gs;   % [96×1]
P_net_bess_2d = P_dch - P_ch;  % [96×1] >0 discharge, <0 charge

t_2day = Dt * (1 : T_total).';

fprintf('=== MPC Summary ===\n');
fprintf('  Total realised cost : %.2f Yen\n', realized_cost);
fprintf('  Total NCL shed      : %.3f kWh\n', sum(P_shed) * Dt);
fprintf('  Final SOC           : %.4f\n\n', SOC_2d(end));

% =========================================================================
% PLOTS – 48-Hour Horizon Results
% =========================================================================

figure('Name','Microgrid EMS – MPC Optimisation Results', ...
       'Position', [100 100 1100 750], 'Color','w');

% ── Common formatting closure ─────────────────────────────────────────────
fmt = @(ax, ylbl) set(ax, 'XLim',[0 48], 'XTick',0:6:48, ...
    'XGrid','on', 'YGrid','on', 'FontSize',11, 'Box','on');

x_tick_labels = {'0','6','12','18','24','30','36','42','48'};

% ── Subplot 1 : BESS State-of-Charge ────────────────────────────────────
ax1 = subplot(3, 1, 1);

% Pre-pend the initial conditions so the curve starts at t = 0
t_plot  = [0; t_2day];
SOC_plot = [SOC_init; SOC_2d];

area(t_plot, SOC_plot * 100, 'FaceColor',[0.18 0.55 0.80], ...
    'FaceAlpha',0.35, 'EdgeColor',[0.10 0.35 0.65], 'LineWidth',1.5);
hold on;
yline(SOC_min * 100, '--r', 'LineWidth',1.4, 'Label','SOC_{min}=10%', ...
    'LabelHorizontalAlignment','left');
yline(SOC_max * 100, '--', 'Color',[0.2 0.7 0.2], 'LineWidth',1.4, ...
    'Label','SOC_{max}=90%', 'LabelHorizontalAlignment','left');
xline(24, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5, ...
    'Label','Day boundary', 'LabelVerticalAlignment','bottom');
% No fixed handoff marker is shown for MPC results
fmt(ax1, 'SOC [%]');
ylim([0 100]); ylabel('SOC [%]', 'FontWeight','bold');
title('BESS State-of-Charge over 48-Hour Horizon', 'FontSize',13);
legend({'SOC trajectory','SOC_{min}','SOC_{max}','Day boundary','SOC handoff target'}, ...
    'Location','southeast','FontSize',9);
xticklabels(x_tick_labels);

% ── Subplot 2 : Net Grid Power ───────────────────────────────────────────
ax2 = subplot(3, 1, 2);

P_pos = max(P_net_grid_2d,  0);   % import  (positive part)
P_neg = min(P_net_grid_2d,  0);   % export  (negative part)

bar(t_2day, P_pos, 'FaceColor',[0.85 0.33 0.10], ...
    'EdgeColor','none', 'BarWidth',1);
hold on;
bar(t_2day, P_neg, 'FaceColor',[0.47 0.67 0.19], ...
    'EdgeColor','none', 'BarWidth',1);
yline(0, 'k', 'LineWidth',0.8);
xline(24, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

fmt(ax2, 'P_{grid} [kW]');
ylabel('Net Grid Power [kW]', 'FontWeight','bold');
title('Net Grid Power (positive = import, negative = export)', 'FontSize',13);
legend({'Grid Import','Grid Export'}, 'Location','northeast','FontSize',9);
xticklabels(x_tick_labels);

% ── Subplot 3 : Net BESS Power ───────────────────────────────────────────
ax3 = subplot(3, 1, 3);

B_dch = max(P_net_bess_2d,  0);  % discharge (positive)
B_ch  = min(P_net_bess_2d,  0);  % charge    (negative)

bar(t_2day, B_dch, 'FaceColor',[0.93 0.69 0.13], ...
    'EdgeColor','none', 'BarWidth',1);
hold on;
bar(t_2day, B_ch,  'FaceColor',[0.30 0.75 0.93], ...
    'EdgeColor','none', 'BarWidth',1);
yline(0, 'k', 'LineWidth',0.8);
xline(24, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

fmt(ax3, 'P_{BESS} [kW]');
ylabel('Net BESS Power [kW]', 'FontWeight','bold');
title('Net BESS Power (positive = discharge, negative = charge)', 'FontSize',13);
legend({'BESS Discharge','BESS Charge'}, 'Location','northeast','FontSize',9);
xlabel('Hour of Horizon [h]', 'FontSize',11, 'FontWeight','bold');
xticklabels(x_tick_labels);

% ── Shared x-axis alignment ───────────────────────────────────────────────
linkaxes([ax1 ax2 ax3], 'x');
xlim([0 48]);

sgtitle(sprintf('Resilient Microgrid EMS – MPC Optimisation\nTotal Cost: %.1f Yen', realized_cost), ...
    'FontSize', 14, 'FontWeight', 'bold');

% ── Save figure ───────────────────────────────────────────────────────────
set(gcf, 'Toolbar', 'none');
print(gcf, 'ems_results_mpc.png', '-dpng', '-r150');
fprintf('\nFigure saved as ems_results_mpc.png\n');

% =========================================================================
% PV vs Load vs Wind Power Comparison
% =========================================================================
figure('Name','PV vs Load vs Wind Power', 'Position', [150 150 900 420], 'Color', 'w');
plot(t_2day, PV_2d, '-o', 'LineWidth', 1.6, 'MarkerSize', 3, 'DisplayName', 'PV generation');
hold on;
plot(t_2day, WTG_2d, '-^', 'LineWidth', 1.6, 'MarkerSize', 3, 'DisplayName', 'Wind turbine power');
plot(t_2day, P_CL_2d + P_NCL_2d, '-s', 'LineWidth', 1.6, 'MarkerSize', 3, 'DisplayName', 'Total load');
xline(24, ':k', 'Day boundary', 'LabelHorizontalAlignment','left');

grid on;
xlabel('Hour of Horizon [h]');
ylabel('Power [kW]');
title('PV Generation vs Wind Power vs Total Load');
legend('Location','best');
set(gca, 'XLim', [0 48], 'XTick', 0:6:48, 'XTickLabel', x_tick_labels, 'FontSize', 11);
set(gcf, 'Toolbar', 'none');
print(gcf, 'ems_results_pv_vs_load_vs_wind.png', '-dpng', '-r150');
fprintf('Figure saved as ems_results_pv_vs_load_vs_wind.png\n');

% =========================================================================
% OPTIONAL: Print per-interval decision table (first 10 intervals)
% =========================================================================
fprintf('\n── MPC Schedule (first 10 intervals) ───────────────────────────\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-6s\n', ...
    'Step','P_gb','P_gs','P_ch','P_dch','P_shed','SOC');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-6s\n', ...
    '','[kW]','[kW]','[kW]','[kW]','[kW]','[-]');
for t = 1:min(10, T_total)
    fprintf('%-6d %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-6.4f\n', ...
        t, P_gb(t), P_gs(t), P_ch(t), P_dch(t), ...
        P_shed(t), SOC_2d(t));
end
fprintf('...  (%d total steps)\n', T_total);
