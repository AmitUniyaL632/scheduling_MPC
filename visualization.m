% clear all; clc;
% =========================================================================
%Plotting the Load and grid economics
% =========================================================================
Fig_1 = figure(Name = "Load & Grid Economics");
y1 = Load_profile;
yyaxis left
% area(y1, 'FaceColor',[0.18 0.55 0.80], ...
%     'FaceAlpha',0.35, 'EdgeColor',[0.10 0.35 0.65], 'LineWidth',1.5);
plot(y1,'LineWidth', 1, 'LineStyle','-','color', 'b',...
    'Marker','square','MarkerFaceColor',[0,1,1], ...
    'MarkerEdgeColor',[0,0,1], 'MarkerSize', 3);
ylabel("Load Profile (kW)");
xlabel('Sample Instant');

y2 = grid_import_price;
y3 = grid_export_price;
yyaxis right
hold on
plot(y2,'LineWidth', 1, 'LineStyle','-','color', 'r',...
    'Marker','square','MarkerFaceColor',[245,112,125]./255, ...
    'MarkerEdgeColor',[255,0,25]./255, 'MarkerSize', 3);
plot(y3, 'LineWidth', 1, 'LineStyle','-','color', 'g',...
    'Marker','square','MarkerFaceColor',[231,255,235]./255, ...
    'MarkerEdgeColor',[0,79,22]./255, 'MarkerSize', 3);

ylabel('Price (₹/kW)');
hold off;
xline(48, ':k', 'Day boundary', 'LabelHorizontalAlignment','left');

title('Load and Grid Prices');
legend('Load Profile', 'Grid import prices', 'Grid export prices');
grid on;


% =========================================================================
% PLOTS – 48-Hour Horizon Results
% =========================================================================

Fig_2 = figure('Name','Microgrid EMS – MPC Optimisation Results');

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

sgtitle(sprintf('Resilient Microgrid EMS – MPC Optimisation\nTotal Cost: %.1f Rupees', realized_cost), ...
    'FontSize', 14, 'FontWeight', 'bold');

% ── Save figure ───────────────────────────────────────────────────────────
% set(gcf, 'Toolbar', 'none');
% print(gcf, 'ems_results_mpc.png', '-dpng', '-r150');
% fprintf('\nFigure saved as ems_results_mpc.png\n');

% =========================================================================
% PV vs Load vs Wind Power Comparison
% =========================================================================
Fig_3 = figure('Name','PV vs Load vs Wind Power');
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
% set(gcf, 'Toolbar', 'none');
% print(gcf, 'ems_results_pv_vs_load_vs_wind.png', '-dpng', '-r150');
% fprintf('Figure saved as ems_results_pv_vs_load_vs_wind.png\n');

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
