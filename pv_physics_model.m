function P_pv = pv_physics_model(config, GHI, Temp_amb)
% PV_PHYSICS_MODEL - Calculates PV power output from GHI and Ambient Temperature
%
% Inputs:
%   config   - Struct containing PV physical parameters (Area, Efficiency, etc.)
%   GHI      - [W/m^2] Global Horizontal Irradiance (Vector)
%   Temp_amb - [°C] Ambient air temperature (Vector)
%
% Outputs:
%   P_pv     - [kW] Calculated PV power output (Vector)

% 1. Calculate Cell Temperature (T_cell)
% Uses the standard NOCT (Nominal Operating Cell Temperature) model
% T_cell = T_amb + (NOCT - 20) / 800 * GHI
T_cell = Temp_amb + ((config.NOCT - 20) / 800) .* GHI;

% 2. Calculate Temperature Degradation Factor
% Factors in power loss due to temperatures exceeding the 25°C STC
% Derating_Factor = 1 + Beta_Temp * (T_cell - Temp_Ref)
Derating_Factor = 1 + config.Beta_Temp .* (T_cell - config.Temp_Ref);

% 3. Calculate Final PV Power Output
% P_pv = GHI * Area * Base_Efficiency * Derating_Factor
% Divided by 1000 to convert Watts to kiloWatts (kW)
P_pv = (GHI .* config.PV_Area .* config.PV_Efficiency .* Derating_Factor) / 1000;

% Ensure power output does not go negative
P_pv = max(0, P_pv);

end
