function P_wtg = wtg_physics_model(config, wind_speed, wind_dir)
% Inputs:
%   config     - struct with turbine parameters
%   wind_speed - [N×1] wind speed (m/s)
%   wind_dir   - [N×1] wind direction (degrees), optional

rho    = config.rho;          % air density (kg/m³) ~1.225
R      = config.R;            % rotor radius (m)
A      = pi * R^2;            % swept area
P_rated = config.P_rated;     % 2 kW
v_ci   = config.v_cut_in;     % ~3 m/s
v_co   = config.v_cut_out;    % ~25 m/s
omega  = config.omega;        % rotor angular speed (rad/s)

N = length(wind_speed);
P_wtg = zeros(N, 1);

for i = 1:N
    v = wind_speed(i);

    % --- Yaw control: align with wind ---
    if nargin == 3
        theta_yaw = wind_dir(i);  % yaw error
    else
        theta_yaw = 0;
    end
    v_eff = v * cosd(theta_yaw);

    % --- Operating region check ---
    if v_eff < v_ci || v_eff >= v_co
        P_wtg(i) = 0;
        continue;
    end

    % --- Tip speed ratio ---
    lambda = (omega * R) / v_eff;

    % --- Available power at beta=0 (no pitching) ---
    Cp_max = compute_Cp(lambda, 0);
    P_avail = 0.5 * rho * A * Cp_max * v_eff^3 / 1000; % kW

    if P_avail <= P_rated
        % Below rated: run at max Cp
        P_wtg(i) = P_avail;
    else
        % Above rated: pitch to clamp at P_rated
        % Solve for beta such that Cp gives exactly P_rated
        beta_opt = fzero(@(b) ...
            0.5*rho*A*compute_Cp(lambda,b)*v_eff^3/1000 - P_rated, ...
            [0, 30]);
        P_wtg(i) = P_rated;
    end
end
end

function Cp = compute_Cp(lambda, beta)
    c = [0.5176, 116, 0.4, 5, 21, 0.0068];
    li = 1/(1/(lambda + 0.08*beta) - 0.035/(beta^3 + 1));
    Cp = c(1)*(c(2)/li - c(3)*beta - c(4))*exp(-c(5)/li) + c(6)*lambda;
    Cp = max(0, Cp); % Cp can't be negative
end