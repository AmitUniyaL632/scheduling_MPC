cd 'e:\PEDES work scheduling day ahead'
clear all
% Test wind model function directly
config = struct('rho', 1.225, 'R', 1.5, 'omega', 12.57, 'P_rated', 2, 'v_cut_in', 3, 'v_cut_out', 25);
wind_speed_test = [5; 8; 10; 15; 20];
wind_dir_test = [0; 0; 0; 0; 0];
try
    P_out = wtg_physics_model(config, wind_speed_test, wind_dir_test);
    disp('SUCCESS: wtg_physics_model function works correctly')
    disp(P_out)
catch ME
    disp('ERROR:')
    disp(ME.message)
end
