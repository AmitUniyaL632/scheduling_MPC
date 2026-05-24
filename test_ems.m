cd 'e:\PEDES work scheduling day ahead'
clear all
try
    ems_main
    disp('SUCCESS: ems_main.m executed without errors')
catch ME
    disp('ERROR in ems_main.m:')
    disp(ME.message)
    disp('Stack trace:')
    disp(ME.stack)
end
