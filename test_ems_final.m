cd 'e:\PEDES work scheduling day ahead'
clear all
try
    ems_main
    disp('===SUCCESS===')
    disp('ems_main.m executed without errors')
    disp('===SUCCESS===')
catch ME
    disp('===ERROR===')
    disp('Error in ems_main.m:')
    disp(ME.message)
    disp('Line information:')
    for i = 1:length(ME.stack)
        fprintf('  File: %s, Line: %d, Function: %s\n', ME.stack(i).file, ME.stack(i).line, ME.stack(i).name)
    end
    disp('===ERROR===')
end
