function d = define_dataset
% Define the basic properties of dataset
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

%% Define the subject and task
d.sub_list = {'002', '005', '006', '093'}; % 095 does not have OPM
d.task_list = {'Auditory', 'Motor', 'Somatosensory'};

%% Define the number of run
% 'Auditory', 'Motor', 'Somatosensory', 'Rest'
num_run_list_opm   = [2 3 2 1;  % 002
                   2 3 2 1;  % 005
                   2 3 2 1;  % 006
                   2 3 2 1;  % 093
                   0 0 0 0]; % 095

d.num_run_table_opm = array2table(num_run_list_opm,...
    'VariableNames', {'Auditory', 'Motor', 'Somatosensory', 'Rest'}, ...
    'RowNames', {'002', '005', '006', '093', '095'});

