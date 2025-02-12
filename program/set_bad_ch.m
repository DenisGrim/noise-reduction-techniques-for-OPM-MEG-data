function p = set_bad_ch(p)
% Append manually-selected bad channels to parameters.
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

%% Define the table for bad channels
badch_table = cell2table(cell(5,4),...
    'VariableNames', {'Auditory', 'Motor', 'Somatosensory', 'Rest'}, ...
    'RowNames', {'002', '005', '006', '093', '095'});

%% Set noisy channels to the table
% --- OPM ---
% 002
badch_table.Auditory{'002'}      = {'M005_y', 'M012_y', 'M012_z', 'M014_z', 'M015_z'};
badch_table.Motor{'002'}         = {'M005_y', 'M014_z', 'M015_z'};
badch_table.Somatosensory{'002'} = {'M005_y', 'M014_z', 'M015_z'};
% 005
badch_table.Auditory{'005'}      = {'M015_y', 'M015_z', 'M012_z'};
badch_table.Motor{'005'}         = {'M012_z', 'M013_y', 'M013_z'};
badch_table.Somatosensory{'005'} = {'M015_y', 'M015_z', 'M014_z'};
% 006
badch_table.Auditory{'006'}      = {'M014_y', 'M015_z'};
badch_table.Motor{'006'}         = {'M012_z', 'M015_z', 'M014_z'};
badch_table.Somatosensory{'006'} = {'M011_z', 'M012_z', 'M009_z'};
% 093
badch_table.Auditory{'093'}      = {'M005_y', 'M011_y', 'M012_z', 'M013_z', 'M015_z'};
badch_table.Motor{'093'}         = {'M014_z', 'M015_z'};
badch_table.Somatosensory{'093'} = {'M005_y', 'M015_y', 'M015_z', 'M014_z'};

%% Extract param for specified subject and task
badch = badch_table.(p.task){p.sub};

% Append to the parameter
p.badch = badch;