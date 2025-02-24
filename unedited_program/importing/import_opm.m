function import_opm(p)
% Import OPM-MEG data
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Since OPM data is already imported as VBMEG format,
% just copy data from dataset.
data_orig = fullfile(p.raw_data_root, 'OPM', p.task);
data_copy = fullfile(p.proj_root, p.dirname.load, p.task);
if ~exist(data_copy, 'dir')
    mkdir(data_copy)
end

copyfile(data_orig, data_copy);
disp('OPM-MEG data in')
disp(data_orig)
disp('were copied in')
disp([data_copy '.'])
