function prefix_out = combine_trial(p, prefix_in)
% Combine trials across runs
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% No additional prefix

% Set prefix of input file
if isempty(prefix_in)
    prefix_in_ = [];
else
    prefix_in_ = [prefix_in '_'];
end

% Set prefix of output file
prefix_out = prefix_in;
if isempty(prefix_out)
    prefix_out_ = [];
else
    prefix_out_ = [prefix_out '_'];
end

% Get paths of files to be combined
data_files = cell(p.num_run, 1);
for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    data_files{run} = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_in_ file_name '.meg.mat']);
end

% Set output file
combined_file = fullfile(p.proj_root, p.dirname.trial, [prefix_out_ p.task '.info.mat']);

% Combine trials
vb_make_fileinfo(data_files, combined_file);


