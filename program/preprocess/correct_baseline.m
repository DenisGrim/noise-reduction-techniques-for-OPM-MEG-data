function prefix_out = correct_baseline(p, prefix_in)
% Correct baseline of MEG data
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Define additional prefix of this process
prefix_add = 'b';

% Set prefix of input file
if isempty(prefix_in)
    prefix_in_ = [];
else
    prefix_in_ = [prefix_in '_'];
end

% Set prefix of output file
prefix_out = [prefix_add prefix_in];
prefix_out_ = [prefix_out '_'];

% Set parameters
twin_base_sec = p.time_base_sec; % time window of baseline

% Load time information
file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_in_ 'run01.meg.mat']);
[~, ~, time_info] = vb_load_meg_data(file);
[~, from] = min(abs(time_info.time-twin_base_sec(1)));
[~, to] = min(abs(time_info.time-twin_base_sec(2)));

for run = 1:p.num_run    
    file_name = sprintf('run%02d', run);
    input_file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_in_ file_name '.meg.mat']);
    output_file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_out_ file_name '.meg.mat']);

    % Load MEG/EEG data (only active ch and trial)
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = true; % Load active channel
    loadspec.ActiveTrial   = true; % Load active trial
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);
    
    % Load active trial info
    info = vb_load_measurement_info(input_file, 'MEGINFO');

    % Make average of pre-stimulus  MEG to 0
    for ch = 1:size(bexp, 1)
        for tr = 1:size(bexp, 3)
            bexp_mean = mean(bexp(ch, from:to, tr), 2);
            bexp(ch, :, tr) = bexp(ch, :, tr)-bexp_mean;
        end
    end
    
    % Save baseline-corrected data
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_ch.trial = info.ActiveTrial;
    append_data_ch(input_file, append_ch, output_file);
    disp([output_file ' was overwritten.'])
end

