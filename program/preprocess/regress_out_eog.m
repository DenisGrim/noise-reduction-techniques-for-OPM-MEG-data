function prefix_out = regress_out_eog(p, prefix_in)
% Regress out EOG components from epoched data.
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Define additional prefix of this process
prefix_add = 'e';

% Set prefix of input file
if isempty(prefix_in)
    prefix_in_ = [];
else
    prefix_in_ = [prefix_in '_'];
end

% Set prefix of output file
prefix_out = [prefix_add prefix_in];
prefix_out_ = [prefix_out '_'];

% Set EOG channels
eog_ch = {'EOG1', 'EOG2'};

% Set figure properties
close all
set_fig_property(6, 3, 10, 10);

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    input_file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_in_ file_name '.meg.mat']);
    output_file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_out_ file_name '.meg.mat']);

    % Load OPM-MEG data (only active channels and trials)
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = true; % Load active channel
    loadspec.ActiveTrial   = true; % Load active trial
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);

    % Load active trial info
    info = vb_load_measurement_info(input_file, 'MEGINFO');

    % Load EOG data
    eog_file = input_file;
    load_spec = [];
    load_spec.ChannelName = eog_ch;
    [eog, ~, time_info] = vb_load_meg_data(eog_file, load_spec);
    
    % Flatten data
    eog_all = reshape(eog, [size(eog,1),size(eog,2)*size(eog,3)]);
    bexp_all = reshape(bexp, [size(bexp,1),size(bexp,2)*size(bexp,3)]);

    % Normalize explanatory variables
    eog_all = bsxfun(@minus, eog_all, mean(eog_all,2));
    eog_all = bsxfun(@rdivide, eog_all, std(eog_all,[],2));

    % Reshape data
    eog = reshape(eog_all, [size(eog,1), size(eog,2), size(eog,3)]);

    % Calculate weight, where bexp = weight x eog
    w = zeros(size(bexp, 1), size(eog_all, 1)+1);
    for ch = 1:size(bexp, 1)
        w(ch, :) = vb_multiple_regression(eog_all, bexp_all(ch, :));
    end
    
    % Remove EOG components from bexp
    for tr = 1:size(bexp, 3)
        for ch = 1:size(bexp, 1)
            pmeg = vb_pre_multiple_regression(eog(:,:,tr), w(ch, :));
            bexp(ch, :, tr) = bexp(ch, :, tr)-pmeg;
        end
    end
    
    % Save EOG-free data
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_ch.trial = info.ActiveTrial;
    append_data_ch(input_file, append_ch, output_file); % It makes output_file
    disp([output_file ' was overwritten.'])
end

% Plot before and after regress out
check_removing_eog(p, prefix_in, prefix_out, eog_ch)

