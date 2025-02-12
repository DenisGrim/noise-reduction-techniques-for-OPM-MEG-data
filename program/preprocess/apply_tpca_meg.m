function output_dirname = apply_tpca_meg(p, input_dirname, output_dirname)
% Apply temporal PCA and remove environmental noise from MEG data using reference sonsor data.
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set input directory
input_dir = fullfile(p.proj_root, p.dirname.modality, input_dirname, p.task);

% Set output directory
output_dir = fullfile(p.proj_root, p.dirname.modality, output_dirname, p.task);
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

% Set parameters for temporal PCA
parm_denoise = p.tpca;

for run = 1:p.num_run
    file_name = [sprintf('run%02d', run), '.meg.mat'];
    input_file = fullfile(input_dir, file_name);
    output_file = fullfile(output_dir, file_name);
    
    % Load OPM-MEG data
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = false; % Load all ch
    loadspec.ActiveTrial   = false; % Load all trial
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);
    time = time_info.time;
    freq = time_info.sample_frequency;
    
    % Load reference
    loadspec = [];
    loadspec.ChannelType = 'REFERENCE';
    loadspec.ActiveChannel = false; % Load all ch
    loadspec.ActiveTrial   = false; % Load all trial
    [refmg, ch_info_ref, ~] = vb_load_meg_data(input_file,loadspec);

    % Bias correction before tPCA
    bexp = bsxfun(@minus, bexp, mean(bexp,2));
    refmg = bsxfun(@minus, refmg, mean(refmg,2));
    
    % Remove noise using tPCA
    parm_denoise.sampling_freq = freq;
    [bexp, ~, ~] = remove_by_temporalPCA(bexp, refmg, parm_denoise);

    % Save denoised MEG data to output file
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_data_ch(input_file, append_ch, output_file)

    % Save corrected ref data to output file
    append_ch = [];
    append_ch.data = refmg;
    append_ch.name = ch_info_ref.Name;
    append_ch.type = ch_info_ref.Type;
    append_extra_ch(output_file, append_ch)
end


