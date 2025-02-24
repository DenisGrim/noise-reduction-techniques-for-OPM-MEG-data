function output_dirname = apply_filtering(p, input_dirname, output_dirname)
% Apply filtering
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set input and output directories
input_dir = fullfile(p.proj_root, input_dirname, p.task);
output_dir = fullfile(p.proj_root, output_dirname, p.task);

% EOG signals will be filtered
chname_extra = {'EOG1', 'EOG2'};

% In the case of the Motor task, EMG signals will be filtered as well
if strcmp(p.task, 'Motor')
    % Set reference channel
    chname_extra{end+1} = 'trigger_emg';
end

for run = 1:p.num_run
    % Set input and output files
    file_name = sprintf('run%02d', run);
    input_file = fullfile(input_dir, [file_name '.meg.mat']);
    output_file = fullfile(output_dir, [file_name '.meg.mat']);
    
    % Load OPM-MEG data (includes inactive channels because of down-sampling)
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    loadspec.ActiveTrial   = false; % Load all trial
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);
    freq = time_info.sample_frequency;
    
    % Load reference data
    loadspec = [];
    loadspec.ChannelType = 'REFERENCE';
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    loadspec.ActiveTrial   = false; % Load all trial
    [refmg, ~, ~] = vb_load_meg_data(input_file, loadspec);
    
    % Load extra-channel data (EOG and trigger)
    loadspec = [];
    loadspec.ChannelType = 'EXTRA';
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    loadspec.ActiveTrial   = false; % Load all trial
    [bexp_ext, ch_info_ext, ~] = vb_load_meg_data(input_file,loadspec);
    
    % Filter OPM-MEG data (only active channels are filtered)
    parm_filt = p.filter;
    parm_filt.chanNum = find(ch_info.Active); % Specify active channels
    bexp = filter2d(bexp, freq, parm_filt, 0);
    
    % Filter specified extra-channel data
    parm_filt_ex = p.filter_ex;
    ix_extra = [];
    if ~isempty(chname_extra)
        for ex = 1:length(chname_extra)
            ix_ex = find(strcmp(ch_info_ext.Name, chname_extra{ex}));
            if ~isempty(ix_ex)
                ix_extra = [ix_extra, ix_ex];
            else
                disp(['Extra channel ''' chname_extra{ex} ''' does not exist.'])
            end
        end
        parm_filt_ex.chanNum = ix_extra;
        bexp_ext = filter2d(bexp_ext, freq, parm_filt_ex, 0);
    end
    
    % Down sampling
    if isfield(parm_filt, 'fsamp') && freq ~= parm_filt.fsamp
        bexp     = vb_convert_freq(bexp, freq, parm_filt.fsamp);
        bexp_ext = vb_convert_freq(bexp_ext, freq, parm_filt.fsamp);
        refmg    = vb_convert_freq(refmg, freq, parm_filt.fsamp);
        new_freq = parm_filt.fsamp;
    else
        new_freq = freq;
    end
    
    % Overwrite input with updated data
    B = load(input_file); % To be overwritten
    B.bexp = bexp;
    B.bexp_ext = bexp_ext;
    B.refmg = refmg;
    B.MEGinfo.SampleFreq = new_freq;
    B.MEGinfo.Nsample = size(bexp, 2);
    
    % Save filtered data
    vb_fsave(output_file, '-struct', 'B');
    disp([output_file ' was saved.'])
end



