function import_eog(p, input_dirname)
% Import EOG data from EEG to OPM-MEG file
% EOG data is registered as 'EOG1' and 'EOG2' in extra-channel
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

set_fig_property(4, 2, 15, 15);

% Set input directory
input_dir = fullfile(p.proj_root, input_dirname, p.task);

% Set trigger or sync channel
switch p.task
    case {'Somatosensory', 'Auditory'}
        % Same trigger signal is contained in both MEG and EOG data.
        ext_chname = 'trigger_pulse';
    case {'Motor', 'Rest'}
        % Trigger signal is contained only in MEG data.
        % Instead of that, sync signal is shared.
        ext_chname = 'sync_pulse';
end

% EOG channels will be registered as
eog_ch_name = {'EOG1', 'EOG2'};
eog_ch_type = [-2, -2];

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    meg_file = fullfile(input_dir, [file_name '.meg.mat']);
    eog_file = fullfile(input_dir, [file_name '.eeg.mat']);
    
    %% Load data
    % Load trigger or sync signal
    loadspec = [];
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    loadspec.ChannelType = 'EXTRA';
    loadspec.ChannelName = {ext_chname};
    [sync_meg, channel_info, time_info] = vb_load_meg_data(meg_file, loadspec);
    time = time_info.time;
    freq = time_info.sample_frequency;
    
    % Load EOG data
    loadspec = [];
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    [eog, channel_info_eog, time_info_eog] = vb_load_meg_data(eog_file, loadspec);
    time_eog = time_info_eog.time;
    freq_eog = time_info_eog.sample_frequency;

    % Load trigger or sync signal
    loadspec = [];
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    loadspec.ChannelType = 'EXTRA';
    loadspec.ChannelName = {'EXT1'}; % EXT1 is either trigger or sync
    [sync_eeg, ~, ~] = vb_load_meg_data(eog_file, loadspec);

    % Up-sampling EOG data and trigger(sync) signal from 1000 to 2000 Hz
    eog_up = vb_convert_freq(eog, freq_eog, freq);
    sync_eeg_up = vb_convert_freq(sync_eeg, freq_eog, freq);

    %% Import 
    [eog_new, sync_onset_new, sync_onset_src, sync_onset_dst] = import_signal_using_sync(sync_eeg_up, sync_meg, eog_up);

    % Compare onset index of dst and new (imported src index)
    ix_onset_new = find(sync_onset_new);
    ix_onset_dst = find(sync_onset_dst);
    dif = abs(ix_onset_dst-ix_onset_new);

    dif_ms = dif ./ freq*1000;
    max_diff_ms = num2str(max(dif)/freq*1000);
    
    %% Check result
    h = figure;
    subplot(2,1,1), plot(time, sync_onset_dst - sync_onset_new)
    xlim([-inf inf])
    xlabel('Time [sec]')
    ylabel('Differences [-1/0/1]')
    title('Sync onset (OPM) - Sync onset (imported EOG)')
    subplot(2,1,2), plot(dif_ms), 
    xlim([-inf inf])
    xlabel('#Pulse')
    ylabel('Differences [ms]')
    title(['Diff of pulse timing between EOG and OPM (max ' num2str(max_diff_ms) ' ms)'])
    fig_file = fullfile(p.fig_root, mfilename, [p.sub '_' p.task '_' file_name]);
    vb_savefig_as_shown(h, fig_file)
    disp([fig_file '.png was saved.'])

    % Append EOG ch to input file
    append_ch = [];
    append_ch.data = eog_new;
    append_ch.name = eog_ch_name;
    append_ch.type = eog_ch_type;
    append_extra_ch(meg_file, append_ch)
    disp([meg_file ' was overwritten.'])
end



