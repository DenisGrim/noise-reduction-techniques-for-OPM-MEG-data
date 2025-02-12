function modify_trigger(p, input_dirname, time_shift_msec)
% Modify trigger signal
% If 'time_shift_msec' is inputted, shift the trigger onsets by time_shift_msec.
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

if ~exist('time_shift_msec','var') || isempty(time_shift_msec)
    time_shift_msec = 0;
end

% Set input directory
input_dir = fullfile(p.proj_root, input_dirname, p.task);

% Set trigger channel
loadspec = [];
loadspec.ChannelType = 'EXTRA';
loadspec.ChannelName = {'trigger_pulse'}; % Original trigger name
% Modified trigger is saved as 'trigger_modified'
new_trig_ch_name = 'trigger_modified';
new_trig_ch_type = -1;% Registered as an extra channel (bexp_ext)

% Modify trigger signals so that they have time width of Nsample_trig.
Nsample_trig = 100;

% Make msec margine before the first trigger and after the last trigger
tmargine_first_sec = 0;
tmargine_last_sec  = 0;

% This is add-hoc operation for handling exception
% Because of the last triger of s093 Motor run01
if strcmp(p.sub, '093') && strcmp(p.task, 'Motor')
    tmargine_last_sec  = 3;
end

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    
    % Load trigger channel data from loaded dir
    meg_file = fullfile(input_dir, [file_name '.meg.mat']);    
    [trigger, ~, time_info] = vb_load_meg_data(meg_file, loadspec);

    % Convert time_shift_msec to time_shift_sample
    time_shift_sec = time_shift_msec / 1000;
    time_shift_sample = time_shift_sec * time_info.sample_frequency;

    % Binarize trigger
    trigger = trigger./max(trigger);
    trigger = trigger > 0.5;
    
    % Detect trigger onsets
    bindiff = trigger(1:end-1) - trigger(2:end);
    ix_onset = find(bindiff==-1)+1;

    % Shift trigger onsets
    ix_onset = ix_onset + time_shift_sample;

    % Check if the first trigger is out of bounds
    t0_sec = time_info.time(1);
    trg_time_1st_sec = time_info.time(ix_onset(1));
    pre_trg_sec  = p.Pretrigger_ms/1000;
    if (trg_time_1st_sec-pre_trg_sec-tmargine_first_sec) < t0_sec
        warning(['Triggers at ' num2str(trg_time_1st_sec) ' [sec] is removed because out of time.'])
        ix_onset(1) = [];
    end

    % Check if the last trigger is out of bounds
    te_sec = time_info.time(end);
    trg_time_end_sec = time_info.time(ix_onset(end));
    post_trg_sec = p.Posttrigger_ms/1000;
    if (trg_time_end_sec+post_trg_sec+tmargine_last_sec) > te_sec
        warning(['Triggers at ' num2str(trg_time_end_sec) ' [sec] is removed because out of time.'])
        ix_onset(end) = [];
    end

    % Make onset time-series
    trigger_onset = zeros(size(trigger));
    trigger_onset(ix_onset) = 1;
    
    % Modify trigger signal to square wave with a certain duration
    trigger_sq = zeros(size(trigger));
    for tt = 1:length(ix_onset)
        if ix_onset(tt)+Nsample_trig-1 > length(trigger)
            warning(['Triggers at ' num2str(ix_onset(tt)) ' is removed because out of time-length.'])
        else
            trigger_sq(ix_onset(tt):ix_onset(tt)+Nsample_trig-1) = 1;
        end
        
        if tt ~= length(ix_onset)
            if ix_onset(tt)+Nsample_trig-1 > ix_onset(tt+1)
                warning(['Triggers at ' num2str(ix_onset(tt)) ' and ' num2str(ix_onset(tt+1)) ' are merged into a trigger.'])
            end
        end
    end
    
    % Check results
    close all
    set_fig_property(2, 1, 15, 15)
    h = figure;
    plot(time_info.time, trigger, 'r'), hold on
    plot(time_info.time, trigger_sq, 'b')
    xlim([time_info.time(ix_onset(1))-0.5 time_info.time(ix_onset(1))+0.5])
    xlabel('Time [sec]')
    legend({'Original', 'Modified'})
    title(['Trigger signal around the first onset'])
    fig_file = fullfile(p.fig_root, mfilename, [p.sub '_' p.task '_' file_name]);
    vb_savefig_as_shown(h, fig_file)
    disp([fig_file '.png was saved.'])

    % Append modified trigger to input file
    append_ch = [];
    append_ch.data = trigger_sq;
    append_ch.name = new_trig_ch_name;
    append_ch.type = new_trig_ch_type;
    append_extra_ch(meg_file, append_ch)
    disp([meg_file ' was overwritten.'])
end

