function	vb_job_trial_onset_opm(proj_root,parm)
% get time index for each trial by checking status channel
% --- Usage
%    vb_job_trial_onset(proj_root,parm)
% --- Input
%  parm.data_file : Data file name         [string]
%  parm.trig_file : Trial onset file name  [string]
%  parm.status_ch : status channel name    [string]
%  parm.trig_type = 'integer' or 'analog' or 'voice' or 'emg'
%  parm.slope = 'const_start' or 'const_end'    if trig_type = 'integer'
%               'low_to_high' or 'high_to_low'  if trig_type = 'analog'
%               No meaning                      if trig_type = 'voice','emg'
%  parm.condition : string describing condition [string or cell array]
%  parm.status_level : status level      [1 x Ncomdition]
%  parm.Pretrigger_ms : Pretrigger period   [msec]
%  parm.Posttrigger_ms: Posttrigger period  [msec]
% --- Save variables
% status : status signal
% status_out = (smoothed) status signal used for onset search
% status_val(m) = status value for m-th condition (m=1:Ncomdition)
% trig(n)       : Onset time index for n-th trial 
% cond_id(n)    : Condition ID for n-th trial 
% ix_trial(:,n) : Time index for n-th trial   [Tperiod x Ntrial]
%                 Tperiod : # of time sample in one trial
%                 Ntrial  : # of trials
% parm : parameter setting
% parm.fsamp : Sample Frequency [Hz]
%
% --- Optional parameter for EMG onset
% parm.t_event  : minimum distance from previous onset event [150 ms]
%                 distance from previous onset should be larger than t_event
% parm.p_val : P-value corresponding to the threshold for [EMG, smooth(EMG)] 
%	           [0.0001, 0.0005] or [0.0001, 0.001] 
%     cumulative histgram is used to determine threshold from P-value
%
% --- Usually following parameters need not be changed
% parm.hist_mode : histgram mode [1]
%                = 1: Estimate threshold by gamma distribution approximation
%                = 0: Estimate threshold by histgram
% parm.t_smooth : moving average window length               [25 ms]
% parm.t_slope  : slope estimation period near onset         [25 ms]
%                 if t_slope==0, zero cross point is not calculated
% parm.t_peak   : peak evaluation period                     [100 ms]
% peak_val : EMG value should exceed peak_val within 't_peak' after onset
%          = mean(peak value > threshold) if hist_mode = 1
%          = (max(y) * status_level)      if hist_mode = 0 or 2
% --- Condition for EMG onset
% 1. distance from previous onset should be larger than t_event
% 2. distance between EMG & smoothed EMG onset should be smaller than t_event
% 3. EMG value should exceed peak_val within t_peak after onset
% 4. zero cross point is estimated 
%    by linear fitting around smoothed EMG threshold point
%
%
% 2009-6-14  Masa-aki Sato
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

% 2021-10-26 K.Suzuki 
%            Add saturation detection for opm data
% 2023-05-10 K.Suzuki
%            Add channel jump detection

proj_root = vb_rm_trailing_slash(proj_root);

% Original data file
fname = [proj_root filesep parm.data_file ];
fname_loaded = [proj_root filesep parm.loaded_data_file ];
ftrig = [proj_root filesep parm.trig_file ];

%% Pre-screening1
% Check saturated time-points of unprocessed data 
% Here, rejected channels are ignored
% Load unprocessed MEG data (includes inactive channels and trials)
loadspec = [];
loadspec.ChannelType = 'MEG';
loadspec.ActiveChannel = false; % Load all ch
loadspec.ActiveTrial   = false; % Load all trial
[y, channel_info, time_info] = vb_load_meg_data(fname_loaded, loadspec);
chnames = channel_info.Name;
Nchan = length(chnames);
time = time_info.time;
fs = time_info.sample_frequency; % fsamp before down sampling

% Mark saturated time-period that get over bound of dynamic-range (for each channel)
flag_sat = or(y>parm.dr_ub, y<parm.dr_lb);

% Plot results
flag_sat_plot = flag_sat .* repmat([1:Nchan]', [1,length(time)]);

max_y = Nchan+0.5;
min_y = -0.5;

h = figure;
plot(time, flag_sat_plot, 'LineWidth',1)
xlim([time(1) time(end)])
ylim([min_y max_y])
xlabel('Time [s]')
ylabel('Saturated channel')
% legend(chnames, 'Location','northeastoutside')
title(['Check saturation'])
grid on

% Extract saturated timepoint
t_sat = flag_sat .* repmat(time, [Nchan,1]);

%% Pre-screening2
% Check jumping time-points of processed data 
% Here, rejected channels are ignored
% Load unprocessed MEG data (includes inactive channels and trials)
loadspec = [];
loadspec.ChannelType = 'MEG';
loadspec.ActiveChannel = false; % Load all ch
loadspec.ActiveTrial   = false; % Load all trial
[y, channel_info, time_info] = vb_load_meg_data(fname, loadspec);
chnames = channel_info.Name;
Nchan = length(chnames);
time = time_info.time;
fs = time_info.sample_frequency; % fsamp before down sampling

% In case of applying filter before jumping detection, un-comment below
% freq = time_info.sample_frequency;
% parm_filt = struct(...
%     'bandpass_freq', [6 50],...
%     'bandpass_filtertype' , 'butter',...
%     'bandpass_filterorder', 4,...
%     'bandstop_freq', [57 63; 117 123],...
%     'bandstop_filtertype' , 'butter',...
%     'bandstop_filterorder', 4);
% y = subFilter2d(y, freq, parm_filt, 0);
   
% Calc temporal difference for each ch
y_diff = diff(y, 1, 2);
% Threshold for each ch is defined as
% median + n_sd*SD
% If detection is not good, modify n_sd
n_sd = 15;
y_diff_median = median(y_diff, 2);
y_diff_std = std(y_diff, 0, 2);
thresh_ch = y_diff_median + y_diff_std * n_sd; 
% Take mean of ch-wise threshold
thresh = mean(thresh_ch);
% Detect outlier
y_outlier = abs(y_diff) > thresh;
% Compensate timepoint reduced by temporal difference 
y_outlier = [zeros(size(y_outlier,1),1) y_outlier];
% Convert to time
t_jump = y_outlier .* repmat(time, [Nchan,1]);

% For adjusting n_sd
figure, 
subplot(2,1,1)
plot(time_info.time(2:end), abs(y_diff)'), title('Diff(bexp) and threshold')
hold on
plot(time_info.time(2:end), repmat(thresh, 1,size(y_diff,2)))
subplot(2,1,2)
plot(time_info.time(2:end), y_outlier(:,2:end)'), title('Detected jumping')

%% Detect trial onset
% Load info
info = vb_load_meg_info(fname);
% Sample Frequency [Hz]
parm.fsamp  = info.SampleFreq;

if iscell(parm.status_ch),
	status_ch = parm.status_ch;
else
	status_ch = {parm.status_ch};
end

% Load status channel
loadspec = [];
loadspec.ChannelName = status_ch;
status = vb_load_meg_data(fname, loadspec);

[ix_trial, trig, cond_id, status_val,status_out] = ...
		vb_get_trial_time_index(status,parm);

if isempty(ix_trial), return; end;

% ix_trial(:,n) : Time index for n-th trial   [Tperiod x Ntrial]
tmin = min(ix_trial,[],1);
tmax = max(ix_trial,[],1);

%% Check violations
% Check time is inside the data
ix1 = find( (tmin > 0) & (tmax <= info.Nsample));
% Check trial period is enough
tperiod = tmax-tmin;
tlength = size(ix_trial,1)-1;
ix2 = find(tperiod==tlength); % tperiod must equal to tlength
% Merge indices
ix = intersect(ix1, ix2);
% Reject violated trials
trig     = trig(ix);
ix_trial = ix_trial(:,ix);
cond_id  = cond_id(ix);

% Check time overlapping between trials
% This should be done for after removing violated trials to save as many trials as possible

flag_ol = inner_reject_overlap_trial(ix_trial); % if overlapping, merge latar trial to former one
% flag_ol = inner_check_overlap(ix_trial); 

ix = find(~flag_ol);
% Reject overlapping trials
trig     = trig(ix);
ix_trial = ix_trial(:,ix);
cond_id  = cond_id(ix);

%% Detect saturated trial
[~, ~, tinfo] = vb_load_meg_data(fname);

ix_start = ix_trial(1,:);
ix_end = ix_trial(end,:);
% Index to time
t_start = tinfo.time(ix_start);
t_end = tinfo.time(ix_end);

time_sat = cell(Nchan,1);
trial_sat = cell(Nchan,1);
for cc=1:Nchan
    ts = setdiff(unique(t_sat(cc,:)), 0); % Extract saturation time
    % Find sat time occuring in trial
    % This is not good way (but comparing index is disable due to down-sampling)
    tr = arrayfun(@(x) find(and(t_start<x, t_end>x)), ts, 'UniformOutput',false);
    trial_sat{cc} = unique(cell2mat(tr)); 
    time_sat{cc} = ts;
end

%% Detect jumping trial
[~, ~, tinfo] = vb_load_meg_data(fname);

ix_start = ix_trial(1,:);
ix_end = ix_trial(end,:);
% Index to time
t_start = tinfo.time(ix_start);
t_end = tinfo.time(ix_end);

time_jump = cell(Nchan,1);
trial_jump = cell(Nchan,1);
for cc=1:Nchan
    ts = setdiff(unique(t_jump(cc,:)), 0); % Extract jumping time
    % Find jump time occuring in trial
    tr = arrayfun(@(x) find(and(t_start<x, t_end>x)), ts, 'UniformOutput',false);
    trial_jump{cc} = unique(cell2mat(tr)); 
    time_jump{cc} = ts;
end

%% Save trial file
vb_fsave(ftrig,'ix_trial','trig','status','status_out','status_val','cond_id','parm', 'time_sat', 'trial_sat', 'time_jump', 'trial_jump');

end % End of function

% This flagged later one (skip former)... not good
function flag_overlap = inner_check_overlap(ix_trial)

Ntrial = size(ix_trial,2);
% Init flag of overlapping trial
flag_overlap = zeros(1,Ntrial);
% Checking loop
tt = 1;
while tt<Ntrial
    % Check overlap btw tt and tt+1
    if ~isempty(intersect(ix_trial(:,tt), ix_trial(:,tt+1)))
        flag_overlap(tt+1) = 1; % Mark overlapping trial
        tt = tt+1; % Skip checking for overlapping trial
    end
    tt = tt+1;
end

end

% Seach overlapped trials
function ix_overlap = inner_search_overlap(ix_trial)
ix_end_former   = ix_trial(end, 1:end-1);
ix_start_later  = ix_trial(1, 2:end);

% Overlap check
ix_overlap = find( (ix_end_former - ix_start_later)>=0 );

% Skip later one
if ~isempty(ix_overlap)
    ix_overlap = ix_overlap+1;
end

end

% This function searches overlapped trials from a head trial
% and reject later trial.
% caseA. If trial 2&3, 3&4 are overlapped, then
%           step1. found overlap 2&3
%           step2. reject 3
%           step3. found no overlap, finish
% caseB. If trial 2&3, 4&5 are overlapped, then
%           step1. found overlap 2&3
%           step2. reject 3
%           step3. found overlap 4&5
%           step4. reject 5
%           step5. found no overlap, finish
function flag_overlap = inner_reject_overlap_trial(ix_trial)
Ntrial = size(ix_trial,2);
% Init flag of overlapping trial
flag_overlap = zeros(1,Ntrial);
ix_ol = -1; % init

while ~isempty(ix_ol)
% Seach overlapped trial
ix_ol = inner_search_overlap(ix_trial(:,~flag_overlap));
if ~isempty(ix_ol)
    % Shift with num of overlap
    flag_overlap(ix_ol(1)+sum(flag_overlap)) = 1; 
end

end

end