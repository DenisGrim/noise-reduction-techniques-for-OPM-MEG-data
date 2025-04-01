function p = set_parameters(sub, task, num_run)
% Set parameters and directory names for analyzing OPM-MEG data

%% Set base directories and files
root_raw = fullfile(fileparts(cd), 'OPMdata_from_OSEdataset');
root_analyzed  = fullfile(fileparts(cd), 'analyzed_data');
p = struct;
p.sub = sub;
p.task = task;
p.raw_data_root = fullfile(root_raw, p.sub);% Raw data directory of a subject
p.proj_root = fullfile(root_analyzed, p.sub);% Directory to save analyzed data
p.fig_root = fullfile(root_analyzed, 'figure');% Directory to save figures


if isempty(num_run)
    p.num_run = 1;
else
    p.num_run = num_run;
end

%% set bandpass filter
parm_filter = struct(...
    'bandpass_freq', [6 50],...
    'bandpass_filtertype' , 'butter',...
    'bandpass_filterorder', 4,...
    'bandstop_freq', [57 63; 117 123],...
    'bandstop_filtertype' , 'butter',...
    'bandstop_filterorder', 4,...
    'fsamp', 500);
p.filter = parm_filter;
p.filter_ex = parm_filter;% For extra channels

%% For making trial data
p.Pretrigger_ms  = 300; % Before trigger [msec]
p.Posttrigger_ms = 500; % After trigger [msec]
p.snr_type = 'power'; % 'power' for power SNR. Defaults to average otherwise

%% Dynamic-range of OPM for detecting saturation used in vb_job_trial_onset_opm.m
p.dr_ub = 1.5*10^-9;  % Upper bound
p.dr_lb = -1*p.dr_ub; % Lower bound

%% Threshold for trial and channel rejection used in reject_chtr.m
% Lower/Upper, and standard deviation bounds
p.threshold_rejection = [1, 7, 2];

%% For correcting baseline
p.time_base_sec = [-p.Pretrigger_ms -200] ./ 1000 ; % time window of baseline [s]

%% For visualizing the trial-averaged OPM-MEG data
% Time of interests for detecting a peak
switch task
    case 'Somatosensory'
        toi = [0.01 0.05];
    case 'Auditory'
        toi = [0.05 0.15];
    case 'Motor'
        toi = [0.05 0.15];
end
p.time_of_interest_sec = toi;

% Time window to show the trial-averaged OPM-MEG data
p.time_show_sec = [-p.Pretrigger_ms p.Posttrigger_ms] ./ 1000;

%% File and directory names
parm_dirname = struct;
parm_dirname.modality = 'opm-meg';
parm_dirname.load = 'load';
parm_dirname.hfc = 'hfc';
parm_dirname.detrend = 'detrend';
parm_dirname.filter = 'filter';
parm_dirname.trial = 'trial';
p.dirname = parm_dirname;


