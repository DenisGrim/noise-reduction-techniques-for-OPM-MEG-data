function make_trial(p, input_dirname)
% Segment continuous data into trials using trigger signals
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set figure properties
close all
set_fig_property(4, 2, 15, 15);

% Set parameters
parm.status_ch = {'trigger_modified'};
parm.Pretrigger_ms  = p.Pretrigger_ms; % before trigger [msec]
parm.Posttrigger_ms = p.Posttrigger_ms; % after trigger [msec]
parm.condition = {'stim'};
parm.trig_type = 'analog';
parm.slope = 'low_to_high';
parm.status_level = 0.5;

plot_parm.mode  = 1; % = 1: plot all trial
plot_parm.NXmax = 25; % # of trial in X-axis
plot_parm.NYmax = 5; % # of subplot in Y-axis

extention = '.meg.mat';

% Set dynamic range of the OPM measurement system
parm.dr_ub = p.dr_ub;
parm.dr_lb = p.dr_lb;

% Make directory to save figures
fig_dir = fullfile(p.fig_root, mfilename, p.task, p.sub);
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    
    %% Detect trial onsets
    % Set input file (relative path from proj_root)
    parm.data_file = fullfile(input_dirname, p.task, [file_name extention]);
    
    % Set output trigger file (saved in the same dir as input)
    parm.trig_file = fullfile(input_dirname, p.task, [file_name '.trig.mat']);
    
    % Set raw data, which is  required for checking saturation
    parm.loaded_data_file = fullfile(p.dirname.load, p.task, [file_name extention]);
    
    % Detect trial onsets from trigger signals
    vb_job_trial_onset_opm(p.proj_root, parm);
    disp([fullfile(p.proj_root, parm.trig_file) ' was saved.'])
    
    % Plot and save trial info
    plot_parm.png = fullfile(fig_dir, ['run' num2str(run)]);
    vb_plot_status(fullfile(p.proj_root, parm.trig_file), plot_parm);
    
    %% Segment continuous data into trials using the trial onsets
    % Set input file (absolute path)
    input_file = fullfile(p.proj_root, input_dirname, p.task, [file_name extention]);
    proc_spec.trig_file = fullfile(p.proj_root, input_dirname, p.task, [file_name '.trig.mat']);
    
    % Set output file
    trial_file = fullfile(p.proj_root, p.dirname.trial, p.task, [file_name extention]);
    
    % Make trial data
    vb_msrmnt_make_trial_data(input_file, proc_spec, trial_file);
    disp([trial_file ' was saved.'])
end



