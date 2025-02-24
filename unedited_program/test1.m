
% Disable toolboxes
license('checkout', 'Bioinformatics_Toolbox', 'disable')
license('checkout', 'Computer_Vision_Toolbox', 'disable')
license('checkout', 'Control_System_Toolbox', 'disable')
license('checkout', 'Curve_Fitting_Toolbox', 'disable')
license('checkout', 'Deep_Learning_Toolbox', 'disable')
license('checkout', 'Global_Optimization_Toolbox', 'disable')
license('checkout', 'Image_Processing_Toolbox', 'disable')
license('checkout', 'Optimization_Toolbox', 'disable')
license('checkout', 'Parallel_Computing_Toolbox', 'disable')
%license('checkout', 'Reinforcement_Learning_Toolbox', 'disable')
license('checkout', 'Robotics_System_Toolbox', 'disable')
license('checkout', 'Signal_Processing_Toolbox', 'disable')
license('checkout', 'Simulink_3D_Animation', 'disable')
license('checkout', 'Symbolic_Math_Toolbox', 'disable')
license('checkout', 'Identification_Toolbox', 'disable')
license('checkout', 'Wavelet_Toolbox', 'disable')
%license('checkout', 'Yokogawa_MEG_Reader_toolbox_for_MATLAB', 'disable')
license('checkout', 'Statistics_Toolbox', 'disable')


clear, close all

%% Preparation
% Add VBMEG to search path (modify this for your environment)
path_of_VBMEG = '/home/cbi/takeda/analysis/toolbox/vbmeg3_0_0_a_2';
addpath(path_of_VBMEG);
vbmeg

% Add all the directories in the current directory to search path
addpath(genpath(cd))

% Get properties of dataset
d = define_dataset;

% Set subject
for ss = 1:4
sub = d.sub_list{ss};

% Set task
for tt = 1:3
task = d.task_list{tt};
num_run = d.num_run_table_opm{sub, task};

% Set parameters for analyzing OPM-MEG data
p = set_parameters(sub, task, num_run);

% Add the information of bad channels to the parameters
p = set_bad_ch(p);

%% Processing OPM-MEG data
% Import OPM-MEG data
import_opm(p)

% Import EOG data in EEG data into OPM-MEG data
import_eog(p, p.dirname.load)

% Modify trigger signal
% In case of the Auditory task, trigger onsets will be shifted by 60 msec
if strcmp(p.task, 'Auditory')
    modify_trigger(p, p.dirname.load, 60);
else
    modify_trigger(p, p.dirname.load);
end

% Import the information of bad channels, which were manually selected
import_bad_ch(p, p.dirname.load);

% The order of following processes is easily changable by modifying input/ouput dirname
PREP = {p.dirname.load};

if strcmp(p.task, 'Auditory') || strcmp(p.task, 'Somatosensory')
    % Homogeneous field correction
    PREP{end+1} = apply_hfc(p, PREP{end}, p.dirname.hfc);
end

% Detrend data using spline interpolation
PREP{end+1} = apply_detrending(p, PREP{end}, p.dirname.detrend);

% Filter data
PREP{end+1} = apply_filtering(p, PREP{end}, p.dirname.filter);

% Check processed data
for pp = 1:length(PREP)
    show_processed_data(p, PREP{pp});
end
compare_processed_psd(p, PREP);

% Segment continuous data into trials (epochs)
make_trial(p, PREP{end});

% Reject bad channels and trials
prefix = []; % Initialize prefix
prefix = reject_chtr(p, prefix, PREP{end});

% Regress out EOG component from OPM-MEG data
prefix = regress_out_eog(p, prefix);

% Correct baseline
prefix = correct_baseline(p, prefix);

% Combine trials across runs
prefix = combine_trial(p, prefix); 

% Show trial-averaged OPM-MEG data
show_trial_average(p, prefix);


end
end
