clear, close all

%% Preparation
% Add VBMEG to search path (modify this for your environment)
path_of_VBMEG = '/home/denis/NAIL/vbmeg3_0_0_a_2';
addpath(path_of_VBMEG);
vbmeg

% Add all the directories in the current directory to search path
addpath(genpath(cd))

% Get properties of dataset
dataset = define_dataset;

% Set subject, task, and the number of runs
ss = 1;% 1, 2, 3, or 4
sub = dataset.sub_list{ss};
tt = 1;% 1: Auditory, 2: Motor, 3: Somatosensory
task = dataset.task_list{tt};
num_run = dataset.num_run_table_opm{sub, task};

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

PREP = {p.dirname.load};

% noise correction
if strcmp(p.dirname.noise_technique, 'homfc')
        PREP{end+1} = apply_hfc(p, PREP{end}, p.dirname.noise_technique);
elseif strcmp(p.dirname.noise_technique, 'harfc')
        PREP{end+1} = apply_harmonicfc(p, PREP{end}, p.dirname.noise_technique);
elseif strcmp(p.dirname.noise_technique, 'car')
        PREP{end+1} = apply_car(p, PREP{end}, p.dirname.noise_technique);
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

% Reject bad channels and trials based on data statistics
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

