function output_dirname = apply_harmonicfc(p, input_dirname, output_dirname, doPlot)

disp(mfilename);

% Set input directory
input_dir = fullfile(p.proj_root, input_dirname, p.task);

% Set output directory
output_dir = fullfile(p.proj_root, output_dirname, p.task);
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

% If true, save the figure of results
% Defalt: true
if ~exist('doPlot', 'var') || isempty(doPlot)
    doPlot = true;
end

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    input_file = fullfile(input_dir, [file_name '.meg.mat']);
    output_file = fullfile(output_dir, [file_name '.meg.mat']);

    % Load MEG data
    % (bad channels are ignored because they may affect estimated homogeneous field)
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = true; % Load active channels
    loadspec.ActiveTrial   = false; % Load all the trials
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);
    time = time_info.time;
    freq = time_info.sample_frequency;
    [nch, ntime, ntrial] = size(bexp);
    % If bexp is already epoched, make it flat
    if ntrial~=1
        bexp = reshape(bexp, [nch, ntime*ntrial]);
    end

    % quick data copy for plot
    y = bexp;

    %% harmonic field correction??
    % Load Qpick of active channels
    [pick, Qpick] = vb_load_sensor(input_file, 'MEG');

    args=[];
    args.o= Qpick; % sensor orientation
    args.v = pick; % sensor position
    args.li = 1; % harmonic order
    % use spm to get X
    X = den_spm_opm_vslm(args);

    % get projector matrix
    M = eye(size(X,1))-X*pinv(X);

    % Normalize data
    [meanNormalized, ave1, std1] = normalize_data(y, 'ch_mean');
    [deviationNormalizedData, ave2, std2] = normalize_data(meanNormalized, 'variance');

    % Apply HF
    correctedDataNormalized = M * deviationNormalizedData;
    % Hz = Xi * deviationNormalized; % calculated and not used in reference
    HarmonicField = pinv(X)* meanNormalized; %  not sure if this is correct?

    % Re-normalize data
    correctedDataMeanNormalized = re_normalize_data(correctedDataNormalized, 'variance', ave1, std1);
    correctedData = re_normalize_data(correctedDataMeanNormalized, 'ch_mean', ave2, std2);

    bexp = correctedData;


    if doPlot
        % Show results of HFC
        set_fig_property(4, 3, 15, 15);
        close all
        h = figure; hold on
        subplot(3, 1, 1), plot(time, y'), title('Before Harmonic FC') %title('Before HFC')
        xlim([min(time) max(time)]), ylim([min(y(:)) max(y(:))])
        ylabel('Magnetic field [T]')
        
            subplot(3, 1, 2), plot(time, correctedData'), title('After Harmonic FC')
        xlim([min(time) max(time)]), ylim([min(correctedData(:)) max(correctedData(:))])
        ylabel('Magnetic field [T]')
        subplot(3, 1, 3), plot(time, HarmonicField'), title('Estimated harmonic magnetic field')
        xlim([min(time) max(time)]), ylim([min(HarmonicField(:)) max(HarmonicField(:))])
        xlabel('Time [sec]'), ylabel('Magnetic field [T]'), legend({'x' 'y' 'z'})
        fig_file = fullfile(p.fig_root, mfilename, p.task, [p.sub '_' file_name '_hf']);
        vb_savefig_as_shown(h, fig_file)
        disp([fig_file '.png was saved.'])
        
    end

    %% Save projection matrix M
    
    file_proj = fullfile(output_dir, ['M_' file_name '.mat']);
    save(file_proj, 'M');
    disp(['The projection matrix M was saved in ' file_proj '.'])
    

    %% Save denoised data
    % If original bexp was epoched, re-shape it to original size
    if ntrial~=1
        bexp = reshape(bexp, [nch, ntime, ntrial]);
        % Reshape estiamted harmonic field as well
        HarmonicField = reshape(HarmonicField, [size(HarmonicField,1), ntime, ntrial]);
    end

    % Make output file and append estimated HF as reference ch
    
    append_ch = [];
    append_ch.data = HarmonicField;
    append_ch.name = {'HFx', 'HFy', 'HFz'};
    append_ch.type = [257, 257, 257];
    append_extra_ch(input_file, append_ch, output_file); % It makes output_file
    

    % Overwrite bexp with updated data
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_data_ch(output_file, append_ch)

    disp([output_file ' was overwritten.'])
end



