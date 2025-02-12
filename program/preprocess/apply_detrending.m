function output_dirname = apply_detrending(p, input_dirname, output_dirname, interp_sec, doPlot)
% Remove temporal trends from data by spline interpolation
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set input directory
input_dir = fullfile(p.proj_root, input_dirname, p.task);

% Set output directory
output_dir = fullfile(p.proj_root, output_dirname, p.task);
if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

if ~exist('interp_sec','var') || isempty(interp_sec)
    % Default: thin out (resampling) every 2 sec for interpolation
    interp_sec = 2;
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
    
    % Load OPM-MEG data
    loadspec = [];
    loadspec.ChannelType = 'MEG';
    loadspec.ActiveChannel = false;% Load all channels
    loadspec.ActiveTrial   = false;% Load all trials
    [bexp, ch_info, time_info] = vb_load_meg_data(input_file, loadspec);
    freq = time_info.sample_frequency;

    % Detrend data
    y = bexp;
    time =[0:size(y,2)-1] ./ freq;
    y_det = detrend_spline(time, y, round(freq)*interp_sec);
    bexp = y_det;

    if doPlot
        % Show result of dentrending
        set_fig_property(4, 2, 15, 15);
        close all
        h = figure; hold on
        subplot(2,1,1), plot(time, y'), title('Before detrend')
        xlim([min(time) max(time)]), ylim([min(y(:)) max(y(:))])
        ylabel('Magnetic field [T]')
        subplot(2,1,2), plot(time, y_det'), title('After detrend')
        xlim([min(time) max(time)]), ylim([min(y_det(:)) max(y_det(:))])
        xlabel('Time [sec]'), ylabel('Magnetic field [T]')
        fig_file = fullfile(p.fig_root, mfilename, p.task, [p.sub '_' file_name]);
        vb_savefig_as_shown(h, fig_file)
        disp([fig_file '.png was saved.'])
    end

    % Make output file and overwrite bexp with updated data
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_data_ch(input_file, append_ch, output_file); % It makes output_file
    disp([output_file ' was overwritten.'])
end



