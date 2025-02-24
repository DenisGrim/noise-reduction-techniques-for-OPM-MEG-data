function show_processed_data(p, input_dirname, foi)
% Show processed results
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set figure properties
close all
set_fig_property(3, 2, 15, 15);

% Set freq of interests
if ~exist('foi', 'var') || isempty(foi)
    foi = [1; 100];
end

% Set parameters for OPM-MEG data
chtype = 'MEG';
unit_coef = 10^15; % T -> fT
unit_disp = 'fT';
noise_floor = 15; % Noise floor is 15 [fT/sqrt(HZ)]

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    data_file = fullfile(p.proj_root, input_dirname, p.task, [file_name '.meg.mat']);
    
    %% Load data
    loadspec = [];
    loadspec.ActiveChannel = true;% Load active channels only
    [y, channel_info, time_info] = vb_load_meg_data(data_file, loadspec);
    chnames = channel_info.Name;
    time = time_info.time;
    fs = time_info.sample_frequency;

    
    [nch, ntime, ntrial] = size(y);
    % If data is already epoched, make it flat
    if ntrial~=1
        y = reshape(y, [nch, ntime*ntrial]);
        dt = 1/fs;
        time = 0:dt:(ntime*ntrial-1)*dt; % Make long time
    end
    
    line_colors = my_lines(nch); % Get line colors
    line_colors = min(line_colors+0.0, 1); % Make color lighter with upper bound=1
    
    % Convert unit from T to fT
    y = y.*unit_coef;
    
    %% Plot time series
    h = figure;
    subplot(2, 1, 1), hold on
    max_y = max(y(:));
    min_y = min(y(:));
    for nn = 1:nch
        plot(time, y(nn, :), 'LineWidth', 0.1, 'Color', line_colors(nn,:))
    end
    plot(time, mean(y), 'LineWidth', 0.5, 'Color', 'k')
    xlim([time(1) time(end)])
    xlabel('Time [sec]')
    ylabel(['Magnetic field [' unit_disp ']'])
    title([input_dirname ' Run-' num2str(run)])
    grid on
    
    %% Plot power spectral density (PSD)
    % Set unit
    data_unit = unit_disp;  % Unit of original data
    spect_unit = 'psds'; % Unit of resulted PSD (see calc_psd)
    
    % Set window length [ms]
    len_win_ms = 2000;
    
    % Set window function
    % Format of window functon: func(window_length, options)
    func_win = @hanning;
    
    % Constant line draw in the figure
    constant = noise_floor;
    
    % Calculate PSD
    [freq, psdx] = calc_psd(y, fs, func_win, len_win_ms, spect_unit);
    
    % Set parameters for plotting PSD
    opt = [];
    opt.const = constant;
    opt.dunit = data_unit;
    opt.sunit = spect_unit;
    
    % % Plot PSD
    subplot(2, 1, 2),
    plot_psd(freq, psdx, opt)
    xlim(foi)
    set_axis
    
    %% Save figure
    fig_file = fullfile(p.fig_root, mfilename, p.task, [p.sub '_' input_dirname '_' file_name]);
    vb_savefig_as_shown(h, fig_file, '-dpng')
    disp([fig_file '.png was saved.'])
end
