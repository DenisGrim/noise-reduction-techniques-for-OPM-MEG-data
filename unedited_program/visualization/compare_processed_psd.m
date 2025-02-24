function compare_processed_psd(p, input_dirnames)
% Show power spectral densities (PSDs) of processed results
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set figure properties
close all
set_fig_property(4, 2, 15, 15);

% Set parameters for OPM-MEG data
unit_coef = 10^15; % T -> fT
unit_disp = 'fT';
noise_floor = 15; % In case of OPM, noise floor is 15 [fT/sqrt(HZ)]

Nresult = length(input_dirnames);
line_colors = my_lines(Nresult); % Get line colors

for run = 1:p.num_run
    h = figure;
    
    for re = 1:Nresult
        input_dirname = input_dirnames{re};
        file_name = sprintf('run%02d', run);
        data_file = fullfile(p.proj_root, input_dirname, p.task, [file_name '.meg.mat']);
        
        %% Load data
        % Load OPM-MEG data
        loadspec = [];
        loadspec.ActiveChannel = true;% Load active channels only
        [y, channel_info, time_info] = vb_load_meg_data(data_file, loadspec);
        fs = time_info.sample_frequency;
        
        % Convert unit from T to fT
        y = y.*unit_coef;
        
        %% Plot power spectral density (PSD)
        % Calculate PSD
        spect_unit = 'psds'; % Unit of resulted PSD (see calc_psd)
        [freq, psdx] = calc_psd(y, fs, @hanning, 2000, spect_unit);
        
        % Average PSDs across channels
        psdx_mean = mean(psdx, 2);
        
        % Set parameters for plotting PSD
        opt = [];
        if re == Nresult
            opt.const = noise_floor;
        end
        opt.dunit = unit_disp;
        opt.sunit = spect_unit;
        opt.line_color = line_colors(re, :);
        opt.funit = 'Hz';
        
        % Plot PSD
        plot_psd(freq, psdx_mean, opt)
        xlim([0 70])
        set_axis
    end
    
    for re = 1:Nresult
        qw{re} = plot(nan, 'LineWidth', 2, 'Color', line_colors(re,:));
    end
    lgn = [input_dirnames 'noise floor'];
    qw{Nresult+1} = plot(nan, '--k', 'LineWidth',2);
    legend([qw{:}], lgn , 'location', 'northeastoutside')
    
    %% Save figure
    fig_file = fullfile(p.fig_root, mfilename, p.task, [p.sub '_' file_name]);
    vb_savefig_as_shown(h, fig_file, '-dpng')
    disp([fig_file '.png was saved.'])
end