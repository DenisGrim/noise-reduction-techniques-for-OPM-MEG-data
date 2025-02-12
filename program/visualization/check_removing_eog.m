function check_removing_eog(p, prefix_before, prefix_after, eog_ch)
% Check the result of removing EOG components
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set figure properties
close all
set_fig_property(4, 2, 15, 15);

% Set dir of files before and after EOG removing (same directory)
dir_trial = fullfile(p.proj_root, p.dirname.trial, p.task);

% Check prefix
if isempty(prefix_before)
    prefix_before_ = [];
else
    prefix_before_ = [prefix_before '_'];
end
prefix_after_ = [prefix_after '_'];

% Set EOG channels
% eog_ch = {'EOG1', 'EOG2'};

% Set modality-specific ch parameter
switch p.dirname.modality
    case 'opm-meg'
        extention = '.meg.mat';
        plot_type = [5, 6]; % y and z
        plot_label = {'y', 'z'}; % label of saved fig
    case 'squid-meg'
        extention = '.meg.mat';
        plot_type = 2; % axial channel
        plot_label = {'axial'}; % label of saved fig
    case 'eeg'
        extention = '.eeg.mat';
        plot_type = 1; % eeg channel
        plot_label = {'eeg'}; % label of saved fig
end

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    % Set input files
    file_before = fullfile(dir_trial, [prefix_before_ file_name extention]);
    file_after  = fullfile(dir_trial, [prefix_after_ file_name extention]);

    % Loop for ch type
    for cc = 1:length(plot_type)
        % Load data
        [bexp1, ~, time_info] = vb_load_meg_data(file_before);
        bexp2 = vb_load_meg_data(file_after);

        % Load positions of channels
        [pos, channel_info] = vb_load_channel(file_before);

        % Load EOG
        load_spec = [];
        load_spec.ChannelName = eog_ch;
        [eog, ~, time_info] = vb_load_meg_data(file_before, load_spec);

        % Normalize EOG
        eog_all = reshape(eog, [size(eog,1),size(eog,2)*size(eog,3)]); % flatten
        eog_all = bsxfun(@minus, eog_all, mean(eog_all,2));
        eog_all = bsxfun(@rdivide, eog_all, std(eog_all,[],2));
        eog = reshape(eog_all, [size(eog,1), size(eog,2), size(eog,3)]); % reshape

        % Extract data with specific ch type
        ch_type = plot_type(cc);
        ix_ch = find(channel_info.Type == ch_type);
        pos = pos(ix_ch, :);
        bexp1 = bexp1(ix_ch, :, :);
        bexp2 = bexp2(ix_ch, :, :);

        % Calculate power
        p1 = sum(sum(bexp1.^2, 2), 3);
        p2 = sum(sum(bexp2.^2, 2), 3);
        dp = p1-p2;

        % Select ch and trial indicating maximum difference
        [~, ch_md] = max(dp); % Select channel
        d = squeeze(bexp1(ch_md, :, :)-bexp2(ch_md, :, :));
        pd = sum(d.^2, 1);
        [~, tr_md] = max(pd); % Select trial

        %% Plot sensor space
        % Get minimum and maximum values of power
        ma = max(p1(:));
        mi = min(p2(:));

        % Show power of original MEG data
        close all
        h = figure;
        subplot(2, 2, 1)
        vb_plot_sensor_2d(pos, p1, [mi ma], true);
        axis square off
        colorbar
        title('Power (before)')
        colormap(my_parula())

        % Show power of EOG-removed MEG data
        subplot(2, 2, 2)
        vb_plot_sensor_2d(pos, p2, [mi ma], true);
        axis square off
        colorbar
        title('Power (after)')
        colormap(my_parula())

        % Show an example of original and EOG-removed MEG data
        subplot(2, 2, 3)
        plot(time_info.time, [bexp1(ch_md, :, tr_md)', bexp2(ch_md, :, tr_md)'])
        xlim([time_info.time(1) time_info.time(end)])
        grid on
        legend('Before', 'After', 'Location', 'northwest')
        xlabel('Time [s]')
        ylabel('Magnetic field [T]')
        title(['Channel ' channel_info.Name{ch_md} ', Trial ' num2str(tr_md)], 'interpreter', 'none')

        % Show location of selected channel
        subplot(2, 2, 4)
        x = zeros(length(ix_ch), 1);
        x(ch_md) = 1;
        vb_plot_sensor_2d(pos, x);
        vb_plot_sensor_2d_head_plot_add(gca);
        axis square off
        title(['Channel ' channel_info.Name{ch_md}], 'interpreter', 'none')
        colormap(my_parula())

        % Save figure
        fig_file = fullfile(p.fig_root, mfilename, p.dirname.modality, p.task, p.sub, ['space_' plot_label{cc} '_run' num2str(run)]);
        vb_savefig_as_shown(h, fig_file)
        disp([fig_file '.png was saved.'])
    end
end


