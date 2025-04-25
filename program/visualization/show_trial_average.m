function  show_trial_average(p, prefix_in)
% Show trial-averaged OPM-MEG data

disp(mfilename);

% No additional prefix

% Set prefix of input file
if isempty(prefix_in)
    prefix_in_ = [];
else
    prefix_in_ = [prefix_in '_'];
end

% Set figure properties
close all
set_fig_property(3, 2, 15, 15);

% Set parameters
%twin_to_show_sec = p.time_show_sec; % sec
%twin_of_interest_sec = p.time_of_interest_sec; % sec


%%
data_file = fullfile(p.proj_root, p.dirname.trial, ['ber_' p.task '.info.mat']);
% Load OPM-MEG data, channel, and time information
[data, channel_info, time_info] = vb_load_meg_data(data_file);
time = time_info.time;

[~, from_toi] = min(abs(time-p.time_of_interest_sec(1)));
[~, to_toi] = min(abs(time-p.time_of_interest_sec(2)));


%% calculate and write SNR into file
get_SNR(p, data, time)
%%
% Average OPM-MEG data across trials
mdata = mean(data, 3);

% Detect a peak within p.time_of_interest_sec
power = sum(mdata.^2, 1);
[~, t_peak] = max(power(from_toi:to_toi));
t_peak = t_peak+from_toi-1;

% Load positions of channels
pos = vb_load_channel(data_file);

% Show trial-averaged OPM-MEG data
h = figure;

% Time-series
subplot(2, 1, 1)
plot(time, mdata')
hold on
ma = max(abs(mdata(:)));
plot([time(t_peak) time(t_peak)], [-ma ma], 'k')
axis([p.time_show_sec(1) p.time_show_sec(2) -ma ma])
grid on
xlabel('Time [sec]')
ylabel('Magnetic field [T]')
title(['Trial-averaged OPM-MEG data (' p.task ')'])

% Spatial pattern of Z-axis at the peak
subplot(2, 2, 3)
ix_ch = find(channel_info.Type == 6);
ma = max(abs(mdata(ix_ch, t_peak)));
vb_plot_sensor_2d(pos(ix_ch, :), mdata(ix_ch, t_peak), [-ma ma], true);
axis equal off
colorbar
title(['Z-axis (' num2str(time(t_peak)) ' sec)'])

% Spatial pattern of Y-axis at the peak
subplot(2, 2, 4)
ix_ch = find(channel_info.Type == 5);
ma = max(abs(mdata(ix_ch, t_peak)));
vb_plot_sensor_2d(pos(ix_ch, :), mdata(ix_ch, t_peak), [-ma ma], true);
axis equal off
colorbar
title(['Y-axis (' num2str(time(t_peak)) ' sec)'])
set_axis

% Save figure
fig_file = fullfile(p.fig_root, mfilename, p.task, p.sub);
vb_savefig_as_shown(h, fig_file)
disp([fig_file '.png was saved.'])

