function output_dirname = apply_car(p, input_dirname, output_dirname, doPlot)
% Remove environmental noise using Homogeneous Field Correction (HFC) (Tierney et al., 2021)
%
% [Reference]
% Tierney et al. Modelling optically pumped magnetometer inteference in MEG
% as a spatially homogeneous magnetic field. NeuroImage 244 (2021) 118484.
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

    %copy data for 'before' plot later
    y = bexp;

    bexp = bexp - mean(bexp, 1);

    if doPlot
        % Show results of HFC
        set_fig_property(4, 3, 15, 15);
        close all
        h = figure; hold on
        subplot(3, 1, 1), plot(time, y'), title('Before Common Average')
        xlim([min(time) max(time)]), ylim([min(y(:)) max(y(:))])
        ylabel('Magnetic field [T]')
        subplot(3, 1, 2), plot(time, bexp'), title('After Common Average')
        xlim([min(time) max(time)]), ylim([min(bexp(:)) max(bexp(:))])
        ylabel('Magnetic field [T]')
        fig_file = fullfile(p.fig_root, mfilename, p.task, [p.sub '_' file_name '_car']);
        vb_savefig_as_shown(h, fig_file)
        disp([fig_file '.png was saved.'])
    end
%{
    %% Save projection matrix M
    file_proj = fullfile(output_dir, ['M_' file_name '.mat']);
    save(file_proj, 'M');
    disp(['The projection matrix M was saved in ' file_proj '.'])

    %% Save denoised data
    % If original bexp was epoched, re-shape it to original size
%}
    if ntrial~=1
        bexp = reshape(bexp, [nch, ntime, ntrial]);
        % Reshape estiamted homogeneous field as well
        %H    = reshape(H, [size(H,1), ntime, ntrial]);
    end

    %{%}
    % Make output file and append estimated HF as reference ch
    append_ch = [];
    %append_ch.data = H;
    append_ch.data = zeros(3, size(bexp, 2));
    append_ch.name = {'HFx', 'HFy', 'HFz'};
    append_ch.type = [257, 257, 257];
    append_extra_ch(input_file, append_ch, output_file); % It makes output_file

    % Overwrite bexp with updated data
    append_ch = [];
    append_ch.data = bexp;
    append_ch.name = ch_info.Name;
    append_ch.type = ch_info.Type;
    append_data_ch(output_file, append_ch)
    
%{
    B = load(input_file);
    B.bexp = bexp;
    vb_fsave(output_file, '-struct', 'B');
%}
    disp([output_file ' was overwritten.'])
end



