function import_bad_ch(p, input_dirname)
% Import the information of bad channels that were manually selected
% p.badch is specified by either index or channel name
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set input directory
input_dir = fullfile(p.proj_root, input_dirname, p.task);

for run = 1:p.num_run
    file_name = [sprintf('run%02d', run), '.meg.mat'];
    input_file = fullfile(input_dir, file_name);
    
    % Load MEG data
    loadspec = [];
    loadspec.ActiveChannel = false; % Load all channels including inactive channels
    [~, channel_info] = vb_load_meg_data(input_file, loadspec);
    
    % Find index of bad channels
    % ch_remove is specified double or cell
    if isa(p.badch, 'double')
        % By index (not ID for now)
        ix_badch = p.badch;
    elseif isempty(p.badch)
        ix_badch = [];
    elseif iscell(p.badch)
        % By channel name (assume cell array of char)
        ix_badch = cell2mat( cellfun(@(x) find(strcmp(channel_info.Name,x)), p.badch, 'UniformOutput',false) );
    end
    
    % Deactivate bad channels
    channel_info.Active(ix_badch) = false;
    
    % Load original file
    B = load(input_file); % To be overwritten

    % Overwrite active status with updated one
    B.MEGinfo.ActiveChannel = channel_info.Active;
    B.MEGinfo.ChannelInfo.Active = channel_info.Active;

    % Overwrite data
    save(input_file, '-struct', 'B');
    disp([input_file ' was overwritten.'])
end
