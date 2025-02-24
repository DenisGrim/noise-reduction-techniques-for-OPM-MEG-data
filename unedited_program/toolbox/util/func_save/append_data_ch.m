function append_data_ch(file_input, new_ch, file_output)
% Append new data channels (bexp or eegdata) to output file
% If specified new_ch_name is already exist, the ch is overwritten.
%   file_input : Input file.
%   new_ch  : New channel structure contains below fields.
%    .data  : The data of new channel.
%             It must have the matched sample length.
%             Multiple channels can be specified.
%             [Nch, Nsample]
%    .name  : The name of new channel.
%    .type  : The type of new channel.
%             It must be single or list of integer and follow the rule of channel type.
%    .trial : Logical index of active trial.
%             If it is not specified, ActiveTrial is loaded from file_input.
%
%    [For MEG]
%    .pick  : Position of new channel.
%    .Qpick : Orientation of new channel.
%             If these are not specified, they are not updated (existent ch) or set as [0 0 0] (new ch).
%    [For EEG]
%    .Coord : Position of new channel.
%             If these are not specified, they are not updated (existent ch) or set as [0 0 0] (new ch).
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

%% ----- Check inputs -----
% If output file is not specified, input file is overwritten.
if ~exist('file_output','var') || isempty(file_output)
    file_output = file_input;
else
    % If output directory does not exist, it is made.
    [dir_output,~,~] = fileparts(file_output);
    disp('Output file is specified.')
    if ~exist(dir_output,'dir')
        disp('Output dir does not exist.')
        mkdir(dir_output)
        disp([dir_output ' was made.'])
    end
    % If output file is specified, input file is copied.
    copyfile(file_input, file_output, 'f')
    disp(['Input copied to ' file_output])
end

% Load the device type (EEG or MEG)
load(file_input, 'Measurement'); 

% Deploy fields of new_ch 
if isfield(new_ch, 'data')
    new_ch_data = new_ch.data;
else
    error('''new_ch.data'' must be specified.')
end

if isfield(new_ch, 'name')
    new_ch_name = new_ch.name;
else
    error('''new_ch.name'' must be specified.')
end

if isfield(new_ch, 'type')
    new_ch_type = new_ch.type;
else
    error('''new_ch.type'' must be specified.')
end

if isfield(new_ch, 'trial')
    new_ch_active_trial = new_ch.trial;
    if ~islogical(new_ch_active_trial)
        % Binary check (0/1)
        is_binary = sum(new_ch_active_trial==1|new_ch_active_trial==0)==length(new_ch_active_trial);
        if is_binary
            new_ch_active_trial = logical(new_ch_active_trial);
        else
            error('''new_ch.trial'' must be a logical index.')
        end
    end
else
    % Set same ActiveTrial as input file
    INFO = vb_load_measurement_info(file_input, 'MEGINFO');
    new_ch_active_trial = logical(INFO.ActiveTrial);
    % new_ch_active_trial = logical( ones(size(new_ch_data,3),1) );
    % warning(['''new_ch.trial'' is not specified.' newline ...
    %          '''new_ch.data'' must have all trials.'])
end

Nnewch = length(new_ch_type);

switch Measurement
    case 'MEG'
        if isfield(new_ch, 'pick')
            new_ch_pick = new_ch.pick;
            if size(new_ch_pick,1)~=Nnewch
                error('''new_ch.pick'' must be consistent size with other fields.')
            end
            do_update_pick = true;
        else
            % If not specified, set [0 0 0]
            new_ch_pick = zeros(Nnewch,3);
            do_update_pick = false;
        end

        if isfield(new_ch, 'Qpick')
            new_ch_Qpick = new_ch.Qpick;
            if size(new_ch_Qpick,1)~=Nnewch
                error('''new_ch.Qpick'' must be consistent size with other fields.')
            end
            do_update_Qpick = true;
        else
            % If not specified, set [0 0 0]
            new_ch_Qpick = zeros(Nnewch,3);
            do_update_Qpick = false;
        end

    case 'EEG'
        if isfield(new_ch, 'Coord')
            new_ch_Coord = new_ch.Coord;
            if size(new_ch_Coord,1)~=Nnewch
                error('''new_ch.Coord'' must be consistent size with other fields.')
            end
            do_update_Coord = true;
        else
            % If not specified, set [0 0 0]
            new_ch_Coord = zeros(Nnewch,3);
            do_update_Coord = false;
        end
end

if ~iscell(new_ch_name)
    new_ch_name = {new_ch_name};
end

% Define function to check the number is interger or not
IsInteger = @(x) (abs(round(x)-x)) <= eps('double');

%% ----- Start registration -----
disp('Channel registration starts...')
switch Measurement
    case 'MEG'
        % Load the original data
        B = load(file_input, 'bexp', 'MEGinfo', 'pick', 'Qpick');

        % Check consistency of trial
        if length(new_ch_active_trial)~=B.MEGinfo.Nrepeat
            error('Set appropriate index to ''trial'' field')
        end

        for cc=1:Nnewch
            % Check the length of data (#sample)
            if size(new_ch_data,2)~=B.MEGinfo.Nsample
                error('The #sample of new ch is not consistent.')
            end

            % Check the channel type
            if ~isnumeric(new_ch_type(cc)) || ~IsInteger(new_ch_type(cc)) % It must be interger
                error('''new_ch_type'' must be integer.')
            end

            % Check new ch is already registered or not
            if any(strcmp(B.MEGinfo.ChannelInfo.Name, new_ch_name{cc}))
                % If it is registered, overwrite the channel
                is_exist = true;
                new_info_ix = find(strcmp(B.MEGinfo.ChannelInfo.Name, new_ch_name{cc}));
                new_ch_id = B.MEGinfo.ChannelInfo.ID(new_info_ix);
                new_ch_ix = new_info_ix;
                if new_ch_type(cc)~=B.MEGinfo.ChannelInfo.Type(new_ch_ix)
                    warning(['''' new_ch_name{cc} ''' is overwritten with different channel type'])
                end
                %disp(['''' new_ch_name{cc} ''' is overwritten in ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [overwrite]... '])
            else
                % If new ch is not registerd yet, Channel_id = max(Channel_id)+1
                is_exist = false;
                new_info_ix = length(B.MEGinfo.ChannelInfo.Name) + 1;
                new_ch_id = max([B.MEGinfo.ChannelInfo.ID; B.MEGinfo.ExtraChannelInfo.Channel_id]) + 1;
                new_ch_ix = size(B.bexp,1) + 1;
                %disp(['''' new_ch_name{cc} ''' is appended to ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [append]... '])
            end

            % Append new data to bexp
            B.bexp(new_ch_ix,:,new_ch_active_trial) = new_ch_data(cc,:,:);

            % Append new channel data to ExtraChannelInfo
            B.MEGinfo.ChannelInfo.Name(new_info_ix, 1)   = new_ch_name(cc); % Cell to cell
            B.MEGinfo.ChannelInfo.Active(new_info_ix, 1) = 1;
            B.MEGinfo.ChannelInfo.Type(new_info_ix, 1)   = new_ch_type(cc);
            B.MEGinfo.ChannelInfo.ID(new_info_ix, 1)     = new_ch_id;

            % In case of new ch or pick is specified, register it
            if do_update_pick || ~is_exist
                B.pick(new_ch_ix,:) = new_ch_pick(cc,:);
            end
            % In case of new ch or Qpick is specified, register it
            if do_update_Qpick || ~is_exist
                B.Qpick(new_ch_ix,:) = new_ch_Qpick(cc,:);
            end

            % Make new line every 5 ch
            if rem(cc, 5)==0
                fprintf('\n')
            end
        end % End of ch loop

        % Save appended data to the file
        save(file_output, '-struct', 'B', '-append');

    case 'EEG'
        % Load eeg.mat data and convert to meg.mat format
        E = load(file_input);
        B = convert_eeg2meg(E);

        % Check consistency of trial
        if length(new_ch_active_trial)~=B.MEGinfo.Nrepeat
            error('Set appropriate index to ''trial'' field')
        end

        for cc=1:Nnewch
            % Check the length of data (#sample)
            if size(new_ch_data,2)~=B.MEGinfo.Nsample
                error('The #sample of new ch is not consistent.')
            end

            % Check the channel type
            if ~isnumeric(new_ch_type(cc)) || ~IsInteger(new_ch_type(cc)) % It must be interger
                error('''new_ch_type'' must be integer.')
            end

            % Check new ch is already registered or not
            if any(strcmp(B.MEGinfo.ChannelInfo.Name, new_ch_name{cc}))
                % If it is registered, overwrite the channel
                is_exist = true;
                new_info_ix = find(strcmp(B.MEGinfo.ChannelInfo.Name, new_ch_name{cc}));
                new_ch_id = B.MEGinfo.ChannelInfo.ID(new_info_ix);
                new_ch_ix = new_info_ix;
                if new_ch_type(cc)~=B.MEGinfo.ChannelInfo.Type(new_ch_ix)
                    warning(['''' new_ch_name{cc} ''' is overwritten with different channel type'])
                end
                %disp(['''' new_ch_name{cc} ''' is overwritten in ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [overwrite]... '])
            else
                % If new ch is not registerd yet, Channel_id = max(Channel_id)+1
                is_exist = false;
                new_info_ix = length(B.MEGinfo.ChannelInfo.Name) + 1;
                new_ch_id = max([B.MEGinfo.ChannelInfo.ID; B.MEGinfo.ExtraChannelInfo.Channel_id]) + 1;
                new_ch_ix = size(B.bexp,1) + 1;
                %disp(['''' new_ch_name{cc} ''' is appended to ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [append]... '])
            end

            % Append new data to bexp
            B.bexp(new_ch_ix,:,new_ch_active_trial) = new_ch_data(cc,:,:);

            % Append new channel data to ExtraChannelInfo
            B.MEGinfo.ChannelInfo.Name(new_info_ix, 1)   = new_ch_name(cc); % Cell to cell
            B.MEGinfo.ChannelInfo.Active(new_info_ix, 1) = 1;
            B.MEGinfo.ChannelInfo.Type(new_info_ix, 1)   = new_ch_type(cc);
            B.MEGinfo.ChannelInfo.ID(new_info_ix, 1)     = new_ch_id;

            % In case of new ch or coord is specified, register it to pick
            % (pick converted to EEGinfo.Coord)
            if do_update_Coord || ~is_exist
                B.pick(new_ch_ix, :) = new_ch_Coord(cc,:);
            end

            % Make new line every 5 ch
            if rem(cc, 5)==0
                fprintf('\n')
            end
        end % End of ch loop

        % Re-convert meg.mat to eeg.mat format
        E = reconvert_meg2eeg(B);

        % Save appended data to the file
        save(file_output, '-struct', 'E', '-append');

end % End of measurement swtich

% End of registration.
if rem(cc, 5)~=0
    fprintf('\n')
end        
disp('Finished.')