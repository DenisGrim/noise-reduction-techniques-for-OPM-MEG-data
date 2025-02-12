function append_extra_ch(file_input, new_ch, file_output)
% Append new extra channels (bexp_ext or refmg) to output file
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
%             >0 is registered as ReferenceChannel (refmg).
%             <=0 is registered as ExtraChannel (bexp_ext).
%    .trial : Logical index of active trial.
%             If it is not specified, ActiveTrial is loaded from file_input.
%
%    [For MEG]
%    .pick  : Position of new channel.
%    .Qpick : Orientation of new channel.
%             For ExtraChannel, these are not refered.
%             For ReferenceChannel, these are registered in 
%             'ref_pick' and 'ref_Qpick'.
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

% Define function to check whether number is interger or not
IsInteger = @(x) (abs(round(x)-x)) <= eps('double');

%% ----- Start registration -----
disp('Channel registration starts...')
switch Measurement
    case 'MEG'
        % Load the original data
        B = load(file_input, 'bexp_ext', 'refmg', 'MEGinfo', 'ref_pick', 'ref_Qpick');

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
            if isnumeric(new_ch_type(cc)) && IsInteger(new_ch_type(cc)) % It must be interger
                if new_ch_type(cc) > 0
                    ch_type = 'ref';
                else
                    ch_type = 'ext';
                end
            else
                error('''new_ch_type'' must be integer.')
            end

            % Check new ch is already registered or not
            if any(strcmp(B.MEGinfo.ExtraChannelInfo.Channel_name, new_ch_name{cc}))
                % If it is registered, overwrite the channel
                is_exist = true;
                new_info_ix = find(strcmp(B.MEGinfo.ExtraChannelInfo.Channel_name, new_ch_name{cc}));
                new_ch_id = B.MEGinfo.ExtraChannelInfo.Channel_id(new_info_ix);
                % Get list of channel index
                switch ch_type
                    case 'ext'
                        % ExtraChannel have Channel_type<=0
                        list_ex = B.MEGinfo.ExtraChannelInfo.Channel_type<=0;
                    case 'ref'
                        % ReferenceChannel have Channel_type>0
                        list_ex = B.MEGinfo.ExtraChannelInfo.Channel_type>0;
                end
                % Get list of channel name
                list_exname = B.MEGinfo.ExtraChannelInfo.Channel_name(list_ex);
                % Find corresponds index in bexp_ext or refmg
                new_ch_ix = vb_util_get_index(list_exname, new_ch_name(cc));
                if isempty(new_ch_ix)
                    error(['Specified channel ''' new_ch_name{cc} ''' is already registered as different channel type.'])
                end
                %disp(['''' new_ch_name{cc} ''' is overwritten in ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [overwrite]... '])
            else
                is_exist = false;
                % If new ch is not registerd yet, Channel_id = max(Channel_id)+1
                new_info_ix = length(B.MEGinfo.ExtraChannelInfo.Channel_name) + 1;
                new_ch_id = max([B.MEGinfo.MEGch_id; B.MEGinfo.ExtraChannelInfo.Channel_id]) + 1;
                % Find corresponds index in bexp_ext or refmg
                switch ch_type
                    case 'ext'
                        new_ch_ix = size(B.bexp_ext,1) + 1;
                    case 'ref'
                        new_ch_ix = size(B.refmg,1) + 1;
                end
                %disp(['''' new_ch_name{cc} ''' is appended to ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [append]... '])
            end

            % Append new data to bexp_ext or refmg
            switch ch_type
                case 'ext'
                    B.bexp_ext(new_ch_ix,:,new_ch_active_trial) = new_ch_data(cc,:,:);
                case 'ref'
                    B.refmg(new_ch_ix,:,new_ch_active_trial) = new_ch_data(cc,:,:);
            end
            % Append new channel data to ExtraChannelInfo
            B.MEGinfo.ExtraChannelInfo.Channel_name(new_info_ix, 1)   = new_ch_name(cc); % Cell to cell
            B.MEGinfo.ExtraChannelInfo.Channel_active(new_info_ix, 1) = 1;
            B.MEGinfo.ExtraChannelInfo.Channel_type(new_info_ix, 1)   = new_ch_type(cc);
            B.MEGinfo.ExtraChannelInfo.Channel_id(new_info_ix, 1)     = new_ch_id;

            % In case of new ref ch or pick is specified, register it
            if strcmp(ch_type,'ref') && (do_update_pick || ~is_exist)
                B.ref_pick(new_ch_ix, :) = new_ch_pick(cc,:);
            end
            % In case of new ref ch or Qpick is specified, register it
            if strcmp(ch_type,'ref') && (do_update_Qpick || ~is_exist)
                B.ref_Qpick(new_ch_ix, :) = new_ch_Qpick(cc,:);
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
            if isnumeric(new_ch_type(cc)) && IsInteger(new_ch_type(cc)) % It must be interger
                ch_type = 'ext';
            else
                error('''new_ch_type'' must be integer.')
            end

            % Check new ch is already registered or not
            if any(strcmp(B.MEGinfo.ExtraChannelInfo.Channel_name, new_ch_name{cc}))
                % If it is registered, overwrite the channel
                is_exist = true;
                new_info_ix = find(strcmp(B.MEGinfo.ExtraChannelInfo.Channel_name, new_ch_name{cc}));
                new_ch_id = B.MEGinfo.ExtraChannelInfo.Channel_id(new_info_ix);
                % Get list of channel index
                switch ch_type
                    case 'ext'
                        list_ex = B.MEGinfo.ExtraChannelInfo.Channel_type>0;
                end
                % Get list of channel name
                list_exname = B.MEGinfo.ExtraChannelInfo.Channel_name(list_ex);
                % Find corresponds index in bexp_ext or refmg
                new_ch_ix = vb_util_get_index(list_exname, new_ch_name(cc));
                if isempty(new_ch_ix)
                    error(['Specified channel ''' new_ch_name{cc} ''' is already registered as different channel type.'])
                end
                %disp(['''' new_ch_name{cc} ''' is overwritten in ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [overwrite]... '])
            else
                is_exist = false;
                % If new ch is not registerd yet, Channel_id = max(Channel_id)+1
                new_info_ix = length(B.MEGinfo.ExtraChannelInfo.Channel_name) + 1;
                new_ch_id = max([B.MEGinfo.MEGch_id; B.MEGinfo.ExtraChannelInfo.Channel_id]) + 1;
                new_ch_ix = size(B.bexp_ext,1) + 1;
                %disp(['''' new_ch_name{cc} ''' is appended to ' file_output])
                fprintf(['''' new_ch_name{cc} ''' [append]... '])
            end

            % Append new data to bexp
            switch ch_type
                case 'ext'
                    B.bexp_ext(new_ch_ix,:,new_ch_active_trial) = new_ch_data(cc,:,:);
            end

            % Append new channel data to ExtraChannelInfo
            B.MEGinfo.ExtraChannelInfo.Channel_name(new_info_ix, 1)   = new_ch_name(cc); % Cell to cell
            B.MEGinfo.ExtraChannelInfo.Channel_active(new_info_ix, 1) = 1;
            B.MEGinfo.ExtraChannelInfo.Channel_type(new_info_ix, 1)   = new_ch_type(cc);
            B.MEGinfo.ExtraChannelInfo.Channel_id(new_info_ix, 1)     = new_ch_id;

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