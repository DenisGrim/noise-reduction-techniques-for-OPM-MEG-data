function prefix_out = reject_chtr(p, prefix_in, trigger_dirname)
% Reject noisy MEG channels and trials
% 'trigger_dirname' is used in case of OPM because trig file contains 
% the information of bad channels and trials.
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

disp(mfilename);

% Set figure properties
close all
set_fig_property(4, 2, 15, 15);

% Define additional prefix of this process
prefix_add = 'r';

% Set prefix of input file
if isempty(prefix_in)
    prefix_in_ = [];
else
    prefix_in_ = [prefix_in '_'];
end

% Set prefix of output file
prefix_out = [prefix_add prefix_in];
prefix_out_ = [prefix_out '_'];

% Threshold for rejecting channels
th_ch = 0.4;% Reject a channel when its ratio of bad trials to all the trials >= th_ch

% Threshold for rejecting trials
th_tr = 1; % Reject a trial when its number of bad channels >= th_tr

for run = 1:p.num_run
    file_name = sprintf('run%02d', run);
    
    % Load OPM-MEG data (includes inactive ch and trials)
    % Pre-defined active info as well
    file_input = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_in_ file_name '.meg.mat']);
    load(file_input, 'bexp', 'MEGinfo')
    data = bexp;
    ActiveChannel = MEGinfo.ActiveChannel;
    ActiveTrial = MEGinfo.ActiveTrial;
    
    Nchan  = size(data, 1);
    Ntrial = size(data, 3);
    
    % Init bad flag
    bad_flag = zeros(Nchan, Ntrial); % #ch x #trial
    
    % Take over the pre-defined bad ch/trial
    bad_flag(~ActiveChannel, :) = 1;
    bad_flag(:, ~ActiveTrial) = 1;
    
    % Load saturation and jump flag in trigger file
    trig_file = fullfile(p.proj_root, trigger_dirname, p.task, [file_name '.trig.mat']);
    load(trig_file, 'trial_sat', 'trial_jump')
    
    % Check bad flags based on sensor saturation
    for cc = 1:Nchan
        bad_flag(cc, trial_sat{cc}) = 1; % trial_sat is saturated flag (generated in make_trial_opm)
        bad_flag(cc, trial_jump{cc}) = 1; % trial_jump is sensor jumping flag (generated in make_trial_opm)
    end
    
    % Check bad flags based on ratio
    [ratio1, ~] = inner_calc_ratio12(data);
    
    % thred is [lb ub] of ratio1
    thred(1) = p.threshold_rejection(1);
    thred(2) = p.threshold_rejection(2);
    thred_std = p.threshold_rejection(3);
    flg = (ratio1 < thred(1)) | (ratio1 > thred(2));
    bad_flag = or(bad_flag, flg);
    
    bad_ratio =  sum(bad_flag, 2) ./ Ntrial;
    bad_ch = bad_ratio >= th_ch; % bad ch for this run
    
    % Bad channel selection based on inter-trial variance of ratio1
    % This is just an add-hoc criterion
    ratio1_std = std(ratio1');
    bad_ch(ratio1_std > thred_std) = 1;
    
    % Find bad channels and trials
    good_ch = find(bad_ch == 0);
    bad_ch = find(bad_ch == 1);
    bad_tr = sum(bad_flag(good_ch , :), 1) >= th_tr;
    bad_tr = find(bad_tr == 1);
    
    disp(['Run-' num2str(run) ': ' num2str(length(bad_ch)) ' channels are rejected.'])
    disp(['Run-' num2str(run) ': ' num2str(length(bad_tr)) ' trials are rejected.'])
    
    % Show result
    h = figure;
    imagesc(bad_flag)
    colormap('gray')
    hold on
    for ch = 1:Nchan
        if ismember(ch, bad_ch)
            plot([0 Ntrial], [ch ch], 'r')
        end
    end
    for tr = 1:Ntrial
        if ismember(tr, bad_tr)
            plot([tr tr], [0 Nchan+1], 'b')
        end
    end
    xlabel('Trial')
    ylabel('Channel')
    title('Bad flag (white area), bad channel (red line), and bad trial (blue line)')
    
    % Save figure
    fig_file = fullfile(p.fig_root, mfilename, p.task, ['Sub-' num2str(p.sub) '_' file_name]);
    vb_savefig_as_shown(h, fig_file)
    disp([fig_file '.png was saved.'])
    
    % Save channel- and trial-rejected data
    MEGinfo.ActiveChannel(bad_ch) = 0;
    MEGinfo.ActiveTrial(bad_tr) = 0;
    MEGinfo.ChannelInfo.Active(bad_ch) = 0;
    new_file = fullfile(p.proj_root, p.dirname.trial, p.task, [prefix_out_ file_name '.meg.mat']);
    copyfile(file_input, new_file);
    vb_save(new_file, 'MEGinfo');
    disp([new_file ' was saved.'])
end

end % End of function

function [ratio1, ratio2] = inner_calc_ratio12(data)
% ratio1 = ymax1./ystd1
% ratio2 = ymax2./ystd2
% ymax1 : max amplitude
% ymax2 : max amplitude of temporally differenced data
% ystd1 : mean amplitude excluding outlier
% ystd2 : mean amplitude of temporally differenced data

[ymax1, ymax2, ystd1, ystd2] = vb_channel_statics(data);

ratio1 = ymax1;
ratio2 = ymax2;

ix1  = find(ystd1 > 0);
ix2  = find(ystd2 > 0);

if size(ystd1, 2) == 1
    ratio1(ix1,:) = vb_repmultiply(ymax1(ix1,:), 1./ystd1(ix1));
    ratio2(ix2,:) = vb_repmultiply(ymax2(ix2,:), 1./ystd2(ix2));
else
    ratio1(ix1) = ymax1(ix1)./ystd1(ix1);
    ratio2(ix2) = ymax2(ix2)./ystd2(ix2);
end


end % End of function
