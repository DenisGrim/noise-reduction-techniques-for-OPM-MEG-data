function [data_z, data_ave, data_std, data_min, data_max] = normalize_data(data, type, data_ave, data_std, data_min, data_max)
% Normalize data with various way
% data : Input data [sensor, time-series]
%
% Type : It decides the normalizetion way.
%        'standardize' (default) : z-scoring (mean=0, var=1)
%        'variance'              : variance norimalization (var=1)
%        'mean'                  : mean normalization (mean=0)
%        'minmax'                : min-max normalization (min=0, max=1)
%        'ch_standardize'        : Channel-wise z-scoring (mean=0, var=1)
%        'ch_variance'           : Channel-wise variance norimalization (var=1)
%        'ch_mean'               : Channel-wise mean normalization (mean=0)
%        'ch_minmax'             : Channel-wise min-max normalization (min=0, max=1)
%
% Return normalized data and average/std.
% If data_ave and data_std are given, data is normalized using them.

if ~exist('type','var')||isempty(type)
    type = 'standardize';
end

if strcmp(type(1:3), 'ch_')
    is_chwise = true;
    type = type(4:end);
else
    is_chwise = false;
end

if ~exist('data_ave','var')||isempty(data_ave)
    if is_chwise
        data_ave = mean(data, 2); % mean of each channel
    else
        data_ave = mean(data(:)); % mean among all-channels
    end
end

if ~exist('data_std','var')||isempty(data_std)
    if is_chwise
        data_std = std(data, [], 2); % std of each channel
    else
        data_std = std(data(:));     % std among all-channels
    end
end

if ~exist('data_min','var')||isempty(data_min)
    if is_chwise
        data_min = min(data, [], 2); % min of each channel
    else
        data_min = min(data(:));     % min among all-channels
    end
end

if ~exist('data_max','var')||isempty(data_max)
    if is_chwise
        data_max = max(data, [], 2); % max of each channel
    else
        data_max = max(data(:)); % max among all-channels
    end
end

switch type
    case 'standardize'
        data_z = (data - data_ave) ./ data_std;
    case 'variance'
        data_z = data ./ data_std;
    case 'mean'
        data_z = data - data_ave;
    case 'minmax'
        data_z = (data - data_min) ./ (data_max - data_min);
end
