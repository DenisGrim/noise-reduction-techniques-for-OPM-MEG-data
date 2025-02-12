function [data] = re_normalize_data(data_z, type, data_ave, data_std, data_min, data_max)
% Re-normalize data with various way
% data_z : Normalized data using normalize_data [sensor, time-series]
%
% Type : Specify the same type used in normalized_data.
%        'standardize' (default) : z-scoring (mean=0, var=1)
%        'variance'              : variance norimalization (var=1)
%        'mean'                  : mean normalization (mean=0)
%        'minmax'                : min-max normalization (min=0, max=1)
%        'ch_standardize'        : Channel-wise z-scoring (mean=0, var=1)
%        'ch_variance'           : Channel-wise variance norimalization (var=1)
%        'ch_mean'               : Channel-wise mean normalization (mean=0)
%        'ch_minmax'             : Channel-wise min-max normalization (min=0, max=1)
%
% Return re-normalized data.

if ~exist('type','var')||isempty(type)
    type = 'standardize';
end

if strcmp(type(1:3), 'ch_')
    is_chwise = true;
    type = type(4:end);
else
    is_chwise = false;
end

switch type
    case 'standardize'
        data = (data_z .* data_std) + data_ave;
    case 'variance'
        data = data_z .* data_std;
    case 'mean'
        data = data_z + data_ave;
    case 'minmax'
        data = data_z .* (data_max-data_min) + data_min;
end