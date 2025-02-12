function preprocFilter2d(filter_parms, postfix, ext, isplot, datadir, filenames)
% FILTERING TO CONTINUOUS 2D DATA
% 
% > preprocFilter2d(filter_parms, postfix, ext, isplot, datadir, filenames)
%
% Input
%   filter_parms : arrays of filtering parameter structs (see setParm*.m)
%                if this is empty, GUI for setting will appear.
%   postfix: postfix which is added to PROJNAME = savefilename {default:'_filt.mat'}  
%   ext    : file extension for load filenames {default:'.mat'}
%   isplot : plot flag {default:0}
%   datadir, filenames : this values must be given only when filenames are provided from outside of this
%                     routine
%
% NOTE :
%   when different filters are applied to different sets of channels, 
%   the struct array of filter_parms must be input (this is not supported by GUI). 
%   
% 2010/05/18 OY
% * add 'ext' as input argument
% 2010/04/05 OY
% * GUI for setParmFilter   
% 2010/03/15 OY 
% 2010/03/10 Modified by O. Yamashita 


% arguments check
if nargin < 6
    filenames = [];
end
if nargin < 4
    isplot = 0;
end
if nargin < 3 | isempty(ext)
    ext = '*.mat';
end
if nargin < 2 | isempty(postfix)
    postfix = '_filt.mat';
end
if nargin < 1 | isempty(filter_parms)
   filter_parms = setParmFilter;
end

% Get filenames
if isempty(filenames), [datadir, filenames] = getfiles(0, ext); end
for ii = 1 : length(filenames),
    tmp = explode(filenames{ii}, '.');
    PROJNAME{ii} = tmp{1};
end

%
% Main
%
fprintf('--------------------------------------------------\n');
fprintf('    Digital Filtering (in preprocFilter2d.m)      \n');
fprintf('--------------------------------------------------\n');

for jj = 1 : length(PROJNAME),
    fprintf('-----------------------------\n');
    fprintf('      Processing File %2d    \n', jj);
    fprintf('-----------------------------\n');
    
    % Filename
    loadfilename = [datadir filesep PROJNAME{jj}, '.mat']
    savefilename = [datadir filesep PROJNAME{jj}, postfix]

    % Load continuous data
    [data2d, time, channel_name, events, PROJPARM] = load2d(loadfilename);
    
    % Sampling Frequency
    sampling_freq = 1/(time(2) - time(1));

    % Check filter spec.
    if jj == 1
        for nn = 1 : length(filter_parms)
            subCheckFilterSpec(filter_parms(nn), sampling_freq);
        end
    end
    
    % Filtering to multiple channel sets (such as joint recording of EEG
    % and EMG)
    for nn = 1 : length(filter_parms)
        filter_parm = filter_parms(nn);
        [data2d, filter_parm] = subFilter2d(data2d, sampling_freq, filter_parm, isplot);
        
        filter_parm.loadfilename = loadfilename;
        PROJPARM.filter_parm(nn) = filter_parm;
        PROJPARM.history = addcomment(PROJPARM.history,...
            sprintf('[data2d, filter_parm] = subFilter2d(data2d, sampling_freq, filter_parm, isplot);'));
    end  %% filter_parms loop end
  
    % %Save data
    fprintf('Save filtered file : %s ...\n', savefilename);
    save(savefilename, 'data2d', 'time', 'channel_name', 'events', 'PROJPARM');
    fprintf('Done ...\n');
    
end




