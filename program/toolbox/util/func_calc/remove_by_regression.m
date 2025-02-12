function datarm = remove_by_regression(data, refdata, parm_denoise, doPlot)
% Regress out noise signal from sensor signals
%   datarm = remove_by_regression(
% --- Input
%   data    : Data matrix you want to denoise [ch time] 
%   refdata : Reference data to compute denoised components 
%   parm_denoise : parameter structure for denoising (optional)
%                field:
%                       'bandpass_freq' : Band-pass frequencies for filtered explanatory variables
%                                         Multiple band-pass frequencies can be specified as
%                                         [high-pass1 low-pass1; high-pass2 low-pass2; high-pass3 low-pass3]     
%                       'sampling_freq' : Sampling frequency (Hz)
%                       'bandpass_filtertype'  : (Default: 'butter')
%                       'bandpass_filterorder' : (Default: 4)
%                       'add_derivative'       : Logical value (Default: OFF)
%                       'add_average'          : Logical value (Default: OFF)
%             
% -- Output
%   datarm : denoised data [ch time]
%
% 2021-10-29 k_suzuki

% Check parameters
if exist('parm_denoise', 'var')
    [add_deriv, add_average, bp, bp_filt, bp_order, fs] = inner_check_parm(parm_denoise);
else
    bp = []; % No filtering
    add_deriv = OFF; % No temporal derivatives
    add_average = OFF; % No averaged signal
end

if exist('doPlot', 'var')
    if isempty(doPlot)
        doPlot = OFF;
    end
else
    doPlot = OFF;
end

%
Nt = size(data,2);
Nch = size(data,1);
Nrch = size(refdata,1);
Y = data';

% Prepare explanatory variables using refdata
if isempty(bp)
    varexp = refdata';
else
    varexp = zeros(Nt, Nrch*size(bp,1));
    for ff=1:size(bp,1)
        [filtered,Bb,Ab] = inner_bandpass(refdata, fs , bp(ff,1), bp(ff,2), bp_filt, bp_order);
        varexp(:,1+Nrch*(ff-1):1+Nrch*(ff-1)+(Nrch-1)) = filtered';
        if doPlot, inner_plot_filter_result(refdata,filtered,Bb,Ab,'band',bp(ff,:),fs,1); end
    end
end

% Add 1st order temporal derivatives to varexp
if(add_deriv)
    dref = diff(refdata');
    dref = [dref(1,:); dref];
    varexp = [varexp dref];
end

% Add averaged signal
if(add_average)
    gs = mean(Y,2); % Averaged signal among sensor ch
    varexp = [varexp gs];
end

% Construct design matrix
X = [ones(Nt,1) varexp]; % Design matrix

% Normalization each column
Xn = bsxfun(@rdivide, X, sqrt(sum(X.^2, 1)));

% Solve least-square
%beta = pinv(X'*X) * X' * Y;
beta = pinv(Xn) * Y;

% Regress-out reference noise
% Residual = meg signal without reference noise
Yres = Y - Xn*beta; 
datarm = Yres';

end % End of Function


function [add_derivative, add_average, bandpass_freq, bandpass_filtertype, bandpass_filterorder, sampling_freq] = inner_check_parm(parm)

if isempty(parm)
    bandpass_freq = []; % No filtering
    bandpass_filtertype=[]; bandpass_filterorder=[]; sampling_freq=[];
    add_derivative = []; % No temporal derivative
    add_average = []; % No averaging
else
    % Extract fields for filtering
    field = {'bandpass_freq', 'sampling_freq', 'bandpass_filtertype', 'bandpass_filterorder'};
    exist_field = isfield(parm, field);
    
    % bandpass_freq and smpling_freq are necessary
    if ~exist_field(1) || ~exist_field(2)
        bandpass_freq = []; % No filtering
        bandpass_filtertype=[]; bandpass_filterorder=[]; sampling_freq=[];
    else
        bandpass_freq = parm.bandpass_freq;
        sampling_freq = parm.sampling_freq;
        
        if exist_field(3)
            bandpass_filtertype = parm.bandpass_filtertype;
            if isempty(bandpass_filtertype), bandpass_filtertype = 'butter'; end % Set default parm
        else
            bandpass_filtertype = 'butter'; % Set default parm
        end
        
        if exist_field(4)
            bandpass_filterorder = parm.bandpass_filterorder;
            if isempty(bandpass_filterorder), bandpass_filterorder = 4; end % Set default parm
        else
            bandpass_filterorder = 4; % Set default parm
        end
    end
    
    % Extract field for temporal derivative
    field = {'add_derivative'};
    exist_field = isfield(parm, field);
    if exist_field(1)
        add_derivative = parm.add_derivative;
    else
        add_derivative = OFF;
    end

    % Extract field for averaged signal
    field = {'add_average'};
    exist_field = isfield(parm, field);
    if exist_field(1)
        add_average = parm.add_average;
    else
        add_average = OFF;
    end
end
end % End of Function


function [data, B, A] = inner_bandpass(data, samplefreq, highpassfreq, lowpassfreq, filtertype, order)
% bandpass filtering
% 
% -- Output
% data : bandpassed filtered data
%  B   : filter weights (numerator parts)
%  A   : filter weights (denominator parts). If IIR filter is used, A is empty.
%
%
% 2009/11/20 OY
% * nothing done when filtertype = [];
% 2009/09/02 OY

[Nch,Nt] = size(data);

if isempty(filtertype),
    data = data;
    B = [];
    A = [];
    return
end

switch filtertype
    case 'eegfilt', 
        [data,B] = eegfilt(data, samplefreq,  highpassfreq, lowpassfreq);
        A = [];
    case 'butter'
        [B,A] = butter(order, [highpassfreq/(samplefreq/2), lowpassfreq/(samplefreq/2)],'bandpass');
        for nch = 1 : Nch,
            data(nch,:) = filtfilt(B,A,data(nch,:));
        end
    case 'onlinebutter'
        [B,A] = butter(order, [highpassfreq/(samplefreq/2), lowpassfreq/(samplefreq/2)],'bandpass');
        for nch = 1 : Nch,
            data(nch,:) = filter(B,A,data(nch,:));
        end
        
    otherwise
        error('No filtertype. [''eegfilt'', ''butter'', ''onlinebutter'']');
end
end % End of Function


function inner_plot_filter_result(databfr,dataaft,B,A,type,fcut,fs,isplot)
% fs : sampling frequency
% fcut : fcut(1) : low or high, [fcut(1) fcut(2)] : band

if isplot & ~isempty(B)
    switch	type
        case	'low'
            Fmax = fcut(1)*3;
            Tmax = 10/fcut(1);
        case	'high'
            Fmax = fcut(1)*5;
            Tmax = 1/fcut(1);
        case	'band'
            Fmax = fcut(2)*3;
            Tmax = 5/fcut(1);
    end
    
    t  = 0:1/fs:Tmax;
    f  = [0:0.001:1]*Fmax;
    X  = zeros(length(t),1);
    X(2) = 1;
    
    H = freqz(B,A,f,fs);
    X = filter(B,A, X );

    figure,
    subplot(3,2,[1 2])
    plot(databfr(1,1:1000),'linewidth',2); hold on; plot(dataaft(1,1:1000),'r-');
    title('timeseries example (ch=1)');
    subplot 323
    plot(f,abs(H))
    xlabel('Freq[Hz]')
    title('Gain of frequency response')
    ylim([ 0  1.1])
    subplot 324
    plot(f,angle(H))
    xlabel('Freq[Hz]')
    title('Angle of frequency response')
    subplot 325
    plot(t,X)
    xlabel('Time[sec]')
    title('Impulse response')
    subplot 326
    plot(B)
    hold on
    plot(A,'--')
    xlabel('Tap number')
    title('Filter coefficient')
end
end % End of Function