function data = filter2d(data, sampling_freq, filter_parm, doPlot)
% Apply bandstop, highpass, lowpass, and/or bandpsss filtering to 2D continuous data 
%
% Input
%   data : 2D data [Nch, Nt]  
%   filter_parm : A struct containing filter pareameters (see setParmFilter.m)
%    
% Output 
%   data : Filtered 2D data
%
% 2010/04/14 Modified by O. Yamashita 
% * add bandstop filter
% 2010/03/09 Modified by O. Yamashita file_dialog
% 2023/02/17 Check input parameters by k_suzuki
% 2023/11/22 Y. Takeda modified small details
%
% parameter stuct --> variables in workspace
% struct2vars(filter_parm, {'chanNum',...
%     'lowpass_freq' , 'lowpass_filtertype' , 'lowpass_filterorder',..._
%     'highpass_freq', 'highpass_filtertype', 'highpass_filterorder',..._
%     'bandpass_freq', 'bandpass_filtertype', 'bandpass_filterorder',...
%     'bandstop_freq', 'bandstop_filtertype', 'bandstop_filterorder'});

% Check input parameteres
% If chanNum is not specified, apply filter to all chanenls.
% If filtertype or filterorder is empty, then default value is used.
% If freq is empty, then skip that filtering.
[chanNum, lowpass_freq, lowpass_filtertype, lowpass_filterorder, ...
        highpass_freq, highpass_filtertype, highpass_filterorder, ...
        bandpass_freq, bandpass_filtertype, bandpass_filterorder, ...
        bandstop_freq, bandstop_filtertype, bandstop_filterorder] = inner_check_param(filter_parm, data);

% Baseline correction
data0 = pre_zeroreset(data(chanNum, :));

% Bandstop filter
if ~isempty(bandstop_freq)
    for ii = 1 : size(bandstop_freq, 1)
        [data1, Bs(:,ii),  As(:,ii)] = pre_bandstop(data0, sampling_freq, bandstop_freq(ii,1),...
            bandstop_freq(ii,2), bandstop_filtertype, bandstop_filterorder);
        inner_plot_filter_result(data0,data1, Bs(:,ii), As(:,ii), 'band', bandstop_freq, sampling_freq, doPlot)
        data0 = data1;
    end
else
    Bs = []; As = [];
end

% Highpass filter
if ~isempty(highpass_freq)
    [data1, Bh, Ah] = pre_highpass(data0, sampling_freq, highpass_freq,...
        highpass_filtertype, highpass_filterorder);
    inner_plot_filter_result(data0, data1, Bh, Ah, 'high', highpass_freq, sampling_freq, doPlot)
    data0 = data1;
else
    Bh = []; Ah = [];
end

% Lowpass filter
if ~isempty(lowpass_freq)
    [data1, Bl, Al] = pre_lowpass(data0, sampling_freq, lowpass_freq,...
        lowpass_filtertype, lowpass_filterorder);
    inner_plot_filter_result(data0, data1, Bl, Al, 'low', lowpass_freq, sampling_freq, doPlot)
    data0 = data1;
else
    Bl = []; Al = [];
end

% Bandpass filter
if ~isempty(bandpass_freq)
    [data1, Bb, Ab] = pre_bandpass(data0, sampling_freq, bandpass_freq(1),...
        bandpass_freq(2), bandpass_filtertype, bandpass_filterorder);
    inner_plot_filter_result(data0, data1, Bb, Ab, 'band', bandpass_freq,sampling_freq, doPlot)
    data0 = data1;
else
    Bb = []; Ab = [];
end

%
data(chanNum,:) = data1;

end % End of main function


%% inner functions
function [chanNum, lowpass_freq, lowpass_filtertype, lowpass_filterorder, ...
        highpass_freq, highpass_filtertype, highpass_filterorder, ...
        bandpass_freq, bandpass_filtertype, bandpass_filterorder, ...
        bandstop_freq, bandstop_filtertype, bandstop_filterorder] = inner_check_param(filter_parm, data)

%---- Check num of channels ----
if ~isfield(filter_parm, 'chanNum') || isempty(filter_parm.chanNum)
    chanNum = 1:size(data,1);
else
    chanNum = filter_parm.chanNum;
end
%---------------------------------------

%---- Check low-pass filter params ----
if ~isfield(filter_parm, 'lowpass_freq') || isempty(filter_parm.lowpass_freq)
    lowpass_freq = [];
else
    lowpass_freq = filter_parm.lowpass_freq;
end
if ~isfield(filter_parm, 'lowpass_filtertype') || isempty(filter_parm.lowpass_filtertype)
    lowpass_filtertype = 'butter';
else
    lowpass_filtertype = filter_parm.lowpass_filtertype;
end
if ~isfield(filter_parm, 'lowpass_filtertype') || isempty(filter_parm.lowpass_filterorder)
    lowpass_filterorder = 9;
else
    lowpass_filterorder = filter_parm.lowpass_filterorder;
end
%---------------------------------------

%---- Check high-pass filter params ----
if ~isfield(filter_parm, 'highpass_freq') || isempty(filter_parm.highpass_freq)
    highpass_freq = [];
else
    highpass_freq = filter_parm.highpass_freq;
end
if ~isfield(filter_parm, 'highpass_filtertype') || isempty(filter_parm.highpass_filtertype)
    highpass_filtertype = 'butter';
else
    highpass_filtertype = filter_parm.highpass_filtertype;
end
if ~isfield(filter_parm, 'highpass_filterorder') || isempty(filter_parm.highpass_filterorder)
    highpass_filterorder = 7;
else
    highpass_filterorder = filter_parm.highpass_filterorder;
end
%---------------------------------------

%---- Check band-pass filter params ----
if ~isfield(filter_parm, 'bandpass_freq') || isempty(filter_parm.bandpass_freq)
    bandpass_freq = [];
else
    bandpass_freq = filter_parm.bandpass_freq;
end
if ~isfield(filter_parm, 'bandpass_filtertype') || isempty(filter_parm.bandpass_filtertype)
    bandpass_filtertype = 'butter';
else
    bandpass_filtertype = filter_parm.bandpass_filtertype;
end
if ~isfield(filter_parm, 'bandpass_filterorder') || isempty(filter_parm.bandpass_filterorder)
    bandpass_filterorder = 4;
else
    bandpass_filterorder = filter_parm.bandpass_filterorder;
end
%---------------------------------------

%---- Check band-stop filter params ----
if ~isfield(filter_parm, 'bandstop_freq') || isempty(filter_parm.bandstop_freq)
    bandstop_freq = [];
else
    bandstop_freq = filter_parm.bandstop_freq;
end
if ~isfield(filter_parm, 'bandstop_filtertype') || isempty(filter_parm.bandstop_filtertype)
    bandstop_filtertype = 'butter';
else
    bandstop_filtertype = filter_parm.bandstop_filtertype;
end
if ~isfield(filter_parm, 'bandstop_filterorder') || isempty(filter_parm.bandstop_filterorder)
    bandstop_filterorder = 4;
else
    bandstop_filterorder = filter_parm.bandstop_filterorder;
end
%---------------------------------------
end


function inner_plot_filter_result(databfr, dataaft, B, A, type, fcut, fs, doPlot)
% fs : sampling frequency
% fcut : fcut(1) : low or high, [fcut(1) fcut(2)] : band

if doPlot & ~isempty(B)
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
end