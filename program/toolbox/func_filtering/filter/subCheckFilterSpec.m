function subCheckFilterSpec(filter_parm, sampling_freq)

struct2vars(filter_parm, {'chanNum',...
    'lowpass_freq' , 'lowpass_filtertype' , 'lowpass_filterorder',..._
    'highpass_freq', 'highpass_filtertype', 'highpass_filterorder',..._
    'bandpass_freq', 'bandpass_filtertype', 'bandpass_filterorder',...
    'bandstop_freq', 'bandstop_filtertype', 'bandstop_filterorder'});

% Bandstop filter
for ii = 1 : size(bandstop_freq, 1)
    [tmp,Bs, As] = pre_bandstop([], sampling_freq , bandstop_freq(ii,1),...
        bandstop_freq(ii,2), bandstop_filtertype, bandstop_filterorder);
    if Bs~=1 | As ~= 1,
        check_filter_spec(Bs,As,sampling_freq);
        set(gcf, 'name', 'Bandstop filter');
    end
end

% Highpass filter
if strcmp(highpass_filtertype, 'eegfilt') ~=1
    [tmp, Bh, Ah] = pre_highpass(tmp, sampling_freq , highpass_freq,...
        highpass_filtertype, highpass_filterorder);
    if ~isempty(Bh),
        check_filter_spec(Bh,Ah,sampling_freq,512,100,highpass_freq*10);
        set(gcf, 'name', 'Highpass filter');
    end
end

% Lowpass filter
if strcmp(lowpass_filtertype, 'eegfilt') ~=1
    [tmp, Bl, Al] = pre_lowpass(tmp, sampling_freq , lowpass_freq,...
        lowpass_filtertype, lowpass_filterorder);
    if ~isempty(Bl),
        check_filter_spec(Bl,Al,sampling_freq);
        set(gcf, 'name', 'Lowpass filter');
    end
end

% Bandpass filter
if strcmp(bandpass_filtertype, 'eegfilt') ~=1
    [tmp,Bb, Ab] = pre_bandpass(tmp, sampling_freq , bandpass_freq(1),...
        bandpass_freq(2), bandpass_filtertype, bandpass_filterorder);
    if ~isempty(Bb),
        check_filter_spec(Bb,Ab,sampling_freq);
        set(gcf, 'name', 'Bandpass filter');
    end
end

pause(0.5);