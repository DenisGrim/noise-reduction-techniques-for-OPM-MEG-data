function [data, B, A] = pre_bandpass(data, samplefreq, highpassfreq, lowpassfreq, filtertype, order)
% bandpass filtering
% 
% -- Output
% data : bandpassed filtered data
%  B   : filter weights (numerator parts)
%  A   : filter weights (denominator parts)
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
