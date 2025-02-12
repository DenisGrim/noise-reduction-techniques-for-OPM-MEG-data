function [data, B, A] = pre_bandstop(data, samplefreq, highpassfreq, lowpassfreq, filtertype, order)
% bandstop filtering
% 
% -- Output
% data : bandstop filtered data
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
    B = 1;
    A = 1;
    return
end

switch filtertype
    case 'butter'
        [B,A] = butter(order, [highpassfreq/(samplefreq/2), lowpassfreq/(samplefreq/2)],'stop');
        for nch = 1 : Nch,
            data(nch,:) = filtfilt(B,A,data(nch,:));
        end
    
        
    otherwise
        error('No filtertype. [''butter'']');
end
