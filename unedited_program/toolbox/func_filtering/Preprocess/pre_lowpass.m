function [data,B,A] = pre_lowpass(data, samplefreq, lowpassfreq, filtertype, order)
% Lowpass filtering
% 
% -- Example
% > [data,B,A] = pre_lowpass(data, 1024, 100, 'butter', 5)
%
% -- Output
% data : lowpassed filtered data
%  B   : filter weights (numerator parts)
%  A   : filter weights (denominator parts)
%
% 2009/11/20 OY
% * nothing done when filtertype = [];
% 2009/10/15 OY

[Nch,Nt] = size(data);
if isempty(filtertype),
    data = data;
    B = [];
    A = [];
    return
end

switch filtertype
    case 'eegfilt', 
        [data,B] = eegfilt(data, samplefreq,  0, lowpassfreq);
        A = 1;
    case 'butter'
        [B,A] = butter(order, lowpassfreq/(samplefreq/2),'low');
        for nch = 1 : Nch,
            data(nch,:) = filtfilt(B,A,data(nch,:));
        end
    case 'onlinebutter',
%         [A,B,C,D] = butter(order, lowpassfreq/(samplefreq/2),'low');
%         A = A'; B = B'; C = C'; D = D';
%         Y = zeros( Nch, Nt);
%         X = zeros( Nch, order );
%         for nt=1:Nt
%             Y(:,nt) = X * C + data(:,nt) * D;
%             X      = X * A + data(:,nt) * B;
%         end
%         data=Y;

        [B,A] = butter(order, lowpassfreq/(samplefreq/2),'low');
        for nch = 1 : Nch,
            data(nch,:) = filter(B,A,data(nch,:));
        end
        
    otherwise
        error('No filtertype. {''eegfilt'', ''butter'', ''onlinebutter''}');
end