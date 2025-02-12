function [data, B, A] = pre_highpass(data, samplefreq, highpassfreq, filtertype, order)
% Highpass filtering
% 
% -- Output
% data : bandpassed filtered data
%  B   : filter weights (numerator parts)
%  A   : filter weights (denominator parts)
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
        [data,B] = eegfilt(data, samplefreq,  highpassfreq, 0);
        A = [];
    case 'butter'
        [B,A] = butter(order, highpassfreq/(samplefreq/2),'high');
        for nch = 1 : Nch,
            data(nch,:) = filtfilt(B,A,data(nch,:));
        end
    case 'onlinebutter',
        %         [A,B,C,D] = butter(order, highpassfreq/(samplefreq/2),'high');
        %         A = A'; B = B'; C = C'; D = D';
        %         Y = zeros( Nch, Nt);
        %         X = zeros( Nch, order );
        %         for nt=1:Nt
        %             Y(:,nt) = X * C + data(:,nt) * D;
        %             X      = X * A + data(:,nt) * B;
        %         end
        %         data=Y;
        [B,A] = butter(order, highpassfreq/(samplefreq/2),'high');
        for nch = 1 : Nch,
            data(nch,:) = filter(B,A,data(nch,:));
        end

    case 'onlineexp',
        [data, a] = online_highpass_cut(data, samplefreq, highpassfreq);
        
        % state-space representation
        As = a(1); 
        Bs = 1-a(1);
        Cs = -1;
        Ds = 1;
        
        if exist('ss2tf','file') == 2,
            [B,A] = ss2tf(As,Bs,Cs,Ds);
        else
            warning('ss2tf.m is not found. fiter weights are not returned');
            B = [];
            A = []; 
        end
        
    otherwise
        error('No filtertype. [''eegfilt'', ''butter'', ''onlinebutter'', ''onlineexp'']');
end