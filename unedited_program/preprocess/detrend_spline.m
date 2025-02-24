function yrm = detrend_spline(time, y, Nstep)
% Detrend data by spline interpolation

yrm = zeros(size(y));

opt = 'add_end';
doPlot = false;

for ch = 1 : size(y,1)
    ych = y(ch,:);

    switch opt
        case 'mirror'
            dt = time(2)-time(1);
            % Mirroring
            ych_mirror = [ych flip(ych)];
            time_mirror = time(1):dt:(length(time)*2-1)*dt;
            timel = time_mirror(1:Nstep:end);
            % Make spline for mirrored data (twice longer)
            yy = spline(time_mirror,ych_mirror,timel);
            ytrend = interp1(timel,yy,time_mirror,'spline');
            % Remove mirrored part
            ytrend = ytrend(1:length(ych));

        case 'add_end'
            timel = [time(1:Nstep:end) time(end)]; % To prevent the edge-effect
            % timel = [time(1:Nstep:end) time(end)+time(2)-time(1)];
            yy = spline(time,ych,timel);
            ytrend = interp1(timel,yy,time,'spline');

        case 'normal'
            timel = time(1:Nstep:end);
            yy = spline(time,ych,timel);
            ytrend = interp1(timel,yy,time,'spline');

        case 'normal_replace'
            timel = time(1:Nstep:end);
            yy = spline(time,ych,timel);
            ytrend = interp1(timel,yy,time,'spline');
            % Replace end of ytrend with original y to prevent the edge-effect
            % But it contains not detrended data
            ytrend(end-Nstep:end) = ych(end-Nstep:end);
    end

    if doPlot
        figure;
        subplot(2,1,1)
        plot(time,ych,'.',time,ytrend);
        subplot(2,1,2)
        plot(time,ych-ytrend)
    end

    % Store detrended signal
    yrm(ch,:) = ych - ytrend;
end