function get_SNR(p, data, time)

% averaging trials
mdata = mean(data, 3);

% find index of stimulation
[~, stim_index] = min(abs(time));


signalPower = sqrt(mean(mdata(:, stim_index:end) .^2, 2)); % root mean square of signal
noisePower = std(mdata(:, 1:stim_index), 0, 2); % standard deviation of baseline

%signalPower
%noisePower


SNR = 10 * log10(signalPower ./ noisePower);  % dB SNR


%%
% open file and write all channel's SNRs + mean over channel into it
folderPath = fullfile(p.fig_root, 'SNR', p.sub);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
SNRfile = fopen(fullfile(folderPath, ['SNR_' p.task '.txt']), 'w');
fprintf(SNRfile, '%f\n', SNR);
meanSNR = mean(SNR);
meanSNR
fprintf(SNRfile, '\nMean: %f\n', meanSNR);
fclose(SNRfile);

