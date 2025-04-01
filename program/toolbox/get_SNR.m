function get_SNR(p, data)
% calculate SNR from the 3-dimensional trial data

mdata = mean(data, 3);

%% use mean trial to calculate SNR and write into file
standardDeviation = std(data, 0, 3);


SNRMatrix = abs(mdata ./ standardDeviation); % calculate mean/std for SNR

SNRMatrix = 10 * log10(SNRMatrix); % 10*log_10 for dB

SNRChannelwise = mean(SNRMatrix, 2); % average over time

% open file and write all channel's SNRs + mean over channel into it
folderPath = fullfile(p.fig_root, 'SNR', p.sub);

if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
SNRfile = fopen(fullfile(folderPath, ['SNR_' p.task '.txt']), 'w');
fprintf(SNRfile, '%f\n', SNRChannelwise);
meanSNR = mean(SNRChannelwise);
fprintf(SNRfile, '\nMean: %f\n', meanSNR);
fclose(SNRfile);

