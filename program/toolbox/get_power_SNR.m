function get_power_SNR(p, data)
% calculate SNR from the 3-dimensional trial data

mdata = mean(data, 3);

% use mean trial to calculate SNR and write into file
standardDeviation = std(data, 0, 3);


signalPower = var(mdata, 0, 2);  % Power of the average signal (per channel)
noisePower = mean(standardDeviation .^ 2, 2);  % Mean variance across trials (per channel)

SNR = 10 * log10(signalPower ./ noisePower);  % dB SNR


%%
% open file and write all channel's SNRs + mean over channel into it
folderPath = fullfile(p.fig_root, 'SNR', p.sub);
if ~exist(folderPath, 'dir')
    mkdir(folderPath);
end
SNRfile = fopen(fullfile(folderPath, ['power_SNR_' p.task '.txt']), 'w');
fprintf(SNRfile, '%f\n', SNR);
meanSNR = mean(SNR);
fprintf(SNRfile, '\nMean: %f\n', meanSNR);
fclose(SNRfile);

