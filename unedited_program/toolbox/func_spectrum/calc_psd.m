function [freq, psdx] = calc_psd(data, fsamp, func_win, len_win_ms, sunit, do_basecorrect)
% data  : Data to be converted to PSD.
% fsamp : Sampling frequency.
% func_win : Window function.
%            Input/output specification of window functon: 
%                   window = func(window_length, options)
%            default: hann
%
% len_win_ms : Window length [ms]
%              If data has multiple trials (dim3~=1), len_win_ms is ignored.
%              default: 1000 [ms]
%
% sunit : Unit of output 'psdx'.
%         Choose from below
%              'a'             : amplitude [unit]
%              'ad'            : amplitude dendsity [unit/Hz]
%              'ps'            : power spectrum [unit^2]
%              'psd'           : power spectrum density [unit^2/Hz]
%              'psds'(default) : sqrt(power spectrum density) [unit/sqrt(Hz)]
%              'psddb'         : 10log10(psd) [dB]
%         Below returns complex value (not power)
%              's'             : spectrum [unit]
%              'sd'            : spectrum density [unit/Hz]
%
% do_basecorrect : If true, correct baseline of single epoch (default: false)
%
% This script is based on spm_opm_psd.m written by Tim Tierney
%
% 2022.12.01 (k_suzuki)
%

if ~exist('func_win', 'var')
    func_win = @hann;
end

if ~exist('len_win_ms', 'var')
    len_win_ms = 1000;
end
    
if ~exist('sunit', 'var')
    sunit = 'psds';
end

if ~exist('do_basecorrect', 'var')
    do_basecorrect = false;
end

%% Set window
% If data has multiple trials (dim3~=1), len_win_ms is ignored
nchan = size(data,1);
if size(data,3)>1
    ntime = size(data,2); % num of time points
    nepoch = size(data,3); % num of epochs (trials)
else
    ntime = round(len_win_ms/1000*fsamp); % num of time points
    begs = 1:ntime:size(data,2); % begin index
    ends = begs+ntime-1;          % end index
    
    % Discard a fraction
    if(ends(end)>size(data,2))
        ends(end) = [];
        begs(end) = [];
    end
    nepoch = length(begs); % num of epochs
end

taper  = window(func_win, ntime); % Make a window (taper) with specified time points
coFac = max(taper) / mean(taper);
taper_rep = repmat(taper, 1, nchan);


%% Compute PSD
freq_res = fsamp/ntime; % freq resolution = fsamp/ntime
freq = 0:freq_res:fsamp/2; 
is_odd = mod(ntime,2)==1;
psdx = zeros(length(freq), nchan);

for j = 1:nepoch
    % If data is already epoched then extract epochs
    if (size(data,3)>1)
        data_epoched = data(:,:,j)';
    else % Otherwise, extract data with specified taper length
        inds = begs(j):ends(j);
        data_epoched = data(:,inds,1)';
    end
    
    % Tapering data and correct for amplitude loss;
    data_tapered = data_epoched.*taper_rep*coFac;
    
    % Baseline correction if required 
    if(do_basecorrect)
        mu = median(data_tapered); % correction using meadian (not mean)
        zf = bsxfun(@minus, data_tapered, mu);
        fzf = zf;
    else
        fzf = data_tapered;
    end
    
    % Fourier transform data and get RMS
    xdft = fft(fzf);
    xdft = xdft(1:floor(ntime/2+1),:);
    
    % Calc power spectrum density (psd)
    switch sunit
        case 'a'
            psdtmp = abs(xdft); % [unit]
        case 'ad'
            psdtmp = abs(xdft)./(ntime*fsamp); % [unit/Hz]
        case 'ps'
            psdtmp = abs(xdft).^2; % [unit^2]
        case {'psd', 'psddb'}
            psdtmp = (abs(xdft).^2)./(ntime*fsamp); % [unit^2/Hz]
        case 'psds'
            psdtmp = abs(xdft)./sqrt(ntime*fsamp); % [unit/sqrt(Hz)]
        case 's'
            % Return complex value
            psdtmp = xdft; % [unit]
        case 'sd'
            % Return complex value
            psdtmp = xdft./(ntime*fsamp); % [unit/Hz]
    end
    
    if(is_odd)
        psdtmp(2:end) = psdtmp(2:end);
    else
        psdtmp(2:end-1) = psdtmp(2:end-1);
    end
    
    % Accumulate avearge PSD to limit memory usage
    psdx = psdx + psdtmp/nepoch;
end

if strcmp(sunit, 'psddb')
    % psdx should be divided by noise floor (baseline)?
    psdx = pow2db(psdx);
end