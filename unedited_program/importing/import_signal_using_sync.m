function [sig_new, sync_onset_new, sync_onset_src, sync_onset_dst] = import_signal_using_sync(sync_src, sync_dst, sig, ignore_mismatch)
% Import the signal from source to destination modality using sync pulse.
% Shift the time index of signal by matching the mean diff of sync pulse.
%
% Note : Make sure that sampling rate of both modality must be same.
%
%  Inputs
%   sync_src : Sync signal (must be pulse) of source modality.
%   sync_dst : Sync signal (must be pulse) of destination modality.
%   sig      : Importing signal of source modality.
%   ignore_mismatch : Ignore the mis-matched num of sync signal.
%                     'first' : Discard the first sync with having greater #sync
%                     'end'   : Discard the end sync with having greater #sync
%                     false   : No ignoring (detect mismatch as an error)
%
%  Outputs
%   sig_new        : Imported signal from source to destination modality.
%   sync_onset_new : Imported sync onset from source to destination modality.
%              It is used for checking the importing.
%
%    e.g., In case of importing EOG from EEG to MEG
%          sync_src = sync signal of EEG
%          sync_dst = sync signal of MEG
%          sig      = EOG recorded by EEG
%          sig_new  = EOG imported to MEG
%
% Copyright (C) 2011, ATR All Rights Reserved.
% License : New BSD License(see VBMEG_LICENSE.txt)

if ~exist('ignore_mismatch','var')||isempty(ignore_mismatch)
    ignore_mismatch = false;
end

Nsample_src = size(sync_src, 2);
Nsample_dst = size(sync_dst, 2);
Nch_sig = size(sig, 1);

% --- Check sig is consistent with sync_src ---
if Nsample_src~=size(sig,2)
    error(['Size of importing signal does not consist with source modality.' newline ...
        'Check the order of input variables.'])
end

% --- Binarize sync pulse ---
% - source
sync_src = sync_src./max(sync_src);
sync_src = sync_src > 0.5;
% - destination
sync_dst = sync_dst./max(sync_dst);
sync_dst = sync_dst > 0.5;

% --- Extract onsets of sync pulse ---
% - source
bindiff = sync_src(1:end-1) - sync_src(2:end);
ix_onset_src = find(bindiff==-1)+1; % index of onsets
% - destination
bindiff = sync_dst(1:end-1) - sync_dst(2:end);
ix_onset_dst = find(bindiff==-1)+1; % index of onsets

% --- Check num of sync pulse ---
Nonset_src = length(ix_onset_src);
Nonset_dst = length(ix_onset_dst);
dif_sync = Nonset_src-Nonset_dst;
disp(['Detected pulses: SRC(' num2str(Nonset_src) '), DST(' num2str(Nonset_dst) ')'])

if dif_sync ~= 0
    switch ignore_mismatch
        case 'first'
            % Discard the first sync of src or dst (whichever is greater)
            if dif_sync>0
                ix_onset_src(1:abs(dif_sync)) = [];
            else
                ix_onset_dst(1:abs(dif_sync)) = [];
            end
        case 'end'
            % Discard the end sync of src or dst (whichever is greater)
            if dif_sync>0
                ix_onset_src(end-abs(dif_sync)+1:end) = [];
            else
                ix_onset_dst(end-abs(dif_sync)+1:end) = [];
            end
        case false
            error('Inconsistent sync pulses were detected between source and destination modalities.')
    end
end

% --- Make binary time-series of sync onsets ---
sync_onset_src = zeros(size(sync_src));
sync_onset_dst = zeros(size(sync_dst));
sync_onset_src(ix_onset_src) = 1;
sync_onset_dst(ix_onset_dst) = 1;

% --- Import signal to destination by shifting mean diff of sync pulse ---
% - Init new signal in destination
sig_new = zeros(Nch_sig, Nsample_dst);
sync_onset_new = zeros(1, Nsample_dst); % For checking results

% - Calc mean diff of onset timing
diff_onset = round(mean(ix_onset_dst - ix_onset_src));

if diff_onset > 0 
    % recoding of src started later than dst
    offset = diff_onset + 1;
    if Nsample_dst > offset + Nsample_src
        tmp = 0;
    else
        tmp = offset + Nsample_src - Nsample_dst;
    end
    sig_tmp = sig(:, 1:end-tmp);
    sig_new(:, offset:offset+size(sig_tmp,2)-1) = sig_tmp;
    sync_tmp = sync_onset_src(:, 1:end-tmp);
    sync_onset_new(:, offset:offset+size(sync_tmp,2)-1) = sync_tmp;
else
    % recoding of src started earlier than dst
    offset = abs(diff_onset) + 1;
    if Nsample_src - offset > Nsample_dst
        tmp = Nsample_src - offset - Nsample_dst + 1;
    else
        tmp = 0;
    end
    sig_tmp = sig(:, offset:end-tmp);
    sig_new(:, 1:size(sig_tmp,2)) = sig_tmp;
    sync_tmp = sync_onset_src(:, offset:end-tmp);
    sync_onset_new(:, 1:size(sync_tmp,2)) = sync_tmp;
end



