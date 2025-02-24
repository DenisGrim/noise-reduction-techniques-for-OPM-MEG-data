function plot_psd(freq, spect, opt)
% freq  : Frequencies.
% spect : Data to be plotted.
%         Unit of spect is defined by sunit.
%
% opt is a struct contains
% dunit : Unit of original data (e.g. fT).
% sunit : Unit of spectrum.
%         Choose from below
%              'a'             : amplitude [unit]
%              'ps'            : power spectrum [unit^2]
%              'psd'           : power spectrum density [unit^2/Hz]
%              'psds'(default) : sqrt(power spectrum density) [unit/sqrt(Hz)]
%              'psddb'         : 10log10(psd) [dB]
%
% chname : If exists, add legend of channel names.
% plot_mean : If true, plot mean (default: true)
% const : Constant line (if noise floor with [fT], constant=15)
% 
% 2022.12.01 (k_suzuki)
%

% Set default values
dunit = 'T';
funit = 'Hz';
sunit = 'psds';
chname = [];
plot_mean = true;
constant = [];
line_color = [];

if exist('opt','var') && ~isempty(opt)
    if isfield(opt, 'dunit')
        dunit = opt.dunit;
    end    
    if isfield(opt, 'funit')
        funit = opt.funit;
    end
    if isfield(opt, 'sunit')
        sunit = opt.sunit;
    end
    if isfield(opt, 'chname')
        chname = opt.chname;
    end
    if isfield(opt, 'plot_mean')
        plot_mean = opt.plot_mean;
    end
    if isfield(opt, 'const')
        constant = opt.const;
    end
    if isfield(opt, 'line_color')
        line_color = opt.line_color;
    end
end


switch sunit
    case 'a'
        labY = ['Amplitude [' dunit ']'];
    case 'ps'
        labY = ['Power spectrum [${' dunit '}^2$]'];
    case 'psd'
        labY = ['PSD $[{' dunit '}^2 /Hz$]'];
    case 'psds'
        labY = ['PSD [$' dunit '/\sqrt{Hz}$]'];
    case 'psddb'
        labY = ['PSD [dB]'];
    case ''
        labY = [];
    otherwise
        labY = [];
end

switch funit
    case 'Hz'
        labX = 'Frequency [Hz]';
    case ''
        labX = [];
    otherwise
        labX = [];
end

if ~isempty(chname)&&plot_mean
    chname{end+1} = 'Mean';
end

%%
% Plot spectrum
if isempty(line_color)
    semilogy(freq, spect, 'LineWidth', 2);
else
    semilogy(freq, spect, 'LineWidth', 2, 'Color', line_color);
end
hold on

% Plot constant line (e.g. noise floor)
if ~isempty(constant)
    xconst = 0:round(freq(end));
    yconst = ones(1,round(freq(end))+1)*constant;
    plot(xconst, yconst, '--k', 'LineWidth', 2);
end

% Plot averaged spectrum
if plot_mean
    semilogy(freq, mean(spect'), 'LineWidth', 3, 'Color','k');
end

% Plot xy labels
if ~isempty(labX)
    xlabel(labX);
end
if ~isempty(labY)
    ylabel(labY, 'interpreter', 'latex');
end
grid on

% Plot channel names
if ~isempty(chname)
    legend(chname, 'Location','northeastoutside')
end

