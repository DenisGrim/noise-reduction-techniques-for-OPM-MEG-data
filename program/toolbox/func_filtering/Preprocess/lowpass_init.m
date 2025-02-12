function	[a,z] = lowpass_init(fs,Fc,D)
% Exponential dump lowpass filter coefficient
%   [a,z] = lowpass_init(fs,Fc)
% --- Input
% fs : Sampling frequencies (Hz)
% Fc : Center frequencies (Hz) for lowpass filter bank   [1 x Nbank]
% D  : Dimension of input signal
% --- Output
% a  : Coefficient of lowpass filter bank                [D x Nbank]
% z  : Internal state variable for lowpass online-filter [D x Nbank]
%    = zeros(D,Nbank)
%
% 2007-9-21 Masa-aki Sato

Fc = Fc(:)';

% Lowpass coefficient
a = exp(-(2*pi)*(Fc/fs));
a = repmat( a, [D 1]);

z = zeros(size(a));
