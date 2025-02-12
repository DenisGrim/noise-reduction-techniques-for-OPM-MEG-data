function	[y,z] = online_highpass(a,x,z)
% Online calculation of highpass filter
%   [y,z] = online_highpass(a,x,z)
% --- Input
% a  : Coefficient of highpass filter   [D x Nbank]
% x  : Input signal                     [D x 1]
% z  : Lowpass internal state variable  [D x Nbank]
% --- Output
% y  : highpass signal                    [D x Nbank]
% z  : Updated lowpass internal variable  [D x Nbank]
%
% 2007-9-21 Masa-aki Sato

z = a .* z + x * (1 - a(1,:));
y = repmat(x ,[1 size(z,2)])- z;
