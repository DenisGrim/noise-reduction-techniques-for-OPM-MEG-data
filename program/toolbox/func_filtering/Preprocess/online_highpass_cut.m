function	[y,a] = online_highpass_cut(x,Fs,Fcut)
% Highpass cut filter using Online lowpass filter

[Nch,T] = size(x);

% Initialize lowpass/highpass filter
[a, z] = lowpass_init(Fs,Fcut,Nch);

y = zeros(Nch,T);

for t=1:T
	[y(:,t), z] = online_highpass(a, x(:, t), z);
end
