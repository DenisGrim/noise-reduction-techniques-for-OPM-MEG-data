function data=online_lowpass(data,sam_rate,cutoff_freq,Norder);

[A, B, C, D] = butter(Norder, cutoff_freq/(sam_rate/2),'low');
Norder = size( A, 1 );
A = A'; B = B'; C = C'; D = D';
		
[Nch ,T] = size(data);

Y = zeros( Nch, T);
X = zeros( Nch, Norder );
for t=1:T
Y(:,t) = X * C + data(:,t) * D;
X      = X * A + data(:,t) * B;
end
data=Y;

