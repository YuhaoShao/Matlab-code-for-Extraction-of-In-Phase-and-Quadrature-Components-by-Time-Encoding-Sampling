%This function is to implement Eq.(52) 
function f = gqt_bandpass2(t1,t2,t3,t4,step,f0,B)
%t1:Initial time of reconstruction interval,
%t2:End time of reconstruction interval, 
%t3:Initial time of k-th indicator function 
%t4:End time of k-th indicator function 
%f0:Carrier frequency
%B:Bandwidth
Fs = 1/step; 
t = t1:1/Fs:t2; 
x = zeros(size(t)); 
x(t >= t3 & t <= t4) = 1./(t4-t3); 
x = -x.*sinpi(2.*f0.*t);
x_add0 = [zeros(1,3*length(x)) x zeros(1,3*length(x))];
x = x_add0;
X = fft(x);
N = length(X);
freq = (0:N-1) * Fs / N; 
H = zeros(1,length(X)); 
H(freq <=B/2) = 1;
H(freq >= Fs - B/2) = 1;
Y = X .* H; 
f1 = ifft(Y);
f = f1(3*length(t)+1:3*length(t)+length(t));
end