% This function is to implement the kernel function of bandlimited POCS 

function f = g_bandlimited1(t1,t2,t3,t4,step,f_m)
% 
%$g\left( t \right)=\frac{{{\pi }_{\left[ {{t}_{k}},...
% {{t}_{k+1}} \right]}}}{{{t}_{k+1}}-{{t}_{k}}}*\frac{...
% \sin \left( \pi Bt \right)}{\pi t}$

%t1:Initial time of reconstruction interval,
%t2:End time of reconstruction interval, 
%t3:Initial time of k-th indicator function 
%t4:End time of k-th indicator function 
%f_m: Upper-edge frequency
Fs = 1/step; 
t = t1:1/Fs:t2; 
x = zeros(size(t)); 
x(t >= t3 & t <= t4) = 1./(t4-t3); 
x_add0 = [zeros(1,2*length(x)) x zeros(1,2*length(x))];
x = x_add0;
X = fft(x);
N = length(X);
freq = (0:N-1) * Fs / N; 
H = zeros(1,length(X)); 
H(freq <= f_m) = 1; 
H(freq >= Fs - f_m) = 1; 
Y = X .* H;
f1 = ifft(Y); 
f = f1(2*length(t)+1:2*length(t)+length(t));
end