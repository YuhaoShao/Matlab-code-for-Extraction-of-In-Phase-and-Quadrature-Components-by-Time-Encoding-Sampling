% This program is used to simulate the performance of bandpass TEM in noise
% environments
clc;
clear;
close all;
mtkl = 100;
RMSE = zeros(1,mtkl);
SNR_IN = [0 5 10 15 20 25 30 35 40 45];
SNDR = zeros(1,length(SNR_IN)); 
for iwaiwaiwai = 1:length(SNR_IN)
    for iwaiwai = 1:mtkl
f0 = 600;
w1 = 2.*pi.*10;
w2 = 2.*pi.*2.5;
w0 = 2.*pi.*f0;
f = @(t) (2.*sin(w1.*t)./(w1.*t)).*sin(w0.*t+sin(w2.*t)./(w2.*t));
f_l = f0-15;
f_m = f0+15;
endtime=0.2;
step=1e-7;
t = -endtime:step:endtime;
f_ori = naninterp(f(t));
snr = SNR_IN(iwaiwaiwai);
signal_noise = awgn(f_ori,snr,'measured');
Fs = 1./step;
noise = signal_noise-f_ori;
X = fft(noise);
N = length(X);
 freq = (0:N-1) * Fs / N; %
 H = zeros(1,length(X)); 
 H(freq <= f_m) = 1; 
 H(freq >= Fs - f_m) = 1;
  H(freq <= f_l) = 0; 
 H(freq >= Fs - f_l) = 0;
Y = X .* H; 
noise_lv = ifft(Y); 
SNR_lv = 10.*log10(mean(f_ori.^2)/mean(noise_lv.^2)); 
SNR_target = snr;
noise_scale_factor = 10^((SNR_lv - SNR_target)/20); 
noise_lv = noise_lv.*noise_scale_factor;
signal_noise_lv = f_ori + noise_lv;
f_ori =  signal_noise_lv;
c =  3;   
cita = 1/(240); 
recordTimes=[]; 
recordpo = [];
q_pre = [];
integralValue = -cita;
q_pre_re = 0;
z=zeros(1,length(t)); 
 for i = 1:length(t)-1
         integralValue = integralValue + (f_ori(i)+c)*step;
         q_pre_re = q_pre_re+(f_ori(i))*step;
        z(i) = integralValue ;     
        if integralValue >= cita 
            recordTimes = [recordTimes t(i)]; 
             recordpo = [recordpo i];
             q_pre = [q_pre q_pre_re];
            integralValue = -cita;   
            q_pre_re = 0;
        end       
        if t(i) >= endtime
            break;
        end
 end
tk_1 =  recordTimes;
tk = tk_1;
ii1 = length(tk_1);
q = zeros(1,length(tk_1)-1);
q = q_pre(2:end);
G = zeros(length(q),length(q));
g_bp_max = zeros(length(tk_1)-1,length(t));
for i  = 1:length(tk_1)-1
  g_bp_max(i,:) = g_bandpass1(t(1),t(end),t(recordpo(i)),t(recordpo(i+1)),step,f_l,f_m); 
end
for i_wai = 1:length(tk_1)-1
 for i = 1:length(tk_1)-1
     G(i_wai,i) = trapz(g_bp_max(i,recordpo(i_wai):recordpo(i_wai+1))).*step;
  end           
end   
GG = pinv(G);
qq = q.';
Cjuzhen = GG*qq;
x_huifu=0;
for i = 1:length(tk_1)-1
     y = Cjuzhen(i).*g_bp_max(i,:);
     x_huifu = x_huifu +y;
end   
 f_ori = naninterp(f(t));
datai = x_huifu-(f_ori);
data_length = length(datai);
start_index = round(0.05 * data_length); 
end_index = round(0.95 * data_length); 
middle_90_percent_datai = datai(start_index:end_index);
middle_90_percent_data_t = t(start_index:end_index);
RMSE (iwaiwai) = mean(middle_90_percent_datai.^2)/mean(naninterp(f(middle_90_percent_data_t)).^2);
    end
   SNDR(iwaiwaiwai) = mean(RMSE);
end
plot(SNR_IN,10.*log10(1./SNDR))
