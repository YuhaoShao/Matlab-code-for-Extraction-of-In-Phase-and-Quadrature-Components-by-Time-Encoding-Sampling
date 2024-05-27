% This program is to implement bandpass TEM
% and reconstruct I & Q components
% Input:  signal parameters ,TEM parameters(BP_TEM)
% Output: Reconstructed I components(i_f_hufui) and Q
% components(q_f_hufui), SNDR, curves
%==================================================================
% Setting of signal parameters
clc;
close all;
clear;
f0 = 50;       %Carrier frequency
w1 = 2.*pi.*10;
w2 = 2.*pi.*2.5;
w0 = 2.*pi.*f0;
f = @(t) (2.*sin(w1.*t)./(w1.*t)).*cos(w0.*t+sin(w2.*t)./(w2.*t)); 
f_l = f0-15;   %Upper sideband frequency
f_m = f0+15;   %Lower sideband frequency
B=30;          %Bandwidth
endtime = 0.2;
step = 1e-7;
t = -endtime:step:endtime;
f_ori = naninterp(f(t));   % Signa waveform
%==================================================================

%==================================================================
% Settin of TEM parameters  
c =  3;    
cita = 1/(120*1); 
q_pre = [];
q_pre_re = 0;
recordTimes = [];  
recordpo = [];
integralValue = -cita;   %Initial value for the integrator
z = zeros(1,length(t));  % Integrator output
%==================================================================

%==================================================================
% TEM sampling
 for i = 1:length(t)-1
         if t(i)==0 
              integralValue = integralValue + ((2*cos(1))+c)*step;
              q_pre_re = q_pre_re+ (2*cos(1))*step;
         else
             integralValue = integralValue + (f(t(i))+c)*step;
             q_pre_re = q_pre_re+ (f(t(i)))*step;
         end
        z(i) = integralValue ;
        if integralValue >= cita 
            recordTimes = [recordTimes t(i)];
             recordpo = [recordpo i];
             q_pre = [q_pre q_pre_re];
             q_pre_re = 0;
            integralValue = -cita;    
        end     
        if t(i) >= endtime
            break;
        end
 end
tk_1 =  t(recordpo);
tk =  t(recordpo);
ii1=length(tk_1);
%=================================================================

%=================================================================
% Signal Reconstruction
% ---Computing Matrix q--------
q = q_pre(2:end);

% ---Computing function equation (51).(52)-------------
git_func = zeros(length(q),length(t));
gqt_func = zeros(length(q),length(t));
for wai = 1:length(q)
    git_func(wai,:) = git_bandpass2( t(1), t(end),t(recordpo(wai)),t(recordpo(wai+1)),step,f0,B);
    gqt_func(wai,:) = gqt_bandpass2( t(1), t(end),t(recordpo(wai)),t(recordpo(wai+1)),step,f0,B);
end

% ----- Computing Matrix G ----------------- 
G=zeros(length(q),length(q));
for i2 = 1:length(tk)-1
  t_nei2 =  t(recordpo(i2):recordpo(i2+1));
 for i3 =  1:length(tk)-1  
       git = git_func(i3,recordpo(i2):recordpo(i2+1));
       gqt = gqt_func(i3,recordpo(i2):recordpo(i2+1));
       G_jifen = (git.*cos(w0.*t(recordpo(i2):recordpo(i2) ...
           +length(git)-1))-gqt.*sin(w0.*t(recordpo(i2): ...
           recordpo(i2)+length(gqt)-1)));
       G(i2,i3) = trapz(G_jifen).*step;
  end
end
GG = pinv(G);
cjuzhen = GG*q';
t_huifu= -endtime:step:endtime;
t_huifu_range = [-endtime endtime];

% ----Reconstruction------------
i_f_hufui =0;
q_f_hufui =0;
    for tki= 1:length(tk)-1  
      t_hf1 = tk(tki):step:tk(tki+1); 
      git = git_func(tki,:);
      gqt =  gqt_func(tki,:);
       q_f_hufui = cjuzhen(tki).*gqt+ q_f_hufui;
       i_f_hufui = cjuzhen(tki).*git+ i_f_hufui;  
    end  
% ========================================================================

% ========================================================================  
% --- Plotting --------------------------
subplot(221)
hold on
plot(t_huifu,(2.*sin(w1.*t)./(w1.*t)).*cos(sin(w2.*t)./(w2.*t)),'b');
plot(t_huifu,i_f_hufui,'r--');
xlabel('Time(s)');
ylabel('Amplitude');
box on;
lgd=legend('Original','Reconstructed','location','southwest');
lgd.FontSize = 8;
subplot(222)
hold on
plot(t_huifu,(2.*sin(w1.*t)./(w1.*t)).*sin(sin(w2.*t)./(w2.*t)),'b');
plot(t_huifu,q_f_hufui,'r--');
xlabel('Time(s)');
ylabel('Amplitude');
box on;
lgd=legend('Original','Reconstructed','location','southwest');
lgd.FontSize = 8;
subplot(223)
plot(t_huifu,i_f_hufui-(2.*sin(w1.*t)./(w1.*t)).*cos(sin(w2.*t)./(w2.*t)),'b');
xlabel('Time(s)');
ylabel('Amplitude');
box on;
subplot(224)
plot(t_huifu,q_f_hufui-(2.*sin(w1.*t)./(w1.*t)).*sin(sin(w2.*t)./(w2.*t)),'b');
xlabel('Time(s)');
ylabel('Amplitude');

% Error analyses
data = q_f_hufui-naninterp((2.*sin(w1.*t)./(w1.*t)).*sin(sin(w2.*t)./(w2.*t)));
datai = i_f_hufui-naninterp((2.*sin(w1.*t)./(w1.*t)).*cos(sin(w2.*t)./(w2.*t)));
t_data = t_huifu;

% --- Defining the data interval -------------------
percentile = 90;
data_length = length(data);
start_index = round(0.05 * data_length); 
end_index = round(0.95 * data_length); 
middle_90_percent_data = data(start_index:end_index);
middle_90_percent_datai = datai(start_index:end_index);
middle_90_percent_data_t = t_data(start_index:end_index);

% --- Computing SNDR -------------------------
MSE = mean((middle_90_percent_data).^2)./mean((naninterp((2.* ...
    sin(w1.*middle_90_percent_data_t)./(w1.*middle_90_percent_data_t)).* ...
    sin(sin(w2.*middle_90_percent_data_t)./(w2.*middle_90_percent_data_t)))).^2);
MSEi = mean((middle_90_percent_datai).^2)./mean((naninterp((2.* ...
    sin(w1.*middle_90_percent_data_t)./(w1.*middle_90_percent_data_t)).* ...
    cos(sin(w2.*middle_90_percent_data_t)./(w2.*middle_90_percent_data_t)))).^2);
