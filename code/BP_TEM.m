% This code is to implement POCS algorithm for bandpass signals
% Input:  signal parameters ,TEM parameters(BP_TEM)
% Output: Reconstructed signal, SNDR 
%==================================================================
clc;
clear;
close all;
% Setting of signal parameters 
f0 = 50;        %Carrier frequency
w1 = 2.*pi.*10;
w2 = 2.*pi.*2.5;
w0 = 2.*pi.*f0;
f = @(t) (2.*sin(w1.*t)./(w1.*t)).*sin(w0.*t+sin(w2.*t)./(w2.*t));
f_l = f0-15;     %upper-edge frequency
f_m = f0+15;     %Lower-edge frequency
endtime=0.2;
step=1e-7;
t = -endtime:step:endtime;
f_ori = naninterp(f(t));  %Signal waveform
%==================================================================

%==================================================================
% Settin of TEM parameters
c =  3;           %bias
cita = 1/(120);   %threshold 
recordTimes = []; 
recordpo = [];
q_pre = [];
q_pre_re = 0;
integralValue = -cita;  %Initial value for the integrator
z = zeros(1,length(t)); % Integrator output
%==================================================================

%==================================================================
% TEM sampling
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
tk_1 =  recordTimes; %Pulse trigger time
tk = tk_1;
ii1 = length(tk_1);
%=================================================================

%=================================================================
% Signal Reconstruction
% ---Computing Matrix q--------
q = zeros(1,length(tk_1)-1);
q = q_pre(2:end);
G=zeros(length(q),length(q));

% ---Computing function g_bp(t) ----------------------------------
g_bp_max = zeros(length(tk_1)-1,length(t));
for i  = 1:length(tk_1)-1
  g_bp_max(i,:) = g_bandpass1(t(1),t(end),t(recordpo(i)),...
      t(recordpo(i+1)),step,f_l,f_m); 
end

% ----- Computing Matrix G ----------------- 
for i_wai = 1:length(tk_1)-1
 for i = 1:length(tk_1)-1
     G(i_wai,i) = trapz(g_bp_max(i,recordpo(i_wai):recordpo(i_wai+1))).*step;
  end           
end   
GG = pinv(G);
format long;
qq = q.';
Cjuzhen = GG*qq;

% ----Reconstruction------------
x_huifu = 0;
for i = 1:length(tk_1)-1
     y = Cjuzhen(i).*g_bp_max(i,:);
     x_huifu = x_huifu +y;
end   
% ========================================================================

% ========================================================================
% Error analyses
datai = x_huifu-(f_ori);%Error signal
data_length = length(datai);

% --- Defining the data interval -------------------
start_index = round(0.05 * data_length); 
end_index = round(0.95 * data_length); 
middle_90_percent_datai = datai(start_index:end_index);
middle_90_percent_data_t = t(start_index:end_index);

% --- Computing SNDR -------------------------
rmsei = mean(middle_90_percent_datai.^2)/mean(naninterp(f( ...
    middle_90_percent_data_t)).^2);
rmsei_db = 10*log10(rmsei);



figure;
subplot(211)
hold on;
plot(t,x_huifu,'b-');
plot(t,f_ori,'r--');
xlabel('Time(s)');
ylabel('Amplitude');
box on;
lgd=legend('reconstructed','Original','location','southwest');
axis([-endtime endtime -3 3]);
hold off;
subplot(212)
plot(t,z,'b-');
xlabel('Time(s)');
ylabel('Amplitude');
box on;
xlim([-0.1*endtime 0.1*endtime ]);