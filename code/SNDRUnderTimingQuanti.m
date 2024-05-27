% This program is to simulate bandpass TEM performance under timing
% quantization
clc;
clear;
close all;
mtkl = 100;
NNN1 = [5 6 7 8 9 10 11 12 13 14 15];
MSEI = zeros(length(NNN1),mtkl);
for iwaiwaiwai = 1:length(NNN1)
NNN =NNN1(iwaiwaiwai);
for iwaiwai = 1:mtkl
f0 = 15*90;
w1 = 2.*pi.*10;
w2 = 2.*pi.*2.5;
w0 = 2.*pi.*f0;    
f = @(t) (2.*sin(w1.*t)./(w1.*t)).*cos(w0.*t+sin(w2.*t)./(w2.*t));
f_l = f0-15;
f_m = f0+15;   
endtime = 0.2;
step = 1e-7;
t = -endtime:step:endtime;
f_ori = naninterp(f(t));
c =  3;    
cita = 1/(120*2); 
recordTimes = []; 
recordpo = [];
integralValue = -cita;
z = zeros(1,length(t)); 
 for i = 1:length(t)-1
         integralValue = integralValue + (f_ori(i)+c)*step;
        z(i) = integralValue ;     
        if integralValue >= cita 
            recordTimes = [recordTimes t(i)]; 
             recordpo = [recordpo i]; 
            integralValue =-cita;      
        end       
        if t(i) >= endtime
            break;
        end
 end
tk_ganrao = 4/5*(2*cita);
tk_ganrao = tk_ganrao /(2^NNN);
a_tk = tk_ganrao/4;
data_tk_ganrao = -a_tk + (2*a_tk).*rand(1, length(recordpo));%添加偏置的时间点位置
data_recordpo_ganrao = round(data_tk_ganrao./step);
recordpo = recordpo + data_recordpo_ganrao;
tk_1 =  recordTimes;
tk = tk_1;
ii1 = length(tk_1);
q = zeros(1,length(tk_1)-1);
for i = 1:ii1-1
 q(i) = 2.*cita-c.*(t(recordpo(i+1))-t(recordpo(i)));
end
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
datai = x_huifu-(f_ori);
data_length = length(datai);
start_index = round(0.05 * data_length); 
end_index = round(0.95 * data_length); 
middle_90_percent_datai = datai(start_index:end_index);
middle_90_percent_data_t = t(start_index:end_index);
MSEI(iwaiwaiwai,iwaiwai) = mean(middle_90_percent_datai.^2)/mean(naninterp(f(middle_90_percent_data_t)).^2);
end
end