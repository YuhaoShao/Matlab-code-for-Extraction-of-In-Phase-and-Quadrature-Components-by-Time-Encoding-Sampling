% This code is to valicdate Theorem 1 and Theorem 2. 
%==================================================================
clc;
clear;
close all;
cdcd = 1:0.1:10;
SNDR = zeros(1,length(cdcd));
for iwaiwai = 1:length(cdcd)
 f0 = 15.*cdcd(iwaiwai);
w1 = 2.*pi.*10;
w2 = 2.*pi.*2.5;
w0 = 2.*pi.*f0;
f = @(t) (2.*sin(w1.*t)./(w1.*t)).*cos(w0.*t+sin(w2.*t)./(w2.*t));
f_l = f0-15;
f_m = f0+15;
endtime=0.2;
step=1e-7;
t = -endtime:step:endtime;
f_ori = naninterp(f(t));
c =  3;    
cita = 1/(120*1);  
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
             q_pre =[q_pre q_pre_re];
            integralValue =-cita;   
            q_pre_re = 0;
        end       
        if t(i) >= endtime
            break;
        end
 end
tk_1 =  recordTimes;
tk = tk_1;
ii1=length(tk_1);
q=zeros(1,length(tk_1)-1);
q = q_pre(2:end);
G=zeros(length(q),length(q));
g_bp_max= zeros(length(tk_1)-1,length(t));
for i  = 1:length(tk_1)-1
  g_bp_max(i,:)= g_bandpass1(t(1),t(end),t(recordpo(i)),t(recordpo(i+1)),step,f_l,f_m); 
end
for i_wai = 1:length(tk_1)-1
 for i = 1:length(tk_1)-1
     G(i_wai,i)=trapz(g_bp_max(i,recordpo(i_wai):recordpo(i_wai+1))).*step;
  end           
end   
GG = pinv(G);
format long;
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
SNDR(iwaiwai) = mean(middle_90_percent_datai.^2)/mean(naninterp(f(middle_90_percent_data_t)).^2);
disp(iwaiwai);
end
plot(cdcd,10.*log10(SNDR));