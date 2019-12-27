clear
clc
T =0.1;   %周期1s
fs = 1 / T;
%巴特沃斯滤波器的设计指标参数
Wp = 0.1 * pi / (fs/2);
Ws = 0.3 * pi / (fs/2);
Rp = 1;
Rs = 10;
[n,Wc] = buttord(Wp,Ws,Rp,Rs);  %计算滤波器的阶数n和3dB截止频率
%脉冲响应不变法
[B,A] = butter(n,Wc,'s');  %计算相应的模拟滤波器的系统函数
[Bz,Az] = impinvar(B,A);   %用脉冲响应不变法将模拟滤波器转换为数字滤波器H(s)>>H(z)
figure(1)
freqz(Bz,Az,1024,fs)
title('脉冲响应不变法')
%双线法
[A,B] = butter(n,Wc);
% sos = zp2sos(z,p,G);
figure(2)
freqz(A,B,1024,fs)
title('双线法')
