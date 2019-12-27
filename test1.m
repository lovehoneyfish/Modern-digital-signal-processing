%巴特沃夫滤波器方法
% clear
% clc
% fs=12000;
% N=12000;
% n=0:N-1;
% t = n / fs;  %时间序列
% 
% s = 15 * rand(size(t));  %原始音频信号
% w = 20 * rand(1) * sin(2 * pi * 1500 * t);  %噪声信号
% x = s + w;  %混合信号
% 
% f = fs * (0:(N / 2)) / N;  %频率序列
% 
% %原始信号的频谱分析
% y = fft(s);
% P1 = abs(y/N);
% mag = P1(1:N/2+1);
% % mag(2:end-1) = 2 * mag(2:end-1);
% 
% %混合信号的频谱分析
% y1 = fft(x);
% P2 = abs(y1/N);
% mag1 = P2(1:N/2+1);
% mag1(2:end-1) = 2 * mag1(2:end-1);
% 
% %(带阻)滤波器参数
% Wp = [1400,1600]/(fs/2);
% Ws = [1450,1550]/(fs/2);
% Rp = 0.1;
% Rs = 20;
% 
% %巴特沃思滤波器的设计
% [aa,Wc] = buttord(Wp,Ws,Rp,Rs);
% %[z,p,G] = butter(n,Wc,'stop');
% %figure(1)
% %SOS = zp2sos(z,p,G);       双线性法
% %freqz(SOS,1024,fs);
% [B,A] = butter(aa,Wc,'stop');
% Y = filter(B,A,x);
% 
% %滤波后的频谱分析
% y2 = fft(Y);
% P3 = abs(y2/N);
% mag2 = P3(1:N/2+1);
% % mag2(2:end-1) = 2 * mag2(2:end-1);
% 
% figure(2)
% subplot(2,1,1)
% plot(f,mag,'b')
% ylim([0 0.5]);xlabel('频率/Hz');ylabel('振幅')
% title('原始信号幅频图')
% subplot(2,1,2)
% plot(f,mag1,'b')
% ylim([0 1]);xlabel('频率/Hz');ylabel('振幅')
% title('加入噪声后的幅频图')
% figure(3)
% subplot(2,1,1)
% plot(f,mag,'b')
% ylim([0 0.5]);xlabel('频率/Hz');ylabel('振幅')
% title('原始信号幅频图')
% subplot(2,1,2)
% plot(f,mag2,'b')
% ylim([0 0.5]);xlabel('频率/Hz');ylabel('振幅')
% title('滤波后的幅频图')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%陷波滤波器方法
clear
clc

fs = 12000;  %采样频率
N = 12000;   %采样点数
n = 0:N-1;
t = n / fs;  %时间序列

s = 10000 * rand(size(t));  %原始音频信号 
w = 10000 * rand(1) * sin(2 * pi * 1500 * t);  %噪声信号
x = s + w;  %混合信号

f = fs * (0:(N / 2)) / N;  %0到6khz频率序列

%混合信号的频谱分析
y1 = fft(x);
P2 = abs(y1/N);
mag1 = P2(1:N/2+1);
mag1(2:end-1) = 2 * mag1(2:end-1);

%陷波器设计
r=0.9;
w0=2*pi*1500/fs;   %要滤掉的频率
b=[1 -2*cos(w0) 1];
a=[1 -2*r*cos(w0) r*r];
[H,w]=freqz(b,a,N);
figure(1)
subplot(221);plot(w,abs(H));title('陷波器的幅频响应');
subplot(222);plot(w,angle(H));title('陷波器的相频响应');
subplot(223);zplane(b,a);title('陷波器的零极点图');

%陷波器滤波
y2=filter(b,a,x);
Y=fft(y2);
P3=abs(Y/N);
mag2=P3(1:N/2+1);
mag2(2:end-1) = 2 * mag2(2:end-1);

figure(2);
subplot(2,1,1)
plot(f,mag1)
ylim([0 500]);xlabel('频率/Hz');ylabel('振幅')
title('原始信号频谱图')
subplot(2,1,2)
plot(f,mag2)
ylim([0 500]);xlabel('频率/Hz');ylabel('振幅')
title('陷波器滤波后的信号y(n)')



