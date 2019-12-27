clear
clc
%产生零均值、方差为1的复高斯白噪声序列 
N = 256;
noise = (randn(1,N) + 1i * randn(1,N) / sqrt(2));
%产生三个复正弦信号
f1 = 0.15;
f2 = 0.17;
f3 = 0.26;           %信号的归一化频率
SNR1 = 30;
SNR2 = 30;
SNR3 = 27;          %信号的信噪比
A1 = 10^(SNR1 / 20);
A2 = 10^(SNR2 / 20);
A3 = 10^(SNR3 / 20);    %信号的幅度
signal1 = A1 * exp(1i * 2 * pi * f1 * (0:N-1));
signal2 = A2 * exp(1i * 2 * pi * f2 * (0:N-1));
signal3 = A3 * exp(1i * 2 * pi * f3 * (0:N-1));     %产生复正弦信号
%产生观察样本u(n)
un = signal1 + signal2 + signal3 + noise;
%周期图法
NF = 1024;                  %周期图法中FFT的点数
Spr = fftshift((1/NF) * abs(fft(un,NF)).^2);
A = 10 * log10(Spr);
f = (-length(A)/2 + 1):(length(A)/2);
figure(1)
plot(f/NF,A);
xlabel('w/2\pi');
ylabel('归一化功率谱/dB');
title('周期图法');

%BT法
M = 64;                        %自相关函数的单边长度
r = xcorr(un,M,'biased');           %计算自相关函数
NF = 1024;                      %BT法中FFT的点数
BT = fftshift(fft(r,NF));             %BT法计算功率谱
B = 10 * log10(BT);
f = (-length(B)/2 + 1):(length(B)/2);
figure(2)
plot(f/NF,B);
xlabel('w/2\pi');
ylabel('归一化功率谱/dB');
title('BT法');
