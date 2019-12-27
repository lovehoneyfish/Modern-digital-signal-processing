clear
clc
%产生零均值、方差为1的复高斯白噪声序列 
N = 256;
noise = (randn(1,N) + 1i * randn(1,N) / sqrt(2));
%产生三个复正弦信号
f1 = 0.15;
f2 = 0.17;
f3 = 0.26;              %信号的归一化频率
SNR1 = 30;
SNR2 = 30;
SNR3 = 27;             %信号的信噪比
A1 = 10^(SNR1 / 20);
A2 = 10^(SNR2 / 20);
A3 = 10^(SNR3 / 20);         %信号的幅度
signal1 = A1 * exp(1i * 2 * pi * f1 * (0:N-1));
signal2 = A2 * exp(1i * 2 * pi * f2 * (0:N-1));
signal3 = A3 * exp(1i * 2 * pi * f3 * (0:N-1));     %产生复正弦信号
%产生观察样本u(n)
un = signal1 + signal2 + signal3 + noise;
%计算自相关函数值
p = 16;                            %AR模型的阶数
r0 = xcorr(un,p,'biased');         %直接计算自相关函数
r = r0(p + 1:2 * p + 1);           %提取r(0),r(1),...,r(p)
%计算一阶AR模型的系数与输入方差
a(1,1) = -r(2)/r(1);                       %1阶AR模型的系数
signal(1) = r(1) - (abs(r(2)^2)/r(1));          %1阶AR模型的输入方差
%Levinsion-Durbin迭代算法
for m = 2:p
    k(m) = -(r(m+1) + sum(a(m-1,1:m-1) .* r(m:-1:2)))/signal(m-1);
    a(m,m) = k(m);
    for i = 1:m-1
        a(m,i) = a(m-1,i) + k(m) * conj(a(m-1,m-i));
    end
    signal(m) = signal(m-1) * (1-abs(k(m))^2);
end
%计算十六阶AR模型的功率谱
NF = 1024;                          %AR方法中FFT的点数
Par = signal(p)./fftshift(abs(fft([1,a(p,:)],NF)).^2);  %p阶AR模型的功率谱
C = 10*log10(Par);
f = (-length(C)/2 + 1):(length(C)/2);
figure(1)
plot(f/NF,C);
xlabel('w/(2*pi)');
ylabel('归一化功率谱/dB');
title('16阶AR模型');
