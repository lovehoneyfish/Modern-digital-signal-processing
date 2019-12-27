clear
clc
%产生零均值、方差为1的复高斯白噪声序列 
N = 32;
noise = (randn(1,N) + 1i * randn(1,N))* sqrt(2);
%产生三个复正弦信号
f1 = 0.15;
f2 = 0.17;
f3 = 0.26;        %信号的归一化频率
SNR1 = 30;
SNR2 = 30;
SNR3 = 27;       %信号的信噪比
A1 = 10^(SNR1 / 20);
A2 = 10^(SNR2 / 20);
A3 = 10^(SNR3 / 20);        %信号的幅度
signal1 = A1 * exp(1i * 2 * pi * f1 * (0:N-1));
signal2 = A2 * exp(1i * 2 * pi * f2 * (0:N-1));
signal3 = A3 * exp(1i * 2 * pi * f3 * (0:N-1));     %产生复正弦信号
%产生观察样本u(n)
un = signal1 + signal2 + signal3 + noise;
%基于FFT的自相关函数快速计算方法
Uk = fft(un,2*N);                 %对un进行2N点的FFT
Sk = (1/N)*abs(Uk).^2;            %计算功率谱估计Sk
r0 = ifft(Sk);                    %对功率谱估计Sk求FFT
r1 = [r0(N+2:2*N),r0(1:N)];       %根据教材（3.1.8）求得自相关函数

figure(1)
stem(real(r1));               %提取实部
xlabel('m');
ylabel('实部');
figure(2)
stem(imag(r1));               %提取虚部
xlabel('m');
ylabel('虚部');
%教材式（3.1.2）估计自相关函数
r = xcorr(un,N-1,'biased');
figure(3)
stem(real(r));                %提取实部
xlabel('m');
ylabel('实部');
figure(4)
stem(imag(r));                %提取虚部
xlabel('m');
ylabel('虚部');

