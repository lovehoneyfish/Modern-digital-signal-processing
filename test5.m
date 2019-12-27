clear
clc

%产生零均值、方差为1的复高斯白噪声序列v(n)
N = 1000;
noise = (randn(1,N) + 1i * randn(1,N)) *sqrt(2);

%产生带噪声的信号样本u(n)
signal1 = exp(1i * 0.5 * pi * (0:N-1) + 1i * 2 * pi * rand(1));  %产生第一个信号
signal2 = exp(1i * -0.3 * pi * (0:N-1) + 1i * 2 * pi * rand(1));  %产生第二个信号
un = signal1 + signal2 + noise;                %产生带噪声的信号

%计算自相关矩阵
M = 8;                      %自相关矩阵阶数
for k = 1:N-M
    xs( :,k) = un(k + M - 1: -1 : k);    %构造样本矩阵
end
R = xs * xs' / (N-M);            %计算自相关矩阵

%自相关矩阵特征值分解
[U,E] = svd(R);
ev = diag(E);            %提取对角元素上的特征值

%根据AIC准则进行信号源个数估计
for k = 1:M
    dec = prod(ev(k:M) .^ (1/(M - k + 1)));
    nec = mean(ev(k:M));
    lnv = (dec/nec)^((M - k + 1) * N);
    AIC(k) = -2 * log(lnv) + 2 * (k-1) * (2 * M - k + 1);
end
[Amin,K] = min(AIC);   %寻找使AIC准则最小的索引
% N1 = K - 1;           %寻找信号源个数

%计算MUSIC谱
En = U(:,K:M);   %噪声子空间的向量组成的矩阵
NF = 2048;           %MUSIC的扫描点数
for n = 1:NF
      Aq = exp(-1i * 2 * pi * (n-1) / NF * (0:M-1)');
      Pmusic(n) = 1/(Aq' * En * En' * Aq);         %MUSIC谱
end
S = 10 * log10(Pmusic);
m1 = -1023:1024;
plot(m1/2048,S);
xlabel('w/2*pi');
ylabel('归一化功率谱/dB');
title('MUSIC谱估计');
