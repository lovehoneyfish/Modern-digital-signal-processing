clear
clc
wp = 0.4 * pi;
ws = 0.096 * pi;
rs = 65;
DB = abs(wp - ws);           %计算过渡带宽度
M = ceil(12*pi/ DB);  %计算布莱克曼窗所需阶数 
N=M+mod(M+1,2);   %确保hn长度N为奇数
wc = (wp + ws)/ 2 / pi;    %计算理想低通滤波器通带截止频率
hn = fir1(N-1,wc,blackman(N));  %调用fir1计算低通FIRDF的h(n)


%波形图
figure(1)
stem(hn);
xlabel('n');
ylabel('h(n)');
title('h(n)波形');
% %损耗函数曲线
% [Hn,w] = freqz(hn);        %计算频率响应
% db = 20*log10(abs(Hn));    %化为分贝值
% db1=db';
% figure(2)
% plot(w / pi,db1);
% xlabel('w/ \p');
% ylabel('20lg|Hg(w)|');
% title('损耗函数曲线')
% %相频特性曲线
% figure(3)
% plot(w,angle(Hn));
% xlabel('频率');
% ylabel('相位');
% title('相频特性曲线')
figure(4)
freqz(hn,1,512)

