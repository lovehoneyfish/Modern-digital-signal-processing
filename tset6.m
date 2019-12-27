clear
clc

%产生样本序列

a1 = -0.975;
a2 = 0.95;
s = 0.0731;
trials = 100;        %随机试验次数
data_length = 512;
n = 1:data_length;
v = sqrt(s) * randn(data_length,1);
u0 = [0 0 0];
num = 1;
den = [1 a1 a2];
Zi = filtic(num,den,u0);
u = filter(num,den,v,Zi);
figure(1)
plot(u);
xlabel('n'); ylabel('(n)--y(n)');
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LMS迭代算法

h1 = 0.05;              %步长因子h1  h2
h2 = 0.005;
w1 = zeros(2,data_length);    %在不同步长因子下的w1 w2权向量初始值、存储空间
w2 = zeros(2,data_length);
e1 = zeros(data_length,1);     %估计误差的初始值、存储空间大小
e2 = zeros(data_length,1);
d1 = zeros(data_length,1);
d2 = zeros(data_length,1);      %期望信号的估计同上
for n = 3:data_length-1         %LMS迭代
     w1(:,n+1) = w1(:,n) + h1 * u(n-1:-1:n-2) * conj(e1(n));
     w2(:,n+1) = w2(:,n) + h2 * u(n-1:-1:n-2) * conj(e2(n));
     d1(n+1) = w1(:,n+1)' * u(n:-1:n-1);
     d2(n+1) = w2(:,n+1)' * u(n:-1:n-1);
     e1(n+1) = u(n+1) - d1(n+1);
     e2(n+1) = u(n+1) - d2(n+1);
end
figure(2)
plot(1:512,w1(1,:),'r');
hold on;
plot(1:512,w1(2,:),'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('步长0.05 ');
figure(3)
plot(1:512,w2(1,:),'r');
hold on;
plot(1:512,w2(2,:),'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('步长0.005');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算剩余均方误差、失调参数并画出学习曲线

wopt = zeros(2,trials);     
Jmin = zeros(1,trials);   %最小均方误差初值存储
sum_eig = zeros(trials,1);
%通过维纳霍夫方程计算最小均方误差Jmin
rm(2)=0;
rm(1)=0;
rm(3)=0;
for m=3:trials
    rm(m+1)=xcorr(u(m),'biased');   %计算自相关函数 bisaed 有偏估计
    R=[rm(m) rm(m+1);rm(m-1) rm(m)];
    p=[rm(m-1);rm(m-2)];
%     wopt(m)=R\p;      %单次最佳权值
    [v,d]=eig(R);       %求矩阵R的全部特征值，构成对角阵D，并求R的特征向量构成V的列向量。
    Jmin(m)=rm(m)-p'*wopt(:,m);       %单次维纳误差
    sum_eig(m)=d(1,1)+d(2,2);           %单次特征值之和
end
sJmin = sum(Jmin) / trials;     %100次平均误差
sum_eig_100trials = sum(sum_eig)/100;   %100次平均特征值之和
Jexfin1 = h1 * sJmin * (sum_eig_100trials / (2 - h1 * sum_eig_100trials)); %剩余均方误差稳态取值
Jexfin2 = h2 * sJmin * (sum_eig_100trials /(2 - h2 * sum_eig_100trials));
M1 = Jexfin1 / sJmin;   %计算失调参数
M2 = Jexfin2 / sJmin;

Jex1=Jexfin1;    %剩余均方误差1 2   
Jex2=Jexfin2;


