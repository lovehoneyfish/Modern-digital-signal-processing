clear
clc

%产生AR模型的输出信号
a1=0.99;
seta=0.995;
N=1024;

for i=1:600
y=randn(1,N)*sqrt(seta);  % noise sig 
num=1;  %分子系数
den=[1 a1]; %分母系数
u0=zeros(length(den),1)';  %初始数据
xic=filtic(num,den,u0);  %初始条件
un=filter(num,den,y,xic);   %产生数据

%产生期望响应信号和观测数据矩阵
n0=1;  %需要实现n0步线性预测
m=2;   %滤波器阶数 抽头数
b=un(n0+1:N);  %预测期望响应d(n)
l=length(b);
un1 = [zeros(m-1,1), un];  %扩展数据  %%%%
a=zeros(m,l);      %u(n)
for k=1:l-1
    a(:,k)=un1(m-1+k:-1:k);  %观测数据矩阵A矩阵
end

%应用RLS算法迭代求最有权向量
delta=0.004;    %调整参数
lambda=0.98;       %遗忘因子
w=zeros(m,l+1);     %存储权向量
epsilon=zeros(l,1);  %先验估计误差存储
p1=eye(m)/delta;   %相关矩阵的逆
for k=1:l
    pin=p1*a(:,k);
    denok=1+a(:,k)'*pin;      
    kn=pin/denok;
    epsilon(k)=b(k)-w(:,k)'*a(:,k);
    w(:,k+1) = w(:,k) + kn * conj(epsilon(k));
    p1=p1/lambda-kn*a(:,k)'*p1/lambda;
end
mse=abs(epsilon).^2;
JFWC(i,:) =mse;
W1(i,:) = w(1,:);
W2(i,:) = w(2,:);
end


%求800次均值
mse = mean(JFWC);
w1 = mean(W1);     %两个抽头的权值
w2 = mean(W2);

figure(1)
plot(1:800,mse(1:800),'r');
xlabel('迭代次数');
ylabel('MSE');
title('均方误差');
grid on;
%权值
figure(2)
plot(1:800,w1(1:800),'r');
hold on;
plot(1:800,w2(1:800),'b');
hold off;
xlabel('迭代次数');
ylabel('权值');
title('计算权值');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMS算法 运行节
% %产生样本序列
% clear
% clc
% %产生AR模型的输出信号
% a1=0.99;
% seta=0.995;
% data_length=2000;
% h1 = 0.05;
% for i=1:1000
% y=randn(1,data_length)*sqrt(seta);  % noise sig 
% num=1;  %分子系数
% den=[1 a1]; %分母系数
% u0=zeros(length(den),1)';  %初始数据
% xic=filtic(num,den,u0);  %初始条件
% u=filter(num,den,y,xic);   %产生数据
% end
% 
% %%
% %100次试验
% W1 = zeros(100,4500);
% W2 = zeros(100,4500);
% for m = 1:100
% %产生样本序列
% L = 5000;
% a1 = 0.99;
% s = 0.995;
% n = 1:L;
% v = sqrt(s) * randn(L,1);
%     u(1) = v(1);
%     for k = 2:L
%       u(k) = -a1 * u(k-1) + v(k);
%     end
% u=u(500:end);
% 
% %LMS迭代算法
% M = 2;
% w(1,:) = zeros(1,M);
% e(1) = u(1);
% mu = 0.001;
% uu = zeros(1,M);
% w(2,:) = w(1,:) + mu * e(1) * uu;
% uu = [u(1),uu(1:M-1)];
% dd = (w(2,:)*uu')';
% e(2) = u(3) - dd;
% 
% for k = 3:4501
%   w(k,:) = w(k-1,:) + mu * e(k-1) * uu;
%   uu = [u(k-1),uu(1:M-1)];
%   dd = (w(k,:)*uu')';
%   e(k) = u(k)-dd;
% end
% W1(m,:) = w(1:4500,1)';
% W2(m,:) = w(1:4500,2)';
% end
% 
% W11 = mean(W1);
% W22 = mean(W2);
% 
% plot(1:4500,W11,'r');
% hold on;
% plot(1:4500,W22,'b');
% xlabel('迭代次数');
% ylabel('抽头权值');
% title('步长0.05');
% grid on;
