%进行100次独立实验
trials = 100;
data_length = 512;
wopt = zeros(2,trials);
Jmin = zeros(1,trials);
sum_eig = zeros(trials,1);
w11 = zeros(trials,data_length);
w12 = zeros(trials,data_length);
w21 = zeros(trials,data_length);
w22 = zeros(trials,data_length);
e1 = zeros(trials,1);
e2 = zeros(trials,1);
%步长
h1 = 0.05;
h2 = 0.005;
%计算最小均方误差
for m = 1:trials
    %产生样本序列
    a1 = -0.975;
    a2 = 0.95;
    s = 0.0731;
    n = 1:data_length;
    v = sqrt(s)*randn(data_length,1);
    u0 = [0 0 0];
    num = 1;
    den = [1 a1 a2];
    Zi = filtic(num,den,u0);
    u = filter(num,den,v,Zi);
   %LMS迭代算法
    w1 = zeros(2,data_length);
    w2 = zeros(2,data_length);
    e1 = zeros(data_length,1);
    e2 = zeros(data_length,1);
    d1 = zeros(data_length,1); 
    d2 = zeros(data_length,1);
    for n = 3:data_length-1
      w1(:,n+1) = w1(:,n) + h1 * u(n-1:-1:n-2) * conj(e1(n));
      w11(m,n+1) = w1(1,n+1);
      w12(m,n+1) = w1(2,n+1);
      w2(:,n+1) = w2(:,n) + h2 * u(n-1:-1:n-2) * conj(e2(n));
      w21(m,n+1) = w2(1,n+1);
      w22(m,n+1) = w2(2,n+1);
      d1(n+1) = w1(:,n+1)' * u(n:-1:n-1);
      d2(n+1) = w2(:,n+1)' * u(n:-1:n-1);
      e1(n+1) = u(n+1) - d1(n+1);
      e2(n+1) = u(n+1) - d2(n+1);
    end
    e1(m) = mean(e1);
    e2(m) = mean(e2);
    rm = xcorr(u,'biased');
    R = [rm(512),rm(513);rm(511),rm(512)];
    rm512(m) = rm(512);
    rm513(m) = rm(513);
    rm511(m) = rm(511);
    p = [rm(511);rm(510)];
    wopt(:,m) = R \ p; 
    [v,d] = eig(R);
    Jmin(m) = rm(512) - p' * wopt(:,m);
    sum_eig(m) = d(1,1) + d(2,2);
end
%100次平均误差
sJmin = sum(Jmin) / trials;
%100次平均特征值之和
sum_eig_100trials = sum(sum_eig)/100;
Jexfin1 = h1 * sJmin * (sum_eig_100trials / (2 - h1 * sum_eig_100trials));
Jexfin2 = h2 * sJmin * (sum_eig_100trials /(2 - h2 * sum_eig_100trials));
%计算失调参数
M1 = Jexfin1 / sJmin;
M2 = Jexfin2 / sJmin;
%计算100次的系数平均
q1 = mean(w11(:,:));
q2 = mean(w12(:,:));
q3 = mean(w21(:,:));
q4 = mean(w22(:,:));
rm11 = mean(rm511);
rm22 = mean(rm512);
rm33 = mean(rm513);
figure(3)
subplot(2,1,1)
plot(1:512,q1,'r');
hold on;
plot(1:512,q2,'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('100次步长0.05');

subplot(2,1,2)
plot(1:512,q3,'r');
hold on;
plot(1:512,q4,'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('100次 步长0.005');
w1=mean(wopt(1,:));
w2=mean(wopt(2,:));
%计算均方误差
for t=1:data_length 
    J1(t)=sJmin+([q1(t) q2(t)]-[w1 w2])*[rm22,rm33;rm11,rm22]*([q1(t) q2(t)]-[w1 w2])';
    J2(t)=sJmin+([q3(t) q4(t)]-[w1 w2])*[rm22,rm33;rm11,rm22]*([q1(t) q2(t)]-[w1 w2])';
end
e1_100trials_ave=J1;
e2_100trials_ave=J2;
%计算剩余均方误差
Jex1=e1_100trials_ave-sJmin;
Jex2=e2_100trials_ave-sJmin;

figure(4)
plot(1:512,J1,'r');
hold on;
plot(1:512,J2,'b');
hold off;
xlabel('迭代次数');
ylabel('均方误差');
title('100次 步长 0.05红色 0.005蓝色');
grid on;
