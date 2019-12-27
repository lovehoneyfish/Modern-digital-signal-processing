%作图学习曲线
figure(5)
plot(1:512,J1,'r');
hold on;
plot(1:512,J2,'b');
hold off;
xlabel('迭代次数');
ylabel('均方误差');
title('100次 步长 0.05红色 0.005蓝色(学习曲线)');
grid on;

