%LMS迭代算法
h1 = 0.05;
h2 = 0.005;
w1 = zeros(2,data_length);
w2 = zeros(2,data_length);
e1 = zeros(data_length,1);
e2 = zeros(data_length,1);
d1 = zeros(data_length,1);
d2 = zeros(data_length,1);
for n = 3:data_length-1
     w1(:,n+1) = w1(:,n) + h1 * u(n-1:-1:n-2) * conj(e1(n));
     w2(:,n+1) = w2(:,n) + h2 * u(n-1:-1:n-2) * conj(e2(n));
     d1(n+1) = w1(:,n+1)' * u(n:-1:n-1);
     d2(n+1) = w2(:,n+1)' * u(n:-1:n-1);
     e1(n+1) = u(n+1) - d1(n+1);
     e2(n+1) = u(n+1) - d2(n+1);
end
figure(1)
plot(1:512,w1(1,:),'r');
hold on;
plot(1:512,w1(2,:),'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('步长0.05');
figure(2)
plot(1:512,w2(1,:),'r');
hold on;
plot(1:512,w2(2,:),'b');
hold off;
xlabel('迭代次数');
ylabel('抽头权值');
title('步长0.005');
