clc
clear
p = rand(100,3);
x = zeros(1000,100);
y = zeros(1000,100);
z = zeros(1000,100);
for i = 1:100
    x(1,i) = p(i,1);
    y(1,i) = p(i,2);
    z(1,i) = p(i,3);
end
for iter = 2:1000
    for i = 1:100
        a = p(i,1);
        b = p(i,2);
        c = p(i,3);
        p(i,1) = 0.00002*a*((1-a)*(9866*b+4259*c-5063*b*c))+a;
        p(i,2) = 0.00002*b*((1-b)*(9866*a+621*c-1818*a*c))+b;
        p(i,3) = 0.00002*c*((1-c)*(4259*a+621*b-4307*a*b))+c;
        x(iter,i) = p(i,1);
        y(iter,i) = p(i,2);
        z(iter,i) = p(i,3);
    end
end