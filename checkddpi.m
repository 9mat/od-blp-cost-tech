epsilon = 1e-5;

s = mean(calshare(delta, theta, Data),2);

%% check deriv wrt prices
ddpidp2 = zeros(size(price));

for j = 1:size(s,1)
    p1 = price;
    p2 = price;
    p1(j) = p1(j) - epsilon;
    p2(j) = p2(j) + epsilon;
    
    Data1 = Data;
    Data2 = Data;
    
    Data1.price = p1;
    Data2.price = p2;
    
    s1 = mean(calshare(delta, theta, Data1),2);
    s2 = mean(calshare(delta, theta, Data2),2);
    
    index = iF == firmid(j);
    pi = sum((price(index) - c(index) - (a*e(index) + b*e(index).^2)).*s(index));
    pi1 = sum((p1(index) - c(index) - (a*e(index) + b*e(index).^2)).*s1(index));
    pi2 = sum((p2(index) - c(index) - (a*e(index) + b*e(index).^2)).*s2(index));
    
    ddpidp2(j) = (pi2 - 2*pi + pi1)/(epsilon^2);
end

%% check deriv wrt e
ddpide2 = zeros(size(price));
s = mean(calshare(delta, theta, Data),2);
for j = 1:size(s,1)
    e1 = e;
    e2 = e;
    e1(j) = e1(j) - epsilon;
    e2(j) = e2(j) + epsilon;
    
    Data1 = Data;
    Data2 = Data;
    
    Data1.dpm(j) = Data1.dpm(j) - Data1.pgreal(j)*(e1(j) - e(j));
    Data2.dpm(j) = Data2.dpm(j) - Data2.pgreal(j)*(e2(j) - e(j));
    
    s1 = mean(calshare(delta, theta, Data1),2);
    s2 = mean(calshare(delta, theta, Data2),2);
    
    index = iF == iF(j);
    pi  = sum((price(index) - c(index) - (a*e(index) + b*e(index).^2)).*s(index));
    pi1 = sum((price(index) - c(index) - (a*e1(index) + b*e1(index).^2)).*s1(index));
    pi2 = sum((price(index) - c(index) - (a*e2(index) + b*e2(index).^2)).*s2(index));
    
    ddpide2(j) = (pi2 - 2*pi + pi1)/(epsilon^2);
    fprintf('%d %f\n', j, ddpide2(j));
end


%%

% check cross deriv
k = 5;
r = 3;
assert(iF(k) == iF(r));

p1 = price;
p2 = price;
p3 = price;
p4 = price;

p1(k) = p1(k) + epsilon;
p1(r) = p1(r) + epsilon;

p2(k) = p2(k) - epsilon;
p2(r) = p2(r) + epsilon;

p3(k) = p3(k) + epsilon;
p3(r) = p3(r) - epsilon;

p4(k) = p4(k) - epsilon;
p4(r) = p4(r) - epsilon;


Data1 = Data;
Data2 = Data;
Data3 = Data;
Data4 = Data;

Data1.price = p1;
Data2.price = p2;
Data3.price = p3;
Data4.price = p4;

s1 = mean(calshare(delta, theta, Data1),2);
s2 = mean(calshare(delta, theta, Data2),2);
s3 = mean(calshare(delta, theta, Data3),2);
s4 = mean(calshare(delta, theta, Data4),2);

index = iF == iF(k);
pi1 = sum((p1(index) - c(index) - (a*e(index) + b*e(index).^2)).*s1(index));
pi2 = sum((p2(index) - c(index) - (a*e(index) + b*e(index).^2)).*s2(index));
pi3 = sum((p3(index) - c(index) - (a*e(index) + b*e(index).^2)).*s3(index));
pi4 = sum((p4(index) - c(index) - (a*e(index) + b*e(index).^2)).*s4(index));

ddpidpdp = (pi1+pi4-pi2-pi3)/(4*epsilon^2);

%%

[d, dd] = ddpi(price, e, delta, theta, Data);
ddp = zeros(size(dd));

for j=1:length(price)
    p1 = price;
    p2 = price;
    
    p1(j) = p1(j) - epsilon;
    p2(j) = p2(j) + epsilon;
    
    d1 = ddpi(p1, e, delta, theta, Data);
    d2 = ddpi(p2, e, delta, theta, Data);
    ddp(:,j) = (d2-d1)/(2*epsilon);
end

%
dde = zeros(size(dd));

for j=1:length(price)
    e1 = e;
    e2 = e;
    
    e1(j) = e1(j) - epsilon;
    e2(j) = e2(j) + epsilon;
    
    Data1 = Data;
    Data2 = Data;
    
    Data1.dpm(j) = Data1.dpm(j) - Data1.pgreal(j)*(e1(j) - e(j));
    Data2.dpm(j) = Data2.dpm(j) - Data2.pgreal(j)*(e2(j) - e(j));
    
    d1 = ddpi(price, e1, delta, theta, Data1);
    d2 = ddpi(price, e2, delta, theta, Data2);
    dde(:,j) = (d2-d1)/(2*epsilon);
end
%
dd2 = [ddp(:,1:length(price)) dde(:,1:length(e))];
