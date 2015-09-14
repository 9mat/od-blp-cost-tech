function [ deriv ] = dpi( p, e, delta, theta, Data )
%DPI Summary of this function goes here
%   Detailed explanation goes here

Data.price = p;
s = calshare(delta, theta, Data);

alpha = theta(end-5);
lambda = theta(end-4);
sigmap = theta(end-3);
sigmae = theta(end-2);
a = theta(end-1);
b = theta(end);

margin = p - Data.c - (a*e + b*e.^2);
dcde = a + 2*b*e;

dpidp = zeros(length(p),1);
dpide = zeros(length(e),1);

N = size(s,2);

for f = 1:max(Data.iF);
    index = Data.iF == f;
    si = s(index,:);
    
    pg = mean(Data.pgreal(index));
    alphai = alpha*exp(sigmap*Data.vprice(index,:));
    lambdai = lambda*exp(sigmae*Data.ve(index,:));
    
    dsdp = (diag(sum(si.*alphai,2)) - (si.*alphai)*si')/N;
    dsde = -(diag(sum(si.*lambdai,2)) - (si.*lambdai)*si')*pg/N;
    
    dpidp(index) = mean(si,2) + dsdp*margin(index);
    dpide(index) = -dcde(index).*mean(si,2) + dsde*margin(index);
end

deriv = [dpidp; dpide];

end

