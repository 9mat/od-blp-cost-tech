function [ price, e] = pricing( theta, delta, c, e0, Data, tech )
%PRICING Summary of this function goes here
%   Detailed explanation goes here

alpha = theta(end-5);
lambda = theta(end-4);
sigmap = theta(end-3);
sigmae = theta(end-2);
a = theta(end-1);
b = theta(end);

iF = Data.iF;
price = zeros(size(iF));
gpm = zeros(size(iF));
alphai = alpha*exp(sigmap*Data.vprice);

toler = 1e-6;
converged2 = false;
eb = Data.gpm + e0;
e = e0;
e00 = e0;
price = Data.price;
while ~converged2
    e = e00;
    
    converged = false;
    while ~converged
        s = calshare(delta, theta, Data);
        
        margin1 = zeros(size(c));
        lambdai = lambda*exp(sigmae*Data.ve);
        for f=1:max(iF)
            index = iF == f;
            si = s(index, :);
            sv = si.*lambdai(index,:);
            Delta = (diag(sum(sv,2)) - sv*si')/size(s,2);
            margin1(index) = -Delta\((a + 2*b*e(index)).*mean(si,2)./Data.pgreal(index));
        end
        
        
        margin2 = price - c - (a*e+b*e.^2);
        
        e = e + (margin2 - margin1)/(2*b);
        
        distance2 = max(abs(e(:) - e0(:)));
        fprintf('     distance e = %f\n', distance2);
        converged = (distance2 < toler);
        
        Data.gpm = eb - e;
        Data.dpm = gpm.*Data.pgreal;
        e0 = e;
    end
    
    ce = a*e + b*e.^2;
    s = calshare(delta, theta, Data);
    
    for f=1:max(iF)
        index = iF == f;
        si = s(index, :);
        sv = si.*alphai(index,:);
        Delta  = (diag(sum(sv,2)) - sv*si')/size(s,2);
        price(index) = c(index) + ce(index) -  Delta\mean(si,2);
    end
    
    distance1 = max(abs(price(:) - Data.price(:)));
    fprintf('outer distance                = %f\n', distance1);
    converged2 = distance1 < toler;
    Data.price = price;
    e00 = e;
end
end

