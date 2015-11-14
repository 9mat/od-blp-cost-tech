function [ p, margin, s, iter, flag, distance ] = contraction_bertrand( theta, delta, c, Data, p0 )
%CONTRACTION_BERTRAMD Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta);

alphai = bsxfun(@rdivide, params.alpha*exp(params.sigmap*Data.vprice), Data.income09);

convergence = false;
toler = 1e-6;
iter = 0;
maxiter = 5000;

p = p0;

stepsize = @(r,v) -norm(r)/norm(v);

while ~convergence
    Data.price = p;
    mu = calmu(theta, Data);
    s = calshare(delta, exp(mu), Data.iT);
    margin = calmargin(s, alphai, Data.iF);
    p2 = c + margin;
    
    r = p2 - p;
    
    Data.price = p2;
    mu = calmu(theta, Data);
    s = calshare(delta, exp(mu), Data.iT);
    [margin, flag] = calmargin(s, alphai, Data.iF);
    p3 = c + margin;
    
    v = (p3 - p2) - r;
    a = stepsize(r,v);
    
    p = p - 2*a*r + a^2*v;
    
    distance = max(abs(r));
    convergence = distance < toler;
    iter = iter + 1;
        
    if any(isnan(p)) || (~flag)
        p = p0.*(1 + (rand(size(p))-0.5)*0.2);
    end
    
    if iter > maxiter
        p = NaN(size(p));
        flag = -1;
        return;
    end
end

flag = 0;

end

