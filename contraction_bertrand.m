function [ p, margin, s, iter, flag, distance ] = contraction_bertrand( theta, delta, c, comply_mc, Data, p0 )
%CONTRACTION_BERTRAMD Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta);

alphai = bsxfun(@rdivide, params.alpha*exp(params.sigmap*Data.vprice), Data.income09);

convergence = false;
toler = 1e-6;
iter = 0;
maxiter = 5000;

p = p0;

while ~convergence
    Data.price = p;
    mu = calmu(theta, Data);
    s = calshare(delta, exp(mu), Data.iT);
    
    margin = calmargin(s, alphai, Data.iF);
    
    p = c + margin - comply_mc;
    distance = max(abs(Data.price - p));
    convergence = distance < toler;
    iter = iter + 1;
    
    if any(isnan(p)) || iter > maxiter
        p = NaN(size(p));
        flag = -1;
        return;
    end
end

flag = 0;

end

