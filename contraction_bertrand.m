function [ p, margin, s, iter, flag, distance ] = contraction_bertrand( theta, delta, c, Data, gammaj, p0 )
%CONTRACTION_BERTRAMD Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta);

alphai = bsxfun(@rdivide, params.alpha*exp(params.sigmap*Data.vprice), Data.income09);

convergence = false;
toler = 1e-2;
iter = 0;
maxiter = 5000;

p = p0;

stepsize = @(r,v) -norm(r)/norm(v);

[~, ~, fleet] = unique([Data.iF, Data.fleet], 'rows');

cafe0 = Data.cafe;
cafe20 = Data.cafe2;
comply_mc = calcomply_mc_j(gammaj, Data);

while ~convergence
    Data.price = p;
    mu = calmu(theta, Data);
    s = calshare(delta, exp(mu), Data.iT);
    [margin, flag] = calmargin(s, alphai, Data.iF);

%     share = mean(s,2);
%     newcafe2 = accumarray(fleet, share)./accumarray(fleet, share.*(Data.gpm/100));
%     newcafe2 = newcafe2(fleet);
%     Data.cafe = newcafe2.*cafe0./cafe20;
%     Data.cagpm = 1./Data.cafe*100;
    comply_mc1 = comply_mc;
    comply_mc = calcomply_mc_j(gammaj, Data);
%     comply_mc((Data.comply == 0) & (Data.cafe > cafe0)) = 0;
    dd = max(abs(comply_mc - comply_mc1));

    p2 = c + margin - comply_mc;
    
    r = p2 - p;
    
    Data.price = p2;
    mu = calmu(theta, Data);
    s = calshare(delta, exp(mu), Data.iT);
    [margin, flag] = calmargin(s, alphai, Data.iF);
    
    share = mean(s,2);
    newcafe2 = accumarray(fleet, share)./accumarray(fleet, share.*(Data.gpm/100));
    newcafe2 = newcafe2(fleet);
%     Data.cafe = newcafe2.*cafe0./cafe20;
%     Data.cagpm = 1./Data.cafe*100;
    comply_mc = calcomply_mc(params.gamma, Data);

    p3 = c + margin - comply_mc;
    
    v = (p3 - p2) - r;
    a = stepsize(r,v);
    
    p = p - 2*a*r + a^2*v;

%     p = p + r;
    
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

