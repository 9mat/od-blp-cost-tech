function [ gpm, ps, gammaj, share ] = contraction_tech(theta, deltas, cs0, Data, cce0, coef, ps, gammaj)
%CONTRACTION_TECH Summary of this function goes here
%   Detailed explanation goes here

convergence = false;
toler = 1e-6;

gpm = Data.gpm;

trim05 = prctile(cce0, 3);
trim95 = prctile(cce0, 97);
cce0(cce0>trim95) = trim95;
cce0(cce0<trim05) = trim05;
cce = cce0;

% make sure that coef is for log(gpm) equation, not log(mpg)
pow = numel(coef)-1;
f = @(x) bsxfun(@power, x, 0:pow)*coef;
int_f = @(x) bsxfun(@power, x, 1:pow+1)*(coef./(1:pow+1)');
f0 = f(cce0);
int_f0 = int_f(cce0);

iter = 0;
stepsize = @(r,v) -norm(r)/norm(v);
gpm0 = Data.gpm;

step = 0.1;
while ~convergence
    change_c = (f(cce).*cce - f0.*cce0 - int_f(cce) + int_f0)*10;
    change_c(isnan(change_c)) = 0;
    cs = bsxfun(@plus, cs0, change_c);
    
    Data.gpm = gpm;
    Data.dpm = gpm.*Data.pgreal;

    nnan = sum(any(isnan(ps)));
    ps(:, any(isnan(ps))) = repmat(Data.price, [1 nnan]);
    [cce2, ps, gammaj, share] = calcce(theta, deltas, cs, Data, gammaj, ps);
    
    cce2 = cce + step*(cce2-cce);
    r = cce2 - cce;
    
    gpm2 = gpm0.*exp(f(cce2)-f0);
    gpm2(isnan(gpm2)) = Data.gpm(isnan(gpm2));
    
    
    change_c = (f(cce2).*cce2 - f0.*cce0 - int_f(cce2) + int_f0)*10;
    change_c(isnan(change_c)) = 0;
    cs = bsxfun(@plus, cs0, change_c);

    Data.gpm = gpm2;
    Data.dpm = gpm2.*Data.pgreal;
    [cce3, ps, gammaj] = calcce(theta, deltas, cs, Data, gammaj, ps);
    
    cce3 = cce3 + step*(cce3-cce2);
    v = (cce3 - cce2) - r;
    
    
    a = stepsize(r,v);
    cce = cce - 2*a*r + a^2*v;

    gpm = gpm0.*exp(f(cce) - f0);
%     gpm = gpm.*exp(0.1*r);
    
    distance = max(abs(r/step));
    convergence = distance < toler;
    iter = iter + 1;
    fprintf(' *** tech contraction iter #%d, distance = %f\n', iter, distance);
end

end

