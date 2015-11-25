result_file = 'result-151124-234322.mat';

cf_pgreal(result_file, 2.74);
cf_std(result_file, 27.5*0.7, 14.49);
cf_pgstd(result_file, 2.74, 27.5*0.7, 14.49);

load(result_file);

cce0 = cce;
cce0(cce0>trim95) = trim95;
cce0(cce0<trim05) = trim05;

coef = -eta(end-3:end);
pow = numel(coef)-1;
f = @(x) bsxfun(@power, x, 0:pow)*coef;
int_f = @(x) bsxfun(@power, x, 1:pow+1)*(coef./(1:pow+1)');
f0 = f(cce0);
int_f0 = int_f(cce0);

car = 1 - truck - suv - van - minivan;

%%
cf_pg_carcafe = zeros(max(cdid), 3);
cf_pg_carmpg = zeros(max(cdid), 3);
cf_pg_truckcafe = zeros(max(cdid), 3);
cf_pg_truckmpg = zeros(max(cdid), 3);

settings = loadSettings;

for i=1:3
    if i==1; load(['cf-pg-' runid '.mat']);
    elseif i==2; load(['cf-std-' runid '.mat']);
    elseif i==3; load(['cf-pgstd-' runid '.mat']);
    end
        
    delta1 = bsxfun(@plus, Data.Xv*beta_v, xi);
    change_c = (f(cce1).*cce1 - f0.*cce0 - int_f(cce1) + int_f0)*10;
    change_c(isnan(change_c)) = 0;
    c1 = bsxfun(@plus, c, change_c);
    
    [p, margin, s, gammaj] = contraction_cafe(theta, delta1, c1, Data, gammaj1, ps1, settings);
    share1 = mean(s,2);
    
    for t=1:max(cdid)
        index = (cdid == t) & (car == 1);
        cf_pg_carcafe(t,i) = sum(share(index))/sum(share(index).*gpm1(index));
        cf_pg_carmpg(t,i) = mean(1./gpm1(index)*100);
        
        index = (cdid == t) & (car == 0);
        cf_pg_truckcafe(t,i) = sum(share(index))/sum(share(index).*gpm1(index));
        cf_pg_truckmpg(t,i) = mean(1./gpm1(index)*100);
    end
end

%%
time = 1:max(cdid);
figure; plot(time, cf_pg_carcafe);
figure; plot(time, cf_pg_carmpg);
