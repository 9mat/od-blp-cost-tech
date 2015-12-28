function counterfactual_laststep(trance, cf_type)

%CF_PGSTD Summary of this function goes here
%   Detailed explanation goes here

settings = loadSettings;
datafile = ['trance-' trance '-' settings.result_file];
load(datafile);

if ~ischar(trance); trance = num2str(trance); end

resultfile = ['cf-' cf_type '-trance-' trance '-run-' runid '.mat'];
load(resultfile);

pow = numel(coef)-1;
f = @(x) bsxfun(@power, x, 0:pow)*coef;
int_f = @(x) bsxfun(@power, x, 1:pow+1)*(coef./(1:pow+1)');
f0 = f(cce);
int_f0 = int_f(cce);

Data.gpm = gpm1;
Data.dpm = gpm1.*Data.pgreal;
delta1 = bsxfun(@plus, Data.Xv*beta_v, xi);
change_c = -(f(cce1).*cce1 - f0.*cce - int_f(cce1) + int_f0)*10;
change_c(isnan(change_c)) = 0;
c1 = bsxfun(@plus, c, change_c);

[p2, margin2, s2, gammaj2] = contraction_cafe(theta, delta1, c1, Data, gammaj1, ps1, settings);
share2 = mean(s2,2);

index = (car == 1);
carcafe = sum(share2(index))/sum(share2(index).*gpm1(index));
carmpg = mean(1./gpm1(index)*100);

index = (car == 0);
truckcafe = sum(share2(index))/sum(share2(index).*gpm1(index));
truckmpg = mean(1./gpm1(index)*100);

resultfile = ['cf-' cf_type '-trance-' trance '-run-' runid '.mat'];
save(resultfile);

diary off;

end

