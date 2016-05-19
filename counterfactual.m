function counterfactual(trance, cf_type)

%CF_PGSTD Summary of this function goes here
%   Detailed explanation goes here

if ~ischar(trance); trance = num2str(trance); end

settings = loadSettings;
newpg = settings.newpg;
carstd = settings.carstd;
truckstd = settings.truckstd;

if ischar(newpg); newpg = str2double(newpg); end;
if ischar(carstd); carstd = str2double(carstd); end;
if ischar(truckstd); truckstd = str2double(truckstd); end;

datafile = ['trance-' trance '-' settings.result_file];
load(datafile);

diaryname = ['diary-cf-' cf_type '-trance-' trance '-run-' runid '.txt'];

oldpg = Data.pgreal;
oldstd = Data.cafestd;
oldmampg = mampg;

nsteps = 30;
cce0 = cce;

for step = 1:nsteps
    prop = step/nsteps;
    if ~isempty(strfind(cf_type, 'ma'))
        mampg = prop*mampg06 + (1-prop)*oldmampg;
    end
    
    if ~isempty(strfind(cf_type, 'pg'))
        Data.pgreal = prop*newpg + (1-prop)*oldpg;
    end
    
    if ~isempty(strfind(cf_type, 'ma')) || ~isempty(strfind(cf_type, 'pg'))
        madpm = 1./mampg*100.*Data.pgreal;
        
        Xrc = Data.Xrc;
        Xv = Data.Xv;
        Xrc(:, madpm_idx_rc) = madpm./Data.income09;
        Xv(:, madpm_idx_v) = madpm./Data.income09;
        Data.Xv = Xv;
        Data.Xrc = Xrc;
        Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);
        
        deltas = bsxfun(@plus, Data.Xv*beta_v, xis);
    end
    
    if ~isempty(strfind(cf_type, 'std'))
        Data.cafestd(car==1) = prop*carstd + (1-prop)*oldstd(car==1);
        Data.cafestd(car==0) = prop*truckstd + (1-prop)*oldstd(car==0);
        Data.cagpmstd = 1./Data.cafestd*100;
    end

    coef = -eta(end-3:end);
    [gpm1, ps1, gammaj1, cce1, share1] = contraction_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce0, cce, coef, ps, gammaj0, diaryname);
% [gpm1, ps1, gammaj1, cce1, share1] = broyden_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce, cce, coef, ps, gammaj0);
    cce = cce1;
    gammaj0 = gammaj1;
    ps = ps1;
    
    fprintf(' *** Finish round %d of counterfactual\n', step);
end


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

