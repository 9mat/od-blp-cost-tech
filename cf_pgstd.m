function cf_pgstd( datafile, newpg, carstd, truckstd )
%CF_PGSTD Summary of this function goes here
%   Detailed explanation goes here

if ischar(newpg); newpg = str2double(newpg); end;
if ischar(carstd); carstd = str2double(carstd); end;
if ischar(truckstd); truckstd = str2double(truckstd); end;

load(datafile);

diaryname = ['diary-cf-pgstd-' runid '.txt'];

Data.pgreal = newpg;
madpm = 1./mampg*100*newpg;

Xrc(:, madpm_idx_rc) = madpm./income09;
Xv(:, madpm_idx_v) = madpm./income09;
Data.Xv = Xv;
Data.Xrc = Xrc;
Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);

deltas = bsxfun(@plus, Data.Xv*beta_v, xis);

car = 1-suv-truck-van-minivan;
Data.cafestd(car==1) = carstd;
Data.cafestd(car==0) = truckstd;
Data.cagpmstd = 1./Data.cafestd*100;

coef = -eta(end-3:end);
[gpm1, ps1, gammaj1, cce1, share1] = contraction_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce, coef, ps, gammaj0, diaryname);

resultfile = ['cf-pgstd-' runid '.mat'];
save(resultfile);

diary off;

end

