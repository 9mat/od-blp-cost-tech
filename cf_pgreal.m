function cf_pgreal
%CF_PGREAL Summary of this function goes here
%   Detailed explanation goes here

if exist('parpool') > 0
    parpool;
end

settings = loadSettings;
datafile = settings.result_file;
newpg = settings.newpg;

if ischar(newpg); newpg = str2double(newpg); end;
disp(newpg)
load(datafile);

diaryname = ['diary-cf-pgreal-' runid '.txt'];

Data.pgreal = newpg;
madpm = 1./mampg*100*newpg;

Xrc(:, madpm_idx_rc) = madpm./income09;
Xv(:, madpm_idx_v) = madpm./income09;
Data.Xv = Xv;
Data.Xrc = Xrc;
Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);

deltas = bsxfun(@plus, Data.Xv*beta_v, xis);

coef = -eta(end-3:end);
[gpm1, ps1, gammaj1, cce1, share1] = contraction_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce, cce, coef, ps, gammaj0, diaryname);

resultfile = ['cf-pg-' runid '.mat'];
save(resultfile);

diary off;
end

