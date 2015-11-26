function cf_std
%CF_STD Summary of this function goes here
%   Detailed explanation goes here

settings = loadSettings;
datafile = settings.result_file;
carstd = settings.carstd;
truckstd = settings.truckstd;


if ischar(carstd); carstd = str2double(carstd); end;
if ischar(truckstd); truckstd = str2double(truckstd); end;

load(datafile);

diaryname = ['diary-cf-std-' runid '.txt'];


car = 1-suv-truck-van-minivan;
Data.cafestd(car==1) = carstd;
Data.cafestd(car==0) = truckstd;
Data.cagpmstd = 1./Data.cafestd*100;

coef = -eta(end-3:end);
[gpm1, ps1, gammaj1, cce1, share1] = contraction_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce, cce, coef, ps, gammaj0, diaryname);

resultfile = ['cf-std-' runid '.mat'];
save(resultfile);

diary off;

end

