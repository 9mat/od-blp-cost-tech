function cf_std( datafile, carstd, truckstd )
%CF_STD Summary of this function goes here
%   Detailed explanation goes here

diary
if ischar(carstd); carstd = str2double(carsrd); end;
if ischar(truckstd); truckstd = str2double(truckstd); end;

load(datafile);

diary(['diary-cf-std-' runid '.txt']); diary on;


car = 1-suv-truck-van-minivan;
Data.cafestd(car==1) = carstd;
Data.cafestd(car==0) = truckstd;
Data.cagpmstd = 1./Data.cafestd*100;

coef = -eta(end-3:end);
[gpm1, ps1, gammaj1, share1] = contraction_tech(theta, deltas(:,1:1), cs(:,1:1), Data, cce, coef, ps, gammaj0);

resultfile = ['cf-std-' runid '.mat'];
save(resultfile);

diary off;

end

