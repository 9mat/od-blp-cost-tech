settings = loadSettings;

datafile = settings.result_file;
load(datafile);

combCarCAFE = NaN(9,4);
combCarMPG = NaN(9,4);
combTruckCAFE = NaN(9,4);
combTruckMPG = NaN(9,4);

for cf_code = 1:3
    if cf_code == 1; cf_type = 'pg';
    elseif cf_code == 2; cf_type = 'std';
    elseif cf_code == 3; cf_type = 'pgstd';
    end
    
    for tranceid = 1:9
        resultfile = ['cf-' cf_type '-trance-' num2str(tranceid) '-run-' runid '.mat'];
        if exist(resultfile, 'file') > 0
            load(resultfile);
            combCarCAFE(tranceid, cf_code) = carcafe;
            combCarMPG(tranceid, cf_code) = carmpg;
            combTruckCAFE(tranceid, cf_code) = truckcafe;
            combTruckMPG(tranceid, cf_code) = truckmpg;
        end
    end
end

share = data(:,4);
mpg = data(:,9);
gpm = 100./mpg;
car = 1 - suv - van - minivan - truck;

for tranceid = 1:9
    index = (cdid == tranceid) & (car == 1);
    combCarCAFE(tranceid, 4) = sum(share(index))/sum(share(index).*gpm(index));
    combCarMPG(tranceid, 4) = mean(100./gpm(index));
    
    index = (cdid == tranceid) & (car == 0);    
    combTruckCAFE(tranceid, 4) = sum(share(index))/sum(share(index).*gpm(index));
    combTruckMPG(tranceid, 4) = mean(100./gpm(index));
end