settings = loadSettings;

datafile = settings.result_file;
load(datafile);

combCarCAFE = NaN(9,4);
combCarMPG = NaN(9,4);
combTruckCAFE = NaN(9,4);
combTruckMPG = NaN(9,4);

cf_types = {'pg', 'std', 'pgstd', 'original'};

for cf_code = 1:4
    cf_type = cf_types{cf_code};
    
    if strcmp(cf_type, 'original') == 1
        share = data(:,4);
        mpg = data(:,9);
        gpm = 100./mpg;
        car = 1 - suv - van - minivan - truck;

        for tranceid = 1:9
            index = (cdid == tranceid) & (car == 1);
            combCarCAFE(tranceid, cf_code) = sum(share(index))/sum(share(index).*gpm(index));
            combCarMPG(tranceid, cf_code) = mean(100./gpm(index));

            index = (cdid == tranceid) & (car == 0);    
            combTruckCAFE(tranceid, cf_code) = sum(share(index))/sum(share(index).*gpm(index));
            combTruckMPG(tranceid, cf_code) = mean(100./gpm(index));
        end
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

