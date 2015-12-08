settings = loadSettings;

datafile = ['trance-' trance '-' settings.result_file];
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
    
    for trance = 1:9
        resultfile = ['cf-' cf_type '-trance-' trance '-run-' runid '.mat'];
        if exist(resultfile, 'file') > 0
            load(resultfile);
            combCarCAFE(trance, cf_code) = carcafe;
            combCarMPG(trance, cf_code) = carmpg;
            combTruckCAFE(trance, cf_code) = truckcafe;
            combTruckMPG(trance, cf_code) = truckmpg;
        end
    end
end
