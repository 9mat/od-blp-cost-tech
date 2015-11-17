function [ cce, ps, share, gammaj ] = contraction_cafe( theta, deltas, cs, Data, gammaj, ps )
%CONTRACTION_CAFE Summary of this function goes here
%   Detailed explanation goes here
cafe2short = Data.cafe2short; % accumarray(fleet, share)./accumarray(fleet, share.*(gpm/100));
fleet = Data.fleet;
cafe2long = cafe2short(fleet);
params = getParams(theta);
gammaf = params.gamma*ones(size(cafe2short));
binding = Data.comply == 0;

for i = 1:100
    [cce, ps, share] = calcce(theta, deltas, cs, Data, gammaj, ps);
    newcafe2 = accumarray(fleet, share)./accumarray(fleet, share.*(Data.gpm/100));
    gammaf = gammaf - 0.1*(newcafe2 - cafe2short);
    gammajj = gammaf(fleet);
    gammaj(binding) = gammajj(binding);
    newcafe2long = newcafe2(fleet);
    gammaj(gammaj < 0) = 0; % the constraints become non-binding  
    index = binding & (gammaj > 0);
    distance = max(abs(newcafe2long(index) - cafe2long(index)));
    fprintf('****** CAFE iteration %d, distance = %f\n', i, distance);
end


end

