function [ cce ] = calcce(theta, beta_v, beta_c, xis, omegas, Data )
%CALCCE Summary of this function goes here
%   Detailed explanation goes here

deltas = bsxfun(@plus, Data.Xv*beta_v, xis);
cs = exp(bsxfun(@plus, Data.Xc*beta_c, omegas));

binding = Data.comply == 0;
fined = Data.comply == -1;
params = getParams(theta);

gammai = zeros(size(Data.comply));
gammai(binding) = params.gamma;
gammai(fined) = params.gamma + 0.055;

comply_mc1 = (1-Data.gpm./Data.cagpm).*Data.cafe;
comply_mc2 = (1-Data.cagpm./Data.cagpmstd).*Data.cafestd;
comply_mc2(binding) = 0;
comply_mc = gammai.*(comply_mc1 + comply_mc2);
comply_mc(isnan(comply_mc)) = 0;

lambdai = -bsxfun(@rdivide, params.lambda*bsxfun(@times, exp(params.sigmae*Data.ve), Data.dpm), Data.income09);

ns = min(size(xis,2), size(omegas,2));
J = size(xis,1);
ce = zeros(J,ns);
ss = zeros(J,ns);
ps = zeros(J,ns);

for i=1:ns
    [~, mm, s, iter, flag, distance] = contraction_bertrand(theta, deltas(:,i), cs(:,i), comply_mc, Data, Data.price);
    ce(:,i) = caltechmargin(s, mm, lambdai, Data.iF) + gammai.*mean(s,2).*Data.gpm./Data.cagpm.*Data.cafe;
    ss(:,i) = mean(s,2);
    fprintf(' Simulation #%d, #iterations = %d, exit flag = %d, distance = %f\n', i, iter, flag, distance);
end

index = all(~isnan(ps),1);
cce = mean(ce(:,index),2)./mean(ss(:,index),2);

end

