function [ cce, ps, gammajs, share ] = calcce(theta, deltas, cs, Data, gammaj, ps0, settings )
%CALCCE Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta);
lambdai = -bsxfun(@rdivide, params.lambda*bsxfun(@times, exp(params.sigmae*Data.ve), Data.dpm), Data.income09);

% binding = Data.comply == 0;
% fined = Data.comply == -1;
% gammai = zeros(size(Data.comply));
% gammai(binding) = params.gamma;
% gammai(fined) = params.gamma + 0.055;S

ns = min(size(deltas,2), size(cs,2));
J = size(deltas,1);
ce = zeros(J,ns);
ss = zeros(J,ns);
ps = zeros(J,ns);
gammajs = zeros(J,ns);

if nargin < 6
    ps0 = repmat(Data.price, [1 ns]);
end

for i=1:ns
%     [ps(:,i), mm, s, iter, flag, distance] = contraction_bertrand(theta, deltas(:,i), cs(:,i), Data, gammaj, ps0(:,i));
    [ps(:,i), mm, s, gammajs(:,i), cafe, iter, distance] ...
        = contraction_cafe(theta, deltas(:,i), cs(:,i), Data, gammaj, ps0(:,i), settings);
    cagpm = 1./cafe*100;
    ce(:,i) = caltechmargin(s, mm, lambdai, Data.iF) ...
        + gammajs(:,i).*mean(s,2).*Data.gpm./cagpm.*cafe;
    ss(:,i) = mean(s,2);
%     fprintf(' Simulation #% 5d, #iterations = % 3d, distance = %f\n', i, iter, distance);
end

index = all(~isnan(ps),1);
share = mean(ss(:,index),2);
cce = mean(ce(:,index),2)./mean(ss(:,index),2)/10;

% trim outliers
trim05 = prctile(cce, 3);
trim95 = prctile(cce, 97);
cce(cce>trim95) = trim95;
cce(cce<trim05) = trim05;

end

