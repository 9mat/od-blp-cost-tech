function mu = calmu(theta, Data )
%CALMU Summary of this function goes here
%   Detailed explanation goes here

%%
%
% $$\mu_{jti} = \delta_{jt} + \sum_k \sigma_k X_{jk} \nu_{jki} +
% \alpha \exp(\sigma^p \nu^p_{ki} price_{jt}) +
% \lambda \exp(\sigma^e \nu^e_{ki} dpm_{jt})$$
%

% N = size(Data.v, 3);
% J = size(Data.Xrc, 1);

params = getParams(theta);
alphai = params.alpha*exp(params.sigmap*Data.vprice);
lambdai = params.lambda*exp(params.sigmae*Data.ve);

mu = squeeze(sum(bsxfun(@times, params.sigma', Data.XrcV),2)) ...
    + bsxfun(@times, alphai, Data.price) ...
    + bsxfun(@times, lambdai, Data.dpm);



end

