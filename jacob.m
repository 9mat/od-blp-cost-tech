function [jab, delta, s] = jacob(theta, Data)
%JACOB Summary of this function goes here
%   Detailed explanation goes here

[delta, s] = invertshare(theta, Data);

params = getParams(theta);

v = Data.v;
X = Data.Xrc;
vp = permute(Data.vprice, [1 3 2]);
ve = permute(Data.ve, [1 3 2]);
iT = Data.iT;

J = size(s,1);
K = numel(theta);
N = size(v,3);

jab = zeros([J, K]);


for t = 1:max(iT)
    index = iT == t;
    si = s(index,:);
    
    %%
    % $$\frac{\partial s_j}{\partial\delta_l} = \frac{1}{N} \sum_i s_{li} 
    % (1_{j=l} - s_{ji}) = s_j - \frac{1}{N} \sum_i s_{li} s_{ji}$$
    derShareDelta = diag(mean(si,2)) - si*si'/N;
    
    %% 
    % indexing order: product - random coef - person (j,k,i)
    si = permute(si, [1 3 2]);
    vv = v(index,:,:);
    xx = X(index,:);
    
    %% 
    % $$\frac{\partial \mu_{ji}}{\partial \sigma_l} = X_{jl} \nu_{li}$$
    %
    % $$\frac{\partial \mu_{ji}}{\partial \alpha} = e^{\sigma^p \nu^p_i}
    % price_j $$
    %
    % $$\frac{\partial \mu_{ji}}{\partial \lambda} = e^{\sigma^e \nu^e_i}
    % dpm_j $$
    %
    % $$\frac{\partial \mu_{ji}}{\partial \sigma^p} = \alpha \nu^p_i
    % e^{\sigma^p \nu^p_i} price_j =  \alpha \nu^p_j \frac{\partial
    % \mu_{ji}}{\partial \alpha}$$
    %
    % $$\frac{\partial \mu_{ji}}{\partial \sigma^e} = \lambda \nu^e_j
    % \frac{\partial \mu_{ji}}{\partial \lambda} $$
    %
    
    derMuSigma = bsxfun(@times, xx, vv);
    derMuAlpha = bsxfun(@times, exp(vp*params.sigmap), Data.price);
    derMuLambda = bsxfun(@times, exp(ve*params.sigmae), Data.dpm);
    derMuSigmap = params.alpha*bsxfun(@times, derMuAlpha, Data.vprice);
    derMuSigmae = params.lambda*bsxfun(@times, derMuLambda, Data.ve);
    
    derMuTheta = cat(2, derMuSigma, derMuAlpha, derMuLambda, derMuSigmap, derMuSigmae);
    
    %%
    % $$ \frac{\partial s_j}{\partial\theta_l} = \frac{1}{N} \sum_i \sum_r
    % \frac{\partial s_{ji}}{\partial \mu_{ri}} \frac{\partial
    % \mu_{ri}}{\partial \theta_l} =  \frac{1}{N} \sum_i \sum_r
    % s_{ji} (1_{j=r} - s_{ri}) \frac{\partial
    % \mu_{ri}}{\partial \theta_l} =  \frac{1}{N} \sum_i \left(
    % s_{ji}\frac{\partial \mu_{ji}}{\partial \theta_l}  - s_{ji} \sum_r
    % s_{ri} \frac{\partial \mu_{ri}}{\partial \theta_l} \right) $$
    
    siderMuTheta = bsxfun(@times, derMuTheta, si);
    sumsiderMuTheta = sum(siderMuTheta, 1);
    ssderMuTheta = bsxfun(@times, si, sumsiderMuTheta);
    derShareTheta = mean(siderMuTheta - ssderMuTheta, 3);
    
    %%
    % $$ \frac{\partial\delta}{\partial\sigma'} = -\frac{\partial
    % s}{\partial \delta}^{-1} \frac{\partial s}{\partial \sigma'}$$
    jab(index,:) = -derShareDelta\derShareTheta;
end

end

