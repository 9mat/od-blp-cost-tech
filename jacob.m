function [jab, delta, s] = jacob(sigma, share, X, v, iT)
%JACOB Summary of this function goes here
%   Detailed explanation goes here

[delta, s] = invertshare(sigma, share, X,v,iT);

J = size(s,1);
K = size(X,2);
N = size(v,3);

v = v(iT,:,:);
jab = zeros([J, K]);

for t = 1:max(iT)
    index = iT == t;
    ss = s(index,:);
    
    %%
    % $$\frac{\partial s_j}{\partial\delta_l} = \frac{1}{N} \sum_i s_{li} 
    % (1_{j=l} - s_{ji}) = s_j - \frac{1}{N} \sum_i s_{li} s_{ji}$$
    derShareDelta = diag(mean(ss,2)) - ss*ss'/N;
    
    % indexing order: product - random coef - person (j,k,i)
    ss = permute(ss, [1 3 2]);
    vv = v(index,:,:);
    xx = X(index,:);
    
    xv = bsxfun(@times, xx, vv); % x_{jk}*v_{ki}
    sxv = bsxfun(@times, xv, ss);  % s_{ji}*x_{jk}*v_{ki}
    sumsxv = sum(sxv, 1); % sum_j s_{ji}*x_{jk}*v_{ki}
    ssxv = bsxfun(@times, ss, sumsxv); % s_{ji}*sum_r(s_{ri}*x_{rk}*v_{ki})
    
    %%
    % $$\frac{\partial s_j}{\partial\sigma_k} = \frac{1}{N} \sum_i s_{ji}
    % (v_{ki} X_{jk} - \sum_r s_{ri} v_{ki} X_{rk})$$
    derShareSigma = mean(sxv - ssxv, 3);
    
    %%
    % $$ \frac{\partial\delta}{\partial\sigma'} = -\frac{\partial
    % s}{\partial \delta}^{-1} \frac{\partial s}{\partial \sigma'}$$
    jab(index,:) = -derShareDelta\derShareSigma;
end

end

