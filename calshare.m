function s = calshare(delta, theta, Data)

params  = getParams(theta);
sigma   = params.sigma;
alpha   = params.alpha;
lambda  = params.lambda;
sigmap  = params.sigmap;
sigmae  = params.sigmae;

% mu(jti) = delta(jt) + sum_k sigma(k)*X(jk)*v(jki) +
% alpha*exp(sigmap*vp(i))*price(jt)
mu = bsxfun(@times, sigma', Data.Xrc);
mu = bsxfun(@times, mu, Data.v);
mup = alpha*bsxfun(@times,exp(sigmap*Data.vprice),Data.price);
mue = lambda*bsxfun(@times,exp(sigmae*Data.ve),Data.dpm);
mu = squeeze(sum(mu, 2)) + mup + mue;
mu = bsxfun(@plus, delta, mu);

N = size(Data.v,3);
iT = Data.iT;
mxmu = zeros(max(iT), N);
for i=1:N
    mxmu(:,i) = accumarray(iT, mu(:,i), [], @max);
end

emu2 = exp(mu - mxmu(iT,:));
s0 = zeros(max(iT), N);
for i=1:N
    s0(:,i) = exp(-mxmu(:,i)) + accumarray(iT, emu2(:,i));
end

s = emu2./s0(iT,:);

end