function obj = gmmcost(param, Data)
%SSQ Summary of this function goes here
%   Detailed explanation goes here

sigma = param(1:end-1);
alpha = param(end);

% invert for delta
[delta, s] = invertshare(sigma, Data.share, Data.Xrc, Data.v, Data.iT);

if any(delta > 1e30 | isnan(delta))
    fprintf('Overflow in inverting shares!\n');
    obj = 1e30;
    return
end

iF = Data.iF;
c = zeros(size(iF));
N = size(s, 2);
for f=1:max(iF)
    index = iF == f;
    si = s(index, :);
    ss = Data.share(index);
    Delta  = (diag(sum(si,2)) - si*si')*alpha/N;
    c(index) = Data.price(index) + Delta\ss;
end

if any(c <= 0)
    fprintf('Negative cost!\n');
    obj = 1e30;
    return;
end

if any(c > Data.price)
    fprintf('Negative margin!\n');
    obj = 1e30;
    return;
end

y = [delta - alpha*Data.price; log(c)];
theta = Data.A*y;
gmmres = y - Data.X*theta;
moments = gmmres'*Data.Z;
obj = moments/Data.ZZ*moments';

end

