function obj = gmmcost_rcprice(theta, Data)
%SSQ Summary of this function goes here
%   Detailed explanation goes here

alpha = theta(end-5);
lambda = theta(end-4);
sigmap = theta(end-3);
sigmae = theta(end-2);
a = theta(end-1);
b = theta(end);

% invert for delta
[delta, s] = invertshare(theta, Data);

if any(delta > 1e30 | isnan(delta))
    fprintf('Overflow in inverting shares!\n');
    obj = 1e30 + sum(theta)*1e4;
    return
end

iF = Data.iF;
c = zeros(size(iF));
N = size(s, 2);
alphai = alpha*exp(sigmap*Data.vprice);

for f=1:max(iF)
    index = iF == f;
    si = s(index, :);
    ss = Data.share(index);
    sv = si.*alphai(index,:);
    Delta  = (diag(sum(sv,2)) - sv*si')/N;
    c(index) = Data.price(index) + Delta\ss;
end

if any(c <= 0)
    fprintf('Negative cost!\n');
    obj = 1e30 + sum(theta)*1e4;
    return;
end

if any(c > Data.price)
    fprintf('Negative margin!\n');
    obj = 1e30 + sum(theta)*1e4;
    return;
end

e = zeros(size(c));
lambdai = lambda*exp(sigmae*Data.ve);
margin = Data.price - c;
for f=1:max(iF)
    index = iF == f;
    si = s(index, :);
    sv = si.*lambdai(index,:);
    Delta = (diag(sum(sv,2)) - sv*si')/N;
    Delta = bsxfun(@times, Delta, (margin(index).*Data.pgreal(index))');
    e(index) = (Data.gpm(index)./Data.share(index).*sum(Delta,2)/(a*b)).^(1/b);
end

if any(e <= 0)
    fprintf('Negative techonology!\n');
    obj = 1e30 + sum(theta)*1e4;
    return;
end

if any(e > 1)
    fprintf('Damaged goods!\n');
    obj = 1e30 + sum(theta)*1e4;
    return;
end


eb = log(Data.gpm) - log(e);

c = c - a*e.^b + a;
if any(c <= 0)
    fprintf('Negative cost after technology!\n');
    obj = 1e30 + sum(theta)*1e4;
    return;
end

y = [delta; log(c); eb];
beta = Data.A*y;
gmmres = y - Data.X*beta;
moments = gmmres'*Data.Z;
obj = moments/Data.ZZ*moments';

end

