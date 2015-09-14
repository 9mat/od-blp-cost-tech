load BLP_data.mat

global lastdelta

rng('default');

% number of random draws per market
N = 10;


[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');


mpd = 1./mpd;

X  = [const  mpd];

% random coefficients
X1 = [const mpd];

% mean coefficients
X2 = [price X];

for k = 1:size(X,2)
    sum_firm(:,k) = accumarray(iF, X(:,k));
    sum_total(:,k) = accumarray(iT, X(:,k));
end

Z_firm = sum_firm(iF,:) - X;
Z_rival = sum_total(iT, :) - sum_firm(iF, :);
Z = [X Z_firm Z_rival];

% random draws
v = normrnd(0,1, [max(iT), size(X1, 2), N]);
v = bsxfun(@minus, v, mean(v,3));

% IV logit
delta = log(share) - log(outshr);
sigma0 = zeros(1, size(X1,2));

lastdelta = delta;

optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6, 'maxtime', 1e4, 'solver', 'NLOPT');
Opt = opti('fun', @(sigma) gmm(sigma, share, X1, v, iT, X2, Z), 'x0', sigma0, 'options', optObj);

%%
sigma = solve(Opt);

%%
[delta, s] = invertshare(sigma, share, X1, v, iT);
XZ = X2'*Z;
ZZ = Z'*Z;
XZZZ = XZ/ZZ;
Zdelta = Z'*delta;
beta = (XZZZ*XZ')\(XZZZ*Zdelta);

% Results
% sigma =
% 
%     0.4053
%     3.5831
%     1.4683
%     0.0034
%     0.3646
%     0.0000
    
margin = calmargin(s, beta, sigma, price, v, F, iF, iT);