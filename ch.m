df = importdata('D:\Dropbox\PhD Study\research\data\od\od-annual-2.csv');
data = df.data;

%%
global lastdelta count outshr

cdid        = data(:,1);
firmid      = data(:,2);


const       = ones(size(cdid));
mpd         = data(:,7);
price       = data(:,6);
hpwt        = data(:,8);
space       = data(:,9);

% od-annual-2.csv only
suv         = data(:,10);
truck       = data(:,11);
van         = data(:,12);
minivan     = data(:,13);

share       = data(:,4);
outshr      = data(:,5);

count = 0;

rng('default');

% number of random draws per market
N = 100;

[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');

nT = max(iT);

X  = [const mpd hpwt space suv truck van minivan];

% random coefficients
X1 = [price X];
K1 = size(X1,2);

% mean coefficients
X2 = [price X];
K2 = size(X2,2);

for k = 1:size(X,2)
    sum_firm(:,k) = accumarray(iF, X(:,k));
    sum_total(:,k) = accumarray(iT, X(:,k));
end
%%
count_firm = sum_firm(iF,1);
count_total = sum_total(iT,1);
Z_firm = bsxfun(@rdivide, sum_firm(iF,2:end) - X(:,2:end), max(count_firm - 1,1));
Z_rival = bsxfun(@rdivide, sum_total(iT,2:end) - sum_firm(iF,2:end), count_total - count_firm);
Z = [X count_firm Z_firm count_total Z_rival];

%% random draws
halton_dim  = K1*nT;
halton_skip = 1000;
halton_leap = 100;
halton_scramble    = 'RR2';
tempDraw    = haltonset(halton_dim, 'Skip', halton_skip, 'Leap', halton_leap);
tempDraw    = scramble(tempDraw, halton_scramble);
    
% Make uniform draws
draws       = net(tempDraw, N);
v           = reshape(draws, [N nT K1]);
v           = permute(v, [2 3 1]);
v           = norminv(v);

% IV logit
delta = log(share) - log(outshr);
sigma0 = rand(1, size(X1,2))*1e-2;

lastdelta = delta;

%%
err = checkDeriv(sigma0', @(x) gmm(x,share,X1,v,iT,X2,Z), @(x) gmmgrad(x,share,X1,v,iT,X2,Z));
%%
optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'derivCheck', 'on', 'solver', 'IPOPT');
Opt = opti('fun', @(sigma) gmm(sigma, share, X1, v, iT, X2, Z), 'f', @(sigma) gmmgrad(sigma, share, X1, v, iT, X2, Z), 'x0', sigma0, ...
    'lb', zeros(size(sigma0)), 'ub', 15*ones(size(sigma0)), 'options', optObj);

%%
[sigma, fval] = solve(Opt);

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