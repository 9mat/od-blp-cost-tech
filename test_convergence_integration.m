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

% IV logit
delta = log(share) - log(outshr);
sigma = [    0.1132
    0.1418
         0
    1.9103
    0.4034
    1.3266
    3.2403
    2.9384
    4.2682];

lastdelta = delta;


%% pseudo-random draws
rng('default');

for i = 1:10
    N = i*100;
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
%     v = normrnd(0,1,[nT K1 N]);
    gmmvalue(i) = gmm(sigma, share, X1, v, iT, X2, Z);
    fprintf('N=%d\tGMM=%.5f\n',N, gmmvalue(i));
end

% see result in convergence-simulated-gmm.fig
