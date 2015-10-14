df = importdata('D:\Dropbox\PhD Study\research\data\od\od-monthly-3.csv');
data = df.data;

%%
global lastdelta count outshr

cdid        = data(:,1);
firmid      = data(:,2);


const       = ones(size(cdid));
mpd         = 1./data(:,17:28);
pg          = data(:,29:40);
price       = data(:,41);
hpwt        = data(:,42);
space       = data(:,43);

% od-annual-2.csv only
suv         = data(:,44);
truck       = data(:,45);
van         = data(:,46);
minivan     = data(:,47);

share       = data(:,4:15);
outshr      = data(:,16);

count = 0;

rng('default');

% number of random draws per market
N = 50;

[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');

nT = max(iT);

X  = [const hpwt space suv truck van minivan];

% random coefficients
X1 = [price X];

%% reshape

mark = share > 1e-8;
share = share(mark);

mpd = mpd(mark);
pg  = pg(mark);

X1 = repmat(X1, [12 1]);
X1 = X1(mark(:),:);

iT = bsxfun(@plus, (iT-1)*12, 1:12);
iT = iT(mark);

outshr = repmat(outshr, [1 12]);
outshr = outshr(mark);

%% random draws
K1 = size(X1,2) + 2;
nT = max(iT);

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
sigma0 = rand(1, K1)*1e-2;

lastdelta = delta;
% ssqgrad_monthly(sigma0', share, X1, mpd, pg, v, iT, mark);

%% Opti
optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'derivCheck', 'on', 'solver', 'IPOPT');
Opt = opti('fun', @(sigma) ssq_monthly(sigma, share, X1, mpd, pg, v, iT, mark), ...
    'f', @(sigma) ssqgrad_monthly(sigma, share, X1, mpd, pg, v, iT, mark), ...
    'x0', sigma0, ...
    'lb', -15*ones(size(sigma0)), ...
    'ub', 15*ones(size(sigma0)), ...
    'options', optObj);

[sigma, fval] = solve(Opt);

%%
delta = invertshare(sigma, share, [X1 mpd pg], v, iT);

deltar = NaN(size(mark));
deltar(mark) = delta;

mpdr = NaN(size(mark));
mpdr(mark) = mpd;

pgr = NaN(size(mark));
pgr(mark) = pg;

Ddelta = diff(deltar, 1, 2);
Dmpd = diff(mpdr, 1, 2);
Dpg  = diff(pgr, 1, 2);

mark = ~isnan(Ddelta);
XX = [Dmpd(mark) Dpg(mark)];
YY = Ddelta(mark);
gamma = (XX'*XX)\(XX'*YY);


