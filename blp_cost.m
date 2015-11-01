clear;
df = importdata('..\data\od\od-annual-6.csv');
data = df.data;

%%
global lastdelta count outshr

cdid        = data(:,1);
firmid      = data(:,2);


const       = ones(size(cdid));
mpd         = data(:,7);
dpm         = data(:,8);
mpg         = data(:,9);
price       = data(:,6);
hp          = data(:,10)/100;
weight      = data(:,11)/1000;
space       = data(:,12);
pgreal      = data(:,13);

suv         = data(:,14);
truck       = data(:,15);
van         = data(:,16);
minivan     = data(:,17);

comply      = data(:,18);
cafestd     = data(:,19)/1.4;

share       = data(:,4);
outshr      = data(:,5);

gpm         = 1./mpg*100; % gallons per 100 miles
dpm         = gpm.*pgreal;
hpwt        = hp./weight;
trend       = cdid-1;

count = 0;

rng('default');

% number of random draws per market
N = 50;

[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');

nT = max(iT);

%% price rc and dpm rc will be dealt with separately

% random coefficients
Xrc = [const pgreal hpwt weight ]; % suv truck van minivan];
Krc = size(Xrc,2);

% mean utility coefficients
Xv = [const pgreal hpwt weight suv truck van minivan];
Kv = size(Xv,2);

% cost coefficients
Xc = [const trend log(hpwt) suv truck van minivan];
Kc = size(Xc,2);

% fuel-tech frontier
Xe = [const trend log(hp) log(weight) suv truck van minivan];
Ke = size(Xe,2);

Xzv = [const hpwt space suv truck van minivan];
Xzc = [const log(hp) log(weight) suv truck van minivan];
Xze = [const log(hp) log(weight) suv truck van minivan];

for k = 1:size(Xzv,2)
    sum_firm_v(:,k) = accumarray(iF, Xzv(:,k));
    sum_total_v(:,k) = accumarray(iT, Xzv(:,k));
end

for k = 1:size(Xzc,2)
    sum_firm_c(:,k) = accumarray(iF, Xzc(:,k));
    sum_total_c(:,k) = accumarray(iT, Xzc(:,k));
end

for k = 1:size(Xze,2)
    sum_firm_e(:,k) = accumarray(iF, Xze(:,k));
    sum_total_e(:,k) = accumarray(iT, Xze(:,k));
end

count_firm = sum_firm_v(iF,1);
count_total = sum_total_v(iT,1);

Z_firm_v = bsxfun(@rdivide, sum_firm_v(iF,2:end) - Xzv(:,2:end), max(count_firm - 1,1));
Z_rival_v = bsxfun(@rdivide, sum_total_v(iT,2:end) - sum_firm_v(iF,2:end), count_total - count_firm);
Zv = [Xzv count_firm/10 Z_firm_v count_total/1000 Z_rival_v];

Z_firm_c = bsxfun(@rdivide, sum_firm_c(iF,2:end) - Xzc(:,2:end), max(count_firm - 1,1));
Z_rival_c = bsxfun(@rdivide, sum_total_c(iT,2:end) - sum_firm_c(iF,2:end), count_total - count_firm);
Zc = [Xzc trend/10 count_firm/10 Z_firm_c count_total/1000 Z_rival_c];

Z_firm_e = bsxfun(@rdivide, sum_firm_e(iF,2:end) - Xze(:,2:end), max(count_firm - 1,1));
Z_rival_e = bsxfun(@rdivide, sum_total_e(iT,2:end) - sum_firm_e(iF,2:end), count_total - count_firm);
Ze = [Xze trend/10 count_firm/10 Z_firm_e count_total/1000 Z_rival_e];

% Xz = blkdiag(Xv, Xc, Xe);
% Z = blkdiag(Zv, Zc, Ze);

Xz = blkdiag(Xv, Xc);
Z = blkdiag(Zv, Zc);

ZZ = Z'*Z; % GMM weighted matrix
XZ = Xz'*Z;
XZZZ = XZ/ZZ;
A = (XZZZ*XZ')\XZZZ*Z';


%% random draws
halton_dim  = (Krc+2)*nT; % last one is for price
halton_skip = 1000;
halton_leap = 100;
halton_scramble    = 'RR2';
tempDraw    = haltonset(halton_dim, 'Skip', halton_skip, 'Leap', halton_leap);
tempDraw    = scramble(tempDraw, halton_scramble);
    
% Make uniform draws
draws       = net(tempDraw, N);
v           = reshape(draws, [N nT (Krc+2)]);
v           = permute(v, [2 3 1]);
v           = norminv(v);

v = bsxfun(@minus, v, mean(v,3));

% IV logit
delta = log(share) - log(outshr);
sigma0 = rand(1, size(Xrc,2))*1e-1;
% alpha0 = (price'*Zv/(Zv'*Zv)*Zv'*price)\(price'*Zv/(Zv'*Zv)*Zv'*delta);
% sigmaprice0 = 0;

alpha0 = -1;
sigmaprice0 = 1; %0.0951;
lambda0 = -0.01;
sigmae0 = 0.5;
b0 = -2;
% sigma0 = [   -0.8175
% %   -1.1868
%    -0.8073
%     0.2405
%     0.0141
%    -0.1248
%     0.2025
%    -7.7324]';

theta0 =    [ 
    %sigma
    0.01
    0.01
    0.01
    0.01

%    
%     1.5757 %sigma cartyve
%     1.0700
%     2.0370
%    -0.4915
   
    % alpha lambda
   -0.1668
   -0.0673
    0.01
    0.01
    
%     % a b
% %     2.0825
%     1
%     
% %     (55/1000)/((1/27-1/28)*100) % shadow cost
%     0.001
];

lastdelta = delta;

%%
theta02 =    [ 
    %sigma
    0.3362
    0.0258
    0.4837
    0.6828
    0

%    
%     1.5757 %sigma cartyve
%     1.0700
%     2.0370
%    -0.4915
   
    % alpha lambda
   -0.1668
   -0.0673
    0.9607
    0.3301    
];

%% combined data
Data.Xv = Xv;
Data.Zv = Zv;
Data.Xrc = Xrc;
Data.iT = iT;
Data.iF = iF;
Data.price = price;
Data.share = share;
Data.v = v(iT,1:end-2,:);
Data.XrcV = bsxfun(@times, Data.Xrc, Data.v);
Data.vprice = squeeze(v(iT,end-1,:)); % random draws for price
Data.ve = squeeze(v(iT,end,:)); % random draws for dpm
Data.A = A;
Data.X = Xz;
Data.Z = Z;
Data.ZZ = ZZ;
Data.gpm = gpm;
Data.dpm = dpm;
Data.pgreal = pgreal;
Data.comply = comply;
Data.cafestd = cafestd;



%%
% theta0 = [sigma0 alpha0 lambda0 sigmaprice0 sigmae0 a0 b0];
solveropts = nloptset('algorithm', 'LN_BOBYQA');
optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT', 'solverOpts', solveropts);

optObjIP = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT');

% thetalb = -10*ones(size(theta0));
% thetalb(end-1) = -10;
% thetalb(end) = -10;
% thetalb(end-6:end-5) = -15;
% thetalb(end-5) = -20;
% thetaub = 15*ones(size(theta0));
% thetaub(end-1) = 100;
% thetaub(end) = 10;

thetalb = -10*ones(size(theta0));
thetaub = 10*ones(size(theta0));

%% Demand only
% OptIP = opti('fun', @(x) gmm(x, Data), 'x0', theta0(1:end-2), ...
%     'lb', thetalb(1:end-2), 'ub', thetaub(1:end-2), 'options', optObjIP);

% [thetaIPdemand, fvalIP] = solve(OptIP);


%%
OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', theta0, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObjIP);

[thetaIP, fvalIP] = solve(OptIP);

OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', thetaIP, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObjIP);

[thetaIP, fvalIP] = solve(OptIP);

% OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% 
% [thetaIP, fvalIP] = solve(OptIP);
% 
% OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% [thetaIP, fvalIP] = solve(OptIP);
% 
% OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% [thetaIP, fvalIP] = solve(OptIP);


% OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObj);
% [thetaIP, fvalIP] = solve(OptIP);

%%
% [theta, fval] = solve(OptIP);

%%
% Opt = opti('fun', @(theta) gmmcost_rcprice(theta, Data), 'x0', theta, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObj);
% 
% [thetaIP, fval] = solve(Opt);
% %%
% OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);


%%

% thetaold = theta;

theta = thetaIP;
%%
params = getParams(theta);
alpha = params.alpha;
lambda = params.lambda;
sigmap = params.sigmap;
sigmae = params.sigmae;
% b = params.b;
% gamma = params.gamma;

% invert for delta
[delta, s] = invertshare(theta, Data);

% % shadow cost
% gammaj = zeros(size(Data.price));
% gammaj(Data.comply <= 0) = gamma;
% gammaj(Data.comply < 0) = (55/1000)/((1/27-1/28)*100);

alphai = alpha*exp(sigmap*Data.vprice);

% shadow_cost = gammaj.*(Data.gpm - 1./(Data.cafestd/100));

margin = calmargin(s, alphai, Data.iF);
c = Data.price - margin;
markup = margin./Data.price;

% part1 = zeros(size(c));
% part2 = zeros(size(c));
% 
% lambdai = lambda*exp(sigmae*Data.ve);
% margin = Data.price - c - shadow_cost;
% for f=1:max(iF)
%     index = iF == f;
%     si = s(index, :);
%     part1(index) = Data.pgreal(index).*margin(index,:).*sum(si.*lambdai(index,:),2)/N;
%     part2(index) = (lambdai(index,:).*si)*(si'*margin(index)).*Data.pgreal(index)/N;
% end
% 
% c_e = -(part1-part2)./Data.share - gammaj;
% e = c_e/(2*b);
% 
% eb = log(Data.gpm + e);
% 
% c = c - (b*e.^2);

y = [delta; log(c)];
beta = Data.A*y;
% gmmres = y - Data.X*beta;
% moments = gmmres'*Data.Z;
% obj = moments/Data.ZZ*moments';

% ce = b*e.^2;
% margin = price - c - ce - shadow_cost;
% markup = margin./price;


%%

V = cov(theta, beta, Data);
printmat([theta, se(1:numel(theta))], 'theta', '', 'theta se');
