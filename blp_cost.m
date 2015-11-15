diary diary.txt;
diary on;
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

fleet       = 1 - suv - truck - van - minivan; % car or truck;

comply      = data(:,18);
cafestd     = data(:,19)*0.7;
cafe        = data(:,20)*0.7;
torque      = data(:,21)/100;
income09    = data(:,24)/1000;
mampg       = data(:,23);
gdppc       = data(:,22)/1000;


share       = data(:,4);
outshr      = data(:,5);

gpm         = 1./mpg*100; % gallons per 100 miles
dpm         = gpm.*pgreal;
madpm       = 1./mampg*100.*pgreal;
hpwt        = hp./weight;
trend       = cdid-1;

cagpm       = 1./cafe*100;
cagpmstd    = 1./cafestd*100;

count = 0;

rng('default');

% number of random draws per market
N = 50;

[T, ~, iT] = unique(cdid);
[F, ~, iF] = unique([cdid, firmid], 'rows');

[~, ~, fleet] = unique([iF, fleet], 'rows');
cafe2 = accumarray(fleet, share)./accumarray(fleet, share.*(gpm/100));
cafe2 = cafe2(fleet);
nT = max(iT);

%% price rc and dpm rc will be dealt with separately

% random coefficients
Xrc_lb = 'const hpwt weight space';
Xrc = [const hpwt weight space]; % suv truck van minivan];
Krc = size(Xrc,2);

% mean utility coefficients
Xv_lb = 'const madpm hpwt weight space suv truck van minivan';
Xv = [const madpm./income09 hpwt weight space suv truck van minivan];
Kv = size(Xv,2);

% cost coefficients
Xc_lb = 'const trend log(hpwt) log(weight) log(space) log(gpm) suv truck van minivan';
Xc = [const trend log(hpwt) log(weight) log(space) log(gpm) suv truck van minivan];
Kc = size(Xc,2);

% fuel-tech frontier
Xe = [const trend log(hp) log(weight) log(space) suv truck van minivan];
Ke = size(Xe,2);

Xzv = [const dpm hpwt space torque weight suv truck van minivan];
Xzc = [const log(hp) log(weight) log(space) suv truck van minivan];
Xze = [const log(hp) log(weight) log(space) suv truck van minivan];

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

st = 2;
Z_firm_v = bsxfun(@rdivide, sum_firm_v(iF,st:end) - Xzv(:,st:end), max(count_firm - 1,1));
Z_rival_v = bsxfun(@rdivide, sum_total_v(iT,st:end) - sum_firm_v(iF,st:end), count_total - count_firm);
Zv = [Xzv gdppc pgreal cafestd.*comply count_firm/10 Z_firm_v count_total/1000 Z_rival_v];

st = 2;
Z_firm_c = bsxfun(@rdivide, sum_firm_c(iF,st:end) - Xzc(:,st:end), max(count_firm - 1,1));
Z_rival_c = bsxfun(@rdivide, sum_total_c(iT,st:end) - sum_firm_c(iF,st:end), count_total - count_firm);
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
    1
%     0.1
    -1
    -0.2
    0.1
%    
%     1.5757 %sigma cartyve
%     1.0700
%     2.0370
%    -0.4915
   
    % alpha lambda
   -0.4*66
   -0.3*66
    1
    1
    
%     % a b
% %     2.0825
%     1
%     
% %     (55/1000)/((1/27-1/28)*100) % shadow cost
    0.1
];

lastdelta = delta;

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
Data.cafe = cafe;
Data.cagpm = cagpm;
Data.cagpmstd = cagpmstd;
Data.income09 = income09;
Data.fleet = fleet;
Data.cafe2 = cafe2;


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

thetalb = -100*ones(size(theta0));
thetaub = 50*ones(size(theta0));

%% Demand only
% OptIP = opti('fun', @(x) gmm(x, Data), 'x0', theta0(1:end-2), ...
%     'lb', thetalb(1:end-2), 'ub', thetaub(1:end-2), 'options', optObjIP);

% [thetaIPdemand, fvalIP] = solve(OptIP);


%%
OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', theta0, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObjIP);

[thetaIP, fvalIP] = solve(OptIP);

% thetaIP = thetaIP.*(1 + (rand(size(thetaIP))-0.5)*0.2);
% OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% 
% [thetaIP2, fvalIP2] = solve(OptIP);
% if fvalIP2 < fvalIP; thetaIP = thetaIP2; end
% 
% thetaIP = thetaIP.*(1 + (rand(size(thetaIP))-0.5)*0.2);
% OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% 
% [thetaIP2, fvalIP2] = solve(OptIP);
% if fvalIP2 < fvalIP; thetaIP = thetaIP2; end
% 
% thetaIP = thetaIP.*(1 + (rand(size(thetaIP))-0.5)*0.2);
% OptIP = opti('fun', @(x) gmmcost(x, Data), 'x0', thetaIP, ...
%     'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
% 
% [thetaIP2, fvalIP2] = solve(OptIP);
% if fvalIP2 < fvalIP; thetaIP = thetaIP2; end

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

% shadow_cost = gammaj.*(Data.gpm - 1./(Data.cafestd/100));
[comply_mc, gammai] = calcomply_mc(params.gamma, Data);

alphai = bsxfun(@rdivide, params.alpha*exp(params.sigmap*Data.vprice), Data.income09);
margin = calmargin(s, alphai, Data.iF);
c = Data.price - margin + comply_mc;
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

% printmat(theta, 'theta', [Xrc_lb ' alpha lambda sigmap sigmae gamma'], 'theta');
% printmat(beta(1:Kv), 'beta_v', Xv_lb, 'beta_v');
% printmat(beta(1:Kc), 'beta_c', Xc_lb, 'beta_c');

%%

V = cov(theta, beta, Data);
se = sqrt(diag(V));

beta_v = beta(1:Kv);
beta_c = beta(Kv+1:end);

se_theta = se(1:numel(theta));
se_beta = se(numel(theta)+1:end);
se_beta_v = se_beta(1:Kv);
se_beta_c = se_beta(Kv+1:end);

printmat([theta se_theta theta./se_theta], 'theta', [Xrc_lb ' alpha lambda sigmap sigmae gamma'], 'theta se t');
printmat([beta_v se_beta_v beta_v./se_beta_v], 'beta_v', Xv_lb, 'beta_v se t');
printmat([beta_c se_beta_c beta_c./se_beta_c], 'beta_c', Xc_lb, 'beta_c se t');

diary off;
%%
J = length(price);
res = y - Data.X*beta;
xi = res(1:J);
omega = res(J+1:end);

%% simulate
ns = 5;
xis = xi(randi(J, [J,ns]));
omegas = omega(randi(J, [J,ns]));

deltas = bsxfun(@plus, Data.Xv*beta_v, xis);
cs = exp(bsxfun(@plus, Xc*beta_c, omegas));

[cce, ps] = calcce(theta, deltas, cs, Data);

%%
% [~, iii] = sort(cce);
% pct = zeros(size(cce));
% pct(iii) = (1:length(cce))'/length(cce);
% cce = pct;
trim05 = prctile(cce, 3);
trim95 = prctile(cce, 97);
index = (cce > trim05) & (cce < trim95);

torque = data(:,21);
cdiddummies = dummyvar(cdid);
cdiddummies(:,1) = [];

xplot = sort(cce(index));
colors = 'ymcrgbkp';
for pow = 1:7
    poly = bsxfun(@power, cce, 0:pow);
    Xe = [log(hpwt) log(weight) log(space) log(torque) suv minivan van truck poly];
    ye = -log(gpm);
    eta = ols(Xe(index,:),ye(index));
    coef = eta(end-pow:end);
    f = @(x) bsxfun(@power, x, 0:pow)*coef;
    hold on; plot(xplot, f(xplot), colors(pow));
end
legend('1','2','3','4','5','6','7');

poly = bsxfun(@power, cce, 0:3);
Xe_lb = ['log(hpwt) log(weight) log(space) log(torque) suv minivan van truck 2 3 4 5 6 7 8 9 const cce cce^2 cce^3 cce^4 cce^5'];
Xe =  [log(hpwt) log(weight) log(space) log(torque) suv minivan van truck cdiddummies poly];
ye = -log(gpm);
[eta, se] = ols(Xe(index,:), ye(index), Xe_lb);

%%
% Data.pgreal = pgreal*0.99;
% Data.gpm = gpm*0.5;
% Data.dpm = dpm*0.5*0.99;
% [cce, ps] = calcce(theta, deltas, cs, Data, ps);
%%

Data.pgreal = pgreal*0.99;
coef = -eta(end-3:end);
gpm1 = contraction_tech(theta, deltas(:,1:5), cs(:,1:5), Data, cce, coef, ps);

%% hinge function

hinge1 = max(cce-10, 0);
hinge2 = max(cce-30, 0);
hinge3 = max(cce-60, 0);
hinge4 = max(cce-80, 0);
hinge5 = max(cce-6, 0);
hinge6 = max(cce-8, 0);
hinges = [hinge1 hinge2 hinge3];
torque = data(:,21);

Xe_lb = ['log(hpwt) log(weight) log(space) log(torque) suv minivan van truck ' num2str(2:max(cdid)) ' const hinge1 hinge2 hinge3 hinge4 hinge5'];
Xe = [log(hpwt) log(weight) log(space) log(torque) suv minivan van truck cdiddummies const hinges];
ye = -log(gpm);
[eta2, se2] = ols(Xe(index,:), ye(index), Xe_lb);