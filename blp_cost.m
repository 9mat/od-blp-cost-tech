df = importdata('..\data\od\od-annual-4.csv');
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
Xrc = [const pgreal hpwt space suv truck van minivan];
Krc = size(Xrc,2);

% mean utility coefficients
Xv = [const pgreal hpwt space suv truck van minivan];
Kv = size(Xv,2);

% cost coefficients
Xc = [const trend log(hpwt) log(space) suv truck van minivan];
Kc = size(Xc,2);

% fuel-tech frontier
Xe = [const trend log(hp) log(weight) suv truck van minivan];
Ke = size(Xe,2);

Xzv = [const pgreal/10 hp/10 hpwt space/10 suv truck van minivan];
Xzc = [const log(hpwt) log(space) suv truck van minivan];
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

Xz = blkdiag(Xv, Xc, Xe);
Z = blkdiag(Zv, Zc, Ze);

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
   -0.7662
   -0.3460
   -0.0192
   -0.4322
    2.9170
    0.7499
    0.7612
    0.2066
   -0.6033
   -0.0930
    1.0004
    1.2181
    0.4385
    2.6983
];

lastdelta = delta;

%% combined data
Data.Xrc = Xrc;
Data.iT = iT;
Data.iF = iF;
Data.price = price;
Data.share = share;
Data.v = v(iT,1:end-2,:);
Data.vprice = squeeze(v(iT,end-1,:)); % random draws for price
Data.ve = squeeze(v(iT,end,:)); % random draws for dpm
Data.A = A;
Data.X = Xz;
Data.Z = Z;
Data.ZZ = ZZ;
Data.gpm = gpm;
Data.dpm = dpm;
Data.pgreal = pgreal;

%%
% theta0 = [sigma0 alpha0 lambda0 sigmaprice0 sigmae0 a0 b0];
solveropts = nloptset('algorithm', 'LN_BOBYQA');
optObj = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT', 'solverOpts', solveropts);
optObjIP = optiset('display', 'iter', 'tolrfun', 1e-6, 'tolafun', 1e-6,...
    'maxtime', 1e5, 'solver', 'NLOPT');

thetalb = -10*ones(size(theta0));
thetalb(end-1:end) = 0;
% thetalb(end-3:end-2) = 0;
thetalb(end-5:end-4) = -15;
thetalb(end-4) = -20;
thetaub = 15*ones(size(theta0));
thetaub(end) = 10000;
OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', theta0, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObjIP);

%%
[theta, fval] = solve(OptIP);

%%
Opt = opti('fun', @(theta) gmmcost_rcprice(theta, Data), 'x0', theta, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObj);

[thetaIP, fval] = solve(Opt);
%%
OptIP = opti('fun', @(x) gmmcost_rcprice(x, Data), 'x0', thetaIP, ...
    'lb', thetalb, 'ub', thetaub, 'options', optObjIP);
[thetaIP, fvalIP] = solve(OptIP);

%%

thetaold = theta;
theta = thetaIP;
%%
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

part1 = zeros(size(c));
part2 = zeros(size(c));

lambdai = lambda*exp(sigmae*Data.ve);
margin = Data.price - c;
for f=1:max(iF)
    index = iF == f;
    si = s(index, :);
    part1(index) = Data.pgreal(index).*margin(index,:).*sum(si.*lambdai(index,:),2)/N;
    part2(index) = (lambdai(index,:).*si)*(si'*margin(index)).*Data.pgreal(index)/N;
end

c_e = -(part1-part2)./Data.share;
e = (c_e-a)/(2*b);

eb = log(Data.gpm + e);

c = c - (a*e + b*e.^2);

y = [delta; log(c); eb];
beta = Data.A*y;
gmmres = y - Data.X*beta;
moments = gmmres'*Data.Z;
obj = moments/Data.ZZ*moments';

ce = a*e + b*e.^2;
margin = price - c - ce;
markup = margin./price;


%%
% theta001 = [sigma; alpha]; 
% theta002= [sigmap; sigmae];
% 
% Opt2 = opti('fun', @(theta) gmmcost_rcprice([theta001;theta(1);theta002;theta(2:end)], Data), 'x0', [-12;1.6;-2], ...
%     'lb', [-1000;0; -1000], 'ub', [1000;1000;0], 'options', optObj);
% 
% [theta2, fval2] = solve(Opt2);