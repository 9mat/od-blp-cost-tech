function obj = gmmcost(theta, Data)
%SSQ Summary of this function goes here
%   Detailed explanation goes here

params = getParams(theta);

[delta, s] = invertshare(theta, Data);

if any(delta > 1e30 | isnan(delta))
    fprintf('Overflow in inverting shares!\n');
    obj = 1e30;
    return
end

% marginal compliance cost
comply_mc = calcomply_mc(params.gamma, Data);
alphai = bsxfun(@rdivide, params.alpha*exp(params.sigmap*Data.vprice), Data.income09);
margin = calmargin(s, alphai, Data.iF);
c = Data.price - margin + comply_mc;
logc = log(c);
logc(c<=0) = -1e30;

% if any(c <= 0)
%     fprintf('Negative cost!\n');
%     obj = 1e30;
%     return;
% end
% 
% if any(margin < 0)
%     fprintf('Negative margin!\n');
%     obj = 1e30;
%     return;
% end
% 
y = [delta; logc];
theta = Data.A*y;
gmmres = y - Data.X*theta;
moments = gmmres'*Data.Z;
obj = moments/Data.ZZ*moments';

end

