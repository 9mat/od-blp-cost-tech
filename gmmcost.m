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
binding = Data.comply == 0;
fined = Data.comply == -1;

gammai = zeros(size(Data.comply));
gammai(binding) = params.gamma;
gammai(fined) = params.gamma + 0.05;

comply_mc1 = (1-Data.gpm./Data.cagpm).*Data.cafe;
comply_mc2 = (1-Data.cagpm./Data.cagpmstd).*Data.cafestd;
comply_mc2(binding) = 0;
comply_mc = gammai.*(comply_mc1 + comply_mc2);

comply_mc(isnan(comply_mc)) = 0;

alphai = params.alpha*exp(params.sigmap*Data.vprice);
margin = calmargin(s, alphai, Data.iF);
c = Data.price - margin - comply_mc;
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

