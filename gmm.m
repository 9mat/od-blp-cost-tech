function obj = gmm(sigma, share, X, v, iT, X2, Z)
%SSQ Summary of this function goes here
%   Detailed explanation goes here

% invert for delta
%sigma = exp(sigma);
delta = invertshare(sigma, share, X, v, iT);

if any(delta > 1e30)
    obj = 1e30;
    return
end

XZ = X2'*Z;
ZZ = Z'*Z;
XZZZ = XZ/ZZ;
Zdelta = Z'*delta;
beta = (XZZZ*XZ')\(XZZZ*Zdelta);
% count = count + 1
% beta'

xi = delta - X2*beta;

xiZ = xi'*Z;
obj = xiZ/ZZ*xiZ';


% IV
% avgdpm = accumarray(iJ, dpm)./count;
% dpm1 = dpm - avgdpm(iJ);
% 
% Z = [dpm1 pg1];
% xiZ = xi1'*Z;
% obj = xiZ/(Z'*Z)*(xiZ');
% 

end

