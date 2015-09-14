function grad = gmmgrad(sigma, share, X, v, iT, X2, Z )
%GMMGRAD Summary of this function goes here
%   Detailed explanation goes here


% invert for delta
%sigma = exp(sigma);
[jab, delta, ~] = jacob(sigma, share, X, v, iT);

if any(delta > 1e30)
    grad = 1e30*ones(size(sigma));
    return
end

XZ = X2'*Z;
ZZ = Z'*Z;
XZZZ = XZ/ZZ;
Zdelta = Z'*delta;
Zjab = Z'*jab;
xi = delta - X2*((XZZZ*XZ')\(XZZZ*Zdelta)); 
dxi = jab - X2*((XZZZ*XZ')\(XZZZ*Zjab));
% count = count + 1
% beta'

xiZ = xi'*Z;
grad = 2*dxi'*Z/ZZ*xiZ';
grad = grad';

% IV
% avgdpm = accumarray(iJ, dpm)./count;
% dpm1 = dpm - avgdpm(iJ);
% 
% Z = [dpm1 pg1];
% xiZ = xi1'*Z;
% obj = xiZ/(Z'*Z)*(xiZ');
% 


end

