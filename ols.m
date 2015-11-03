function [ beta, se, res, Rsqr ] = ols( X, y, lb )
%OLS Summary of this function goes here
%   Detailed explanation goes here

beta = (X'*X)\(X'*y);
res = y - X*beta;
v = var(res)*inv(X'*X);
se = sqrt(diag(v));
Rsqr = 1 - var(res)/var(y);

if nargin > 2
   fprintf(' ===================================================\n');
   printmat([beta, se, beta./se], ' ols result', lb, 'coef se t-value');
   fprintf(' ---------------------------------------------------\n');
   fprintf(' *** obs      = %d\n', size(X,1));
   fprintf(' *** R-square = %f\n', Rsqr);
   fprintf(' ===================================================\n');
end

end

