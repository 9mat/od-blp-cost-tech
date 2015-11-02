function [ ce ] = caltechmargin( s, margin, lambdai, iF )
%CALTECHMARGIN Summary of this function goes here
%   Detailed explanation goes here

F = max(iF);
J = length(iF);
N = size(s,2);

sa = bsxfun(@times, s, lambdai);
sumsa = sum(sa, 2);
ce = zeros(J,1);

for f=1:F
    index = iF == f;
    Delta = (diag(sumsa(index)) - sa(index,:)*s(index,:)')/N;
    ce(index) = Delta*margin(index);
end

end

