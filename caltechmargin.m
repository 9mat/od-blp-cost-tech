function [ ce ] = caltechmargin( s, margin, lambdai, iF )
%CALTECHMARGIN Summary of this function goes here
%   Detailed explanation goes here

F = max(iF);
J = length(iF);
N = size(s,2);

sa = bsxfun(@times, s, lambdai);
sumsa = sum(sa, 2);
ce = zeros(J,1);

sc = cell(1,F);
sac = cell(1,F);
sumsac = cell(1,F);
cec = cell(1,F);
marginc = cell(1,F);
for f=1:F
    index = iF == f;
    sc{f} = s(index,:);
    sac{f} = sa(index,:);
    marginc{f} = margin(index);
    sumsac{f} = sumsa(index);
end

parfor f=1:F
%     index = iF == f;
    Delta = (diag(sumsac{f}) - sac{f}*sc{f}')/N;
    cec{f} = Delta*marginc{f};
end

for f=1:F
    ce(iF == f) = cec{f};
end

end

