function [margin, flag] = calmargin(s, alphai, iF)

F = max(iF);
J = length(iF);
N = size(s,2);

margin = zeros(J,1);
share = mean(s,2);

sa = bsxfun(@times, s, alphai);
meansa = mean(sa,2);
% sa = bsxfun(@rdivide, sa, price);

sharec = cell(1,F);
sac = cell(1,F);
sc = cell(1,F);
meansac = cell(1,F);
marginc = cell(1,F);
parfor f=1:F
    index = iF == f;
    sharec{f} = share(index);
    sc{f} = s(index,:);
    sac{f} = sa(index,:);
    meansac{f} = meansa(index);
end

flag = true(1,F);
parfor f=1:F
    Delta = diag(meansac{f}) - sac{f}*sc{f}'/N;
    singular = rcond(Delta);
    if (singular < eps) || isnan(singular)
        flag(f) = false;
    end
    marginc{f} = -Delta\sharec{f};
end

for f=1:F
    margin(iF == f) = marginc{f};
end

flag = all(flag);

end