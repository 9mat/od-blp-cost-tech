function [margin, flag] = calmargin(s, alphai, iF)

F = max(iF);
J = length(iF);
N = size(s,2);

margin = zeros(J,1);
share = mean(s,2);

sa = bsxfun(@times, s, alphai);
meansa = mean(sa,2);
% sa = bsxfun(@rdivide, sa, price);

flag = true;
for f=1:F
    index = iF == f;
    Delta = diag(meansa(index)) - sa(index,:)*s(index,:)'/N;
%     singular = rcond(Delta);
%     if (singular < eps) || isnan(singular)
%         flag(f) = false;
%     end
    margin(index) = -Delta\share(index);
end


end