function margin = calmargin(s, alphai, iF)

F = max(iF);
J = length(iF);
N = size(s,2);

margin = zeros(J,1);
share = mean(s,2);

sa = bsxfun(@times, s, alphai);
meansa = mean(sa,2);
% sa = bsxfun(@rdivide, sa, price);

for f=1:F
    index = iF == f;
    Delta = diag(meansa(index,:)) - sa(index,:)*s(index,:)'/N;
    margin(index) = -Delta\share(index);
end

end