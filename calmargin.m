function margin = calmargin(s, alphai, iF)

F = max(iF);
J = length(iF);
N = size(s,2);

margin = zeros(J,1);
share = mean(s,2);

sa = bsxfun(@times, s, alphai);
% sa = bsxfun(@rdivide, sa, price);

for f=1:F
    index = iF == f;
    Delta = diag(mean(sa(index,:),2)) - sa(index,:)*s(index,:)'/N;
    margin(index) = -Delta\share(index);
end

end