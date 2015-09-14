function margin = calmargin(s, beta, sigma, price, v, F, iF, iT)

margin = zeros(size(iF));
share = mean(s,2);
N = size(s, 2);

alpha = beta(1);
for f=1:length(F)
    index = iF == f;
%     alphaa = mean(alpha(index,:), 1)';
    si = s(index, :);
    ss = share(index);
    Delta  = (diag(sum(si,2)) - si*si')*alpha/N;
    margin(index) = -Delta\ss;
end

end