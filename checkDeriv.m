function [ err ] = checkDeriv(x, f, g )
%CHECKDDERIV Summary of this function goes here
%   Detailed explanation goes here

fx = f(x);
df = zeros([numel(fx) numel(x)]);
epsilon = 1e-7;

for i = 1:length(x)
    x1 = x;
    x2 = x;
    dx = epsilon;
    x2(i) = x(i) + dx;
    x1(i) = x(i) - dx;
    df(:,i) = (f(x2) - f(x1))/(2*dx);
end

gx = g(x);
err = gx - df;

end

