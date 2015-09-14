function [ d, dd ] = ddpi_price( p, e, delta, theta, Data  )
%DDPI_PRICE Summary of this function goes here
%   Detailed explanation goes here
J = length(p);

if nargout == 1
    a = ddpi(p, e, delta, theta, Data );
    d = a(1:J);
else
    [a, aa] = ddpi(p, e, delta, theta, Data);
    d = a(1:J);
    dd = aa(1:J, 1:J);
end

end

