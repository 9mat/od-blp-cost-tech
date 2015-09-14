function [ ssq ] = ssq_dpi( p, e, delta, theta, Data )
%SSQ_DPI Summary of this function goes here
%   Detailed explanation goes here

ssq = sum(ddpi(p, e, delta, theta, Data).^2);

end

