function [ grad ] = grad_ssq_dpi( p, e, delta, theta, Data )
%GRAD_SSQ_DPI Summary of this function goes here
%   Detailed explanation goes here

[d,dd] = ddpi(p,e,delta, theta, Data);
grad = (2*dd'*d)';

end

