function [ params ] = getParams( theta )
%GETPARAMS Summary of this function goes here
%   Detailed explanation goes here

% params.sigma = theta(1:end-6);
% params.alpha = theta(end-5);
% params.lambda = theta(end-4);
% params.sigmap = theta(end-3);
% params.sigmae = theta(end-2);
% params.b = theta(end-1);
% params.gamma = theta(end);


params.sigma    = theta(1:end-7);
params.alpha    = theta(end-6);
params.lambda   = theta(end-5);
params.lambda2  = theta(end-4);
params.sigmap   = theta(end-3);
params.sigmae   = theta(end-2);
params.sigmae2  = theta(end-1);
params.gamma    = theta(end);

% params.sigma = theta;
end

