%% pred.m
% *Summary:* Compute predictive (marginal) distributions of a trajecory
%
%   [M S] = pred(policy, plant, dynmodel, m, s, H)
%
% *Input arguments:*
%
%   policy             policy structure
%   plant              plant structure
%   dynmodel           dynamics model structure
%   m                  D-by-1 mean of the initial state distribution
%   s                  D-by-D covariance of the initial state distribution
%   H                  length of prediction horizon
%
% *Output arguments:*
%
%   M                  D-by-(H+1) sequence of predicted mean vectors
%   S                  D-by-D-(H+1) sequence of predicted covariance
%                      matrices
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-01-23
%
%% High-Level Steps
% # Predict successor state distribution

% function [M, S, Mcon, Scon] = my_pred(policy, plant, dynmodel, m, s, H)
% %% Code
% D = length(m); nU = length(policy.maxU);
% S = zeros(D,D,H+1); M = zeros(D,H+1);
% Mcon = zeros(nU,H);
% Scon = zeros(nU,nU,H);
% M(:,1) = m; S(:,:,1) = s;
% for t = 1:H
%     [m, s, mcon, scon, ccon] = plant.prop(m, s, plant, dynmodel, policy, t);
%     
%     if nargout > 2
%         Mcon(:,t) = mcon;
%         Scon(:,:,t) = scon;
%     end
%         
%     M(:,t+1) = m(end-D+1:end);
%     S(:,:,t+1) = s(end-D+1:end,end-D+1:end);
% end

function [M, S] = my_pred(policy, plant, dynmodel, m, s, H)
%% Code
D = length(m);
S = zeros(D,D,H+1); M = zeros(D,H+1);
M(:,1) = m; S(:,:,1) = s;
for t = 1:H
    [m, s] = plant.prop(m, s, plant, dynmodel, policy, t);
    M(:,t+1) = m(end-D+1:end);
    S(:,:,t+1) = s(end-D+1:end,end-D+1:end);
end