% *Input arguments:*
%
%   cost
%     .fcn      @lossHinge - called to get here
%     .a        slope of loss function
%     .b        corner points of loss function                     [1   x    2 ]
%   m             input mean                                       [D   x    1 ]
%   S             input covariance matrix                          [D   x    D ]
%
%   D               = length(dyno)
%
% *Output arguments:*
%
%   L     expected cost                                             [1 x  1 ]
%   dLdm  derivative of expected cost wrt. state mean vector        [1 x  D ]
%   dLds  derivative of expected cost wrt. state covariance matrix  [1 x D^2]
%   S2    variance of cost                                          [1 x  1 ]
%
%   ???????? --> add these variances
%   dSdm            derivative of S wrt input mean                 [1   x   D]
%   dSds            derivative of S wrt input covariance           [1   x   D^2]
%   C               inv(S) times input-output covariance           [D   x   1]
%   dCdm            derivative of C wrt input mean                 [D   x   D]
%   dCds            derivative of C wrt input covariance           [D   x   D^2]
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-08
%
%% High-Level Steps
% # Precomputations
% # Define static penalty as distance from target setpoint
% # Trigonometric augmentation
% # Calculate loss

function [L, dLdm, dLds, S2, dSdm, dSds, C2, dCdm, dCds] = my_lossHinge(cost, m, s)
%% Code
if isfield(cost,'width')
    cw = cost.width;
else
    cw = 1;
end
if ~isfield(cost,'expl') || isempty(cost.expl)
    b = 0;
else
    b =  cost.expl;
end

% 1. Some precomputations
D0 = size(s,2);                           % state dimension

M = zeros(D0,1);    M(1:D0) = m;
S = zeros(D0);      S(1:D0,1:D0) = s;

Mdm = [eye(D0); zeros(D0-D0,D0)];
Sdm = zeros(D0*D0,D0);
Mds = zeros(D0,D0*D0);
Sds = kron(Mdm,Mdm);

% Calculate loss
L    = 0;
dLdm = zeros(1,D0);
dLds = zeros(1,D0*D0);

S2   = 0;
dSds = zeros(1,D0*D0);
dSdm = zeros(1,D0);

C2   = zeros(D0,1);
dCdm = zeros(D0,D0);
dCds = zeros(D0,D0*D0);

if nargout < 5                                                              %(happens in value.m)
    [r, rdM, rdS, s2, s2dM, s2dS]                 = lossHinge(cost, M, S);
else                                                                        %(happens in calcCost.m for loss variance)
    [r, rdM, rdS, s2, s2dM, s2dS, c2, c2dM, c2dS] = lossHinge(cost, M, S);
    % scale mixture of IO variances and derivs
    dSdm = dSdm + s2dM;
    dSds = dSds + reshape(s2dS,1,D0^2);
    % IS SIMPLE ADDITION THE CORRECT WAY TO SCALE MIXTURES HERE?!?!
    % FOR 1 LENGTH SCALE IT SHOULD MAKE NO DIFFERENCE.
    
    % TODO: tested multiple length scales -> results in negative
    % IO variances! how?
    C2 = C2 + c2;
    dCdm = dCdm + c2dM;
    dCds = dCds + c2dS;
end

L = L + r;
S2 = S2 + s2;
dLdm = dLdm + rdM(:)'*Mdm + rdS(:)'*Sdm;
dLds = dLds + rdM(:)'*Mds + rdS(:)'*Sds;

if (b~=0 || ~isempty(b)) && abs(s2)>1e-12
    L = L + b*sqrt(s2);
    dLdm = dLdm + b/sqrt(s2) * ( s2dM(:)'*Mdm + s2dS(:)'*Sdm )/2;
    dLds = dLds + b/sqrt(s2) * ( s2dM(:)'*Mds + s2dS(:)'*Sds )/2;
end

% normalize
n = length(cw);
L = L/n;
dLdm = dLdm/n;
dLds = dLds/n;
S2 = S2/n;