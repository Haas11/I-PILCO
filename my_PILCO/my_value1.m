
%% value.m
% *Summary:* Compute expected (discounted) cumulative cost for a given (set of) initial
% state distributions
%
%     function [J, dJdp] = value(p, m0, S0, dynmodel, policy, plant, cost, H)
%
% *Input arguments:*
%
%   p            policy parameters chosen by minimize
%   policy       policy structure
%     .fcn       function which implements the policy
%     .p         parameters passed to the policy
%   m0           matrix (D by k) of initial state means
%   S0           covariance matrix (D by D) for initial state
%   dynmodel     dynamics model structure
%   plant        plant structure
%   cost         cost function structure
%     .fcn       function handle to the cost
%     .gamma     discount factor
%  H             length of prediction horizon
%
% *Output arguments:*
%
%  J             expected cumulative (discounted) cost
%  dJdp          (optional) derivative of J wrt the policy parameters
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-03-21
%
%% High-Level Steps
% # Compute distribution of next state
% # Compute corresponding expected immediate cost (discounted)
% # At end of prediction horizon: sum all immediate costs up

function [J, dJdp] = my_multiValue(p, m0, S0, dynmodel, policy, plant, cost, H)
%% Code
K = size(m0,2);     % number of initial states for which is optimized

policy.p = p;            % overwrite policy.p with new parameters from minimize
p = unwrap(policy.p); 
dp = zeros(length(p),K);%0*p;
% m = m0;
S = S0;
% L = zeros(1,H);
L = zeros(size(m0,2),H);
poli = plant.poli;

if nargout <= 1                                         % no derivatives required
    for k=1:K
        m = m0(:,k);
        if isfield(cost,'ep') && cost.ep ~= 0               % non-zero energy penalty
            for t = 1:H
                [m, S] = plant.prop(m, S, plant, dynmodel, policy, t);      % get next state
                [ma, ~, ~] = policy.fcn(policy, m, S);                      % precompute predicted action for cost function
                L(k,t) = cost.gamma^t.*cost.fcn(cost, m, S, ma, policy);      % expected discounted cost  w/ energy penalty
            end
        else
            for t = 1:H
                [m, S] = plant.prop(m, S, plant, dynmodel, policy, t);      % get next state
                L(k,t) = cost.gamma^t.*cost.fcn(cost, m, S);                  % expected discounted cost
            end
        end
    end
    
else                                                                    % otherwise, get derivatives
    
    
    
    for k=1:K
        dmOdp = zeros([size(m0,1), length(p)]);
        dSOdp = zeros([size(m0,1)*size(m0,1), length(p)]);        
        m = m0(:,k);
        
        if ~isfield(cost,'ep') || (isfield(cost,'ep') && cost.ep==0)
            for t = 1:H                                                         % for all time steps in horizon
                
                [m, S, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
                    plant.prop(m, S, plant, dynmodel, policy, t);               % get next state
                
                dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp; 	% (2 a)
                dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp; 	% (2 b)
                
                [L(k,t), dLdm, dLdS] = cost.fcn(cost, m, S);                     % predictive cost (1)
                L(k,t) = cost.gamma^t*L(k,t);                                      % discount
                dp(:,k) = dp(:,k) + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )';
                
                dmOdp = dmdp; dSOdp = dSdp;                                     % bookkeeping
            end
        else
            % =================================================================
            % Cost w/ energy penalty:
            m_in = m; S_in = S;
            for t=1:H                                                           % for all time steps in horizon
                [m_out, S_out, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
                    plant.prop(m_in, S_in, plant, dynmodel, policy, t);         % get next state
                
                dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;                        % (2 a) 	[D x P]
                dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;                        % (2 b) 	[D x P]
                
                [ma, ~, ~, dmadm, ~, ~, dmads, ~, ~, madp, ~, ~] = ...
                    policy.fcn(policy, m_out(poli), S_out(poli,poli));          % (recompute) control mean + derivatives
                Madp  = madp(policy.impIdx,:);                                  % derivative of control mean w.r.t. policy params  [nU x P]
                dLadp = bsxfun(@rdivide,Madp,(policy.maxU(policy.impIdx)'.*2)); % derivative of energy cost  w.r.t. policy params  [nU x P]
                dLadp = cost.ep*sum(dLadp, 1);                                  % sum individual terms 	[1 x P]
                
                [L(k,t), dLdm, dLdS] = cost.fcn(cost, m_out, S_out, ...
                    ma, policy, dmadm, dmads, plant);                               % predictive cost (1) 	dLdm = [1 x D], dLdS = [1 x D^2]
                L(k,t) = cost.gamma^t*L(k,t);                                       % discounted predicted cost
                
                dp(:,k) = dp(:,k) + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )'...    % discounted deriv of state  cost w.r.t. policy hyperparams
                    + cost.gamma^t*(dLadp)';                                    % discounted deriv of energy cost w.r.t. policy params
                
                dmOdp = dmdp; dSOdp = dSdp;                                     % bookkeeping
                m_in = m_out; S_in = S_out;
                % =============================================================
            end
        end
    end
end
LMulti = mean(L,1);
dpMulti = mean(dp,2);

J = sum(LMulti);
dJdp = rewrap(policy.p, dpMulti);                                                % rewrapped hyperparameters [1 x P]
