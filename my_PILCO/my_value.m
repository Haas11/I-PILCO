
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

function [J, dJdp] = my_value(p, m0, S0, dynmodel, policy, plant, cost, H)
%% Code

policy.p = p;            % overwrite policy.p with new parameters from minimize
p = unwrap(policy.p); dp = 0*p;
m = m0; S = S0; L = zeros(1,H);
poli = plant.poli; angi = plant.angi;

D0 = length(m);                                        % size of the input mean
D1 = D0 + 2*length(angi);          % length after mapping all angles to sin/cos
D2 = D1 + length(policy.maxU);          % length after computing control signal
D3 = D2 + D0;                                         % length after predicting

% m_ext = zeros(D3,1);    M(1:D0) = m; 
% S = zeros(D3);      S(1:D0,1:D0) = S0;   % init M and S


% if isfield(cost,'refPen') && cost.refPen
%     % also penalize reference changes
%     idx = [policy.impIdx, policy.refIdx];
%     translVec = [ones(size(policy.impIdx))*2, ones(size(policy.refIdx))];
% else
%     % only penalize stiffness magnitude
%     idx = policy.impIdx;
%     translVec = ones(size(policy.impIdx)).*2;
% end

if nargout <= 1                                         % no derivatives required
    
    if isfield(cost,'ep') && cost.ep ~= 0               % non-zero energy penalty
        for t = 1:H
            [m, S] = plant.prop(m, S, plant, dynmodel, policy, t);      % get next state
            [ma, ~, ~] = policy.fcn(policy, m, S);                      % precompute predicted action for cost function
            L(t) = cost.gamma^t.*cost.fcn(cost, m, S, ma, policy);      % expected discounted cost  w/ energy penalty
        end
    else
        for t = 1:H
            [m, S] = plant.prop(m, S, plant, dynmodel, policy, t);      % get next state
            L(t) = cost.gamma^t.*cost.fcn(cost, m, S);                  % expected discounted cost
        end
    end
    
else                                                                    % otherwise, get derivatives
    
    dmOdp = zeros([size(m0,1), length(p)]);
    dSOdp = zeros([size(m0,1)*size(m0,1), length(p)]);
    
    if ~isfield(cost,'ep') || (isfield(cost,'ep') && cost.ep==0)
        for t = 1:H                                                         % for all time steps in horizon
            
            [m, S, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
                plant.prop(m, S, plant, dynmodel, policy, t);               % get next state
            
            dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp; 	% (2 a)
            dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp; 	% (2 b)
            
            [L(t), dLdm, dLdS] = cost.fcn(cost, m, S);                     % predictive cost (1) nargin = 3
            L(t) = cost.gamma^t*L(t);                                      % discount
            dp = dp + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )';
            
            dmOdp = dmdp; dSOdp = dSdp;                                     % bookkeeping
        end
    else
        % =================================================================
        % Cost w/ energy penalty:
        m_in = m; S_in = S;
        iT = cost.iT;
        for t=1:H                                                           % for all time steps in horizon            
            % TODO: CONSIDER TO USE THE ACTION AT PREVIOUS TIMESTEP AS INPUT TO ENERGY PENALTY.            
            [m_out, S_out, dmdmO, dSdmO, dmdSO, dSdSO, dmdp, dSdp] = ...
                plant.prop(m_in, S_in, plant, dynmodel, policy, t);         % get next state           
            
            % for bench: ======================
            Mdm = [eye(D0); zeros(D3-D0,D0)];   Sdm = zeros(D3*D3,D0);
            Mds = zeros(D3,D0*D0);              Sds = kron(Mdm,Mdm);
            
            m_ext = zeros(D3,1);    m_ext(1:D0) = m_out; 
            S_ext = zeros(D3);      S_ext(1:D0,1:D0) = S_out;   % init M and S
            
            i = 1:D0; j = 1:D0; k = D0+1:D1;
            [m_ext(k), S_ext(k,k), C, mdm, sdm, Cdm, mds, sds, Cds] = gTrig(m_out, S_out, angi);            
            [S_ext, Mdm, Mds, Sdm, Sds] = ...
                fillIn(S_ext,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,[ ],[ ],[ ],i,j,k,D3);   
            
            dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;                        % (2 a) 	[D x P]     D is non extended
            dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;                        % (2 b) 	[D x P]
            
            i = poli; j = 1:D1; k = D1+1:D2;            
            [m_ext(k), S_ext(k,k), Ca, dmadm, dsadm, Cadm, dmads, dsads, Cads, dmadp, dsadp, Cadp] = ...
                policy.fcn(policy, m_ext(poli), S_ext(poli,poli));          % (recompute) control mean + derivatives       
            ma = m_ext(k); sa = S_ext(k,k);
            
            [S_ext, Mdm, Mds, Sdm, Sds, Mdp, Sdp] = ...
                fillIn(S_ext,Ca,dmadm,dsadm,Cadm,dmads,dsads,Cads,Mdm,Sdm,Mds,Sds,dmadp,dsadp,Cadp, i,j,k,D3);
                        
            % ========================
            
%             dmdp = dmdmO*dmOdp + dmdSO*dSOdp + dmdp;                        % (2 a) 	[D x P]
%             dSdp = dSdmO*dmOdp + dSdSO*dSOdp + dSdp;                        % (2 b) 	[D x P]
%             
%             [ma, sa, ~, dmadm, dsadm, ~, dmads, dsads, ~, dmadp, dsadp, ~] = ...
%                 policy.fcn(policy, m_out(poli), S_out(poli,poli));          % (recompute) control mean + derivatives            
            
            % ========================

            % Partial Derivatives and Cost:
            [L(t), dLdm, dLdS] = cost.fcn(cost, m_ext(1:D0), S_ext(1:D0,1:D0), ...
                ma, sa, dmadm, dmads, dsadm, dsads, plant);                 % predictive cost (1) 	dLdm = [1 x D], dLdS = [1 x D^2]
            L(t) = cost.gamma^t*L(t);                                       % discounted predicted cost
                        
            % Explicit Derivatives:           
            dLadp = cost.ep*(2*ma'*iT*dmadp + reshape(iT,1,[])*dsadp);
            
            dp = dp + cost.gamma^t*( dLdm(:)'*dmdp + dLdS(:)'*dSdp )'...    % discounted deriv of state  cost w.r.t. policy hyperparams
                + cost.gamma^t*(dLadp)';                                    % discounted (direct) deriv of energy cost w.r.t. policy params since it is directly influenced by it.
            
            dmOdp = dmdp; dSOdp = dSdp;                                     % bookkeeping
            m_in = m_out; S_in = S_out;
            % =============================================================
        end
    end
end

J = sum(L);
dJdp = rewrap(policy.p, dp);                                                % rewrapped hyperparameters [1 x P]


% A2) Apply chain rule and fill out cross covariance terms --------------------
function [S Mdm Mds Sdm Sds Mdp Sdp] = ...
  fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,dCdp,i,j,k,D)

if isempty(k), return; end

X = reshape(1:D*D,[D D]); XT = X';                         % vectorized indices
I=0*X; I(i,i)=1; ii=X(I==1)'; I=0*X; I(k,k)=1; kk=X(I==1)';
I=0*X; I(j,i)=1; ji=X(I==1)'; I=0*X; I(j,k)=1; jk=X(I==1)'; kj=XT(I==1)';

Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:);                           % chainrule
Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);
if isempty(dCdp) && nargout > 5
  Mdp(k,:)  = mdm*Mdp(i,:) + mds*Sdp(ii,:);
  Sdp(kk,:) = sdm*Mdp(i,:) + sds*Sdp(ii,:);
  dCdp      = Cdm*Mdp(i,:) + Cds*Sdp(ii,:);
elseif nargout > 5
  aa = length(k); bb = aa^2; cc = numel(C);
  mdp = zeros(D,size(Mdp,2)); sdp = zeros(D*D,size(Mdp,2));
  mdp(k,:)  = reshape(Mdp,aa,[]); Mdp = mdp;
  sdp(kk,:) = reshape(Sdp,bb,[]); Sdp = sdp;
  Cdp       = reshape(dCdp,cc,[]); dCdp = Cdp;
end

q = S(j,i)*C; S(j,k) = q; S(k,j) = q';                           % off-diagonal
SS = kron(eye(length(k)),S(j,i)); CC = kron(C',eye(length(j)));
Sdm(jk,:) = SS*dCdm + CC*Sdm(ji,:); Sdm(kj,:) = Sdm(jk,:);
Sds(jk,:) = SS*dCds + CC*Sds(ji,:); Sds(kj,:) = Sds(jk,:);
if nargout > 5
  Sdp(jk,:) = SS*dCdp + CC*Sdp(ji,:); Sdp(kj,:) = Sdp(jk,:);
end
