%% propagated.m
% *Summary:* Propagate the state distribution one time step forward
%            with derivatives
%
%  function [Mnext, Snext, dMdm, dSdm, dMds, dSds, dMdp, dSdp] = ...
%    propagated(m, s, plant, dynmodel, policy)
%
% *Input arguments:*
%
%   m                 mean of the state distribution at time t           [D x 1]
%   s                 covariance of the state distribution at time t     [D x D]
%   plant             plant structure
%   dynmodel          dynamics model structure
%   policy            policy structure
%
% *Output arguments:*
%
%   Mnext             predicted mean at time t+1                         [E x 1]
%   Snext             predicted covariance at time t+1                   [E x E]
%   dMdm              output mean wrt input mean                         [E x D]
%   dMds              output mean wrt input covariance matrix         [E  x D*D]
%   dSdm              output covariance matrix wrt input mean        [E*E x  D ]
%   dSds              output cov wrt input cov                       [E*E x D*D]
%   dMdp              output mean wrt policy parameters                  [E x P]
%   dSdp              output covariance matrix wrt policy parameters  [E*E x  P]
%
%   where P is the number of policy parameters.
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, Henrik Ohlsson,
% and Carl Edward Rasmussen.
%
% Last modified: 2013-01-23
%
%% High-Level Steps
% # Augment state distribution with trigonometric functions
% # Compute distribution of the control signal
% # Compute dynamics-GP prediction
% # Compute distribution of the next state
%

function [Mnext, Snext, dMdm, dSdm, dMds, dSds, dMdp, dSdp, Mcon, Scon] = ...
    my_propagated(m, s, plant, dynmodel, policy, varargin)
%% Code
global diffChecks diffTol gpCheck REF_PRIOR FIRST
GPgradTypes = {'dMdm','dMds',       'dSdm','dSds',       'dVdm','dVds'};

% persistent prevRef
% if isempty(prevRef)
%         prevRef = zeros(1,length(policy.refIdx))';
% end
persistent prevRefVel prevRefPos
if nargin == 6
    t = varargin{1};
else
    t = 0;
end

if isempty(prevRefVel) || t==1  || FIRST
        prevRefVel = [zeros(1,length(policy.refIdx)), 0]';
        prevRefPos = [m(1:length(policy.refIdx)); 0];
        FIRST = 0;
end
% if isempty(prevRefPos) || t==1
% %         prevRefPos = zeros(1,length(policy.refIdx))';
%         prevRefPos = m(1:length(policy.refIdx));
% end

if nargout == 10                                  % just predict, no derivatives
    [Mnext, Snext, Mcon, Scon] = my_propagate(m, s, plant, dynmodel, policy, t);
    dMdm=[]; dSdm=[]; dMds=[]; dSds=[]; dMdp=[]; dSdp=[];
    return
end

angi = plant.angi; 
poli = plant.poli; 
dyni = plant.dyni; 
difi = plant.difi;
refi = plant.refi;
if isfield(dynmodel,'ref')
    ref = dynmodel.ref;
end

D0 = length(m);                                        % size of the input mean
D1 = D0 + 2*length(angi);          % length after mapping all angles to sin/cos
D2 = D1 + length(policy.maxU);          % length after computing control signal
D3 = D2 + D0;                                         % length after predicting
M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s;   % init M and S

Mdm = [eye(D0); zeros(D3-D0,D0)]; Sdm = zeros(D3*D3,D0);
Mds = zeros(D3,D0*D0); Sds = kron(Mdm,Mdm);
X = reshape(1:D3*D3,[D3 D3]); XT = X'; Sds = (Sds + Sds(XT(:),:))/2;
X = reshape(1:D0*D0,[D0 D0]); XT = X'; Sds = (Sds + Sds(:,XT(:)))/2;

% 1) Augment state distribution with trigonometric functions ------------------
i = 1:D0; j = 1:D0; k = D0+1:D1;
[M(k) S(k,k) C mdm sdm Cdm mds sds Cds] = gTrig(M(i), S(i,i), angi);

[S Mdm Mds Sdm Sds] = ...
    fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,[ ],[ ],[ ],i,j,k,D3);

sn2 = exp(2*dynmodel.hyp(end,:)); sn2(difi) = sn2(difi)/2;
mm=zeros(D1,1); mm(i)=M(i); ss(i,i)=S(i,i)+diag(sn2);
[mm(k), ss(k,k) C] = gTrig(mm(i), ss(i,i), angi);     % noisy state measurement
q = ss(j,i)*C; ss(j,k) = q; ss(k,j) = q';

% 2) Compute distribution of the control signal -------------------------------
i = poli; j = 1:D1; k = D1+1:D2;
[M(k) S(k,k) C mdm sdm Cdm mds sds Cds Mdp Sdp Cdp] = ...
    policy.fcn(policy, mm(i), ss(i,i));

[S Mdm Mds Sdm Sds Mdp Sdp] = ...
    fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,Cdp,i,j,k,D3);

% 3) Compute distribution of the change in state ------------------------------
ii = [dyni D1+1:D2]; j = 1:D2;
if isfield(dynmodel,'sub'), Nf = length(dynmodel.sub); else Nf = 1; end
for n=1:Nf                               % potentially multiple dynamics models
    [dyn i k] = sliceModel(dynmodel,n,ii,D1,D2,D3); j = setdiff(j,k);
    
    [M(k) S(k,k) C mdm sdm Cdm mds sds Cds] = dyn.fcn(dyn, M(i), S(i,i));
    
    if diffChecks && gpCheck
        for type=1:6
            [dd,~,~] = gpT(GPgradTypes{type}, dynmodel, M(i), S(i,i));    % dynamical model gradients
            failedIdx = find(dd>diffTol)';
            if ~isempty(failedIdx)
                fprintf('\nGP derivative check for %4s failed!\n',...
                    GPgradTypes{type});
                fprintf('Failed output dimensions:'); disp(failedIdx);
                %               fprintf('Timestep = %1i \n', h);
                fprintf('Differences are:'); disp(dd(failedIdx)');
            end
        end
    end
    
    [S Mdm Mds Sdm Sds Mdp Sdp] = ...
        fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,Mdp,Sdp,[ ],i,j,k,D3);
    
    j = [j k];                                   % update 'previous' state vector
end

% 4) Compute distribution of the next state -----------------------------------
P = [zeros(D0,D2) eye(D0)]; P(difi,difi) = eye(length(difi));  P = sparse( P);

% Mean:
Mnext = P*M; 
% if REF_PRIOR
%     if ~isfield(policy,'refIdx') || isempty(policy.refIdx);
%         Mnext = Mnext + REF_DIFF(t,:)';                                     % recenter around reference
%         
%     elseif ~isempty(policy.refIdx)
%         ma = M(D1+1:D2);
%         deltaRef_pos = ma(policy.refIdx,1);                                      % policy reference positions @t+1
%         
%         deltaRef_vel = (deltaRef_pos - prevRef)./plant.dt;                       % computed change in reference velocities
%         deltaRef = [deltaRef_pos; 
%                     deltaRef_vel; 
%                     zeros(mod(length(plant.dyno),2),1)];    % (zeros for force indices)
%         
%         Mnext = Mnext + deltaRef;                                               % Recentered next state
%         
%         prevRef = deltaRef_pos;                                                 % Bookkeeping
%     end
% end
if REF_PRIOR && ~isempty(refi)
    if ~isfield(policy,'refIdx') || isempty(policy.refIdx);
        Mnext(refi) = Mnext(refi) + ref(t,refi)';                                     % recenter around pre-computed reference
        
    elseif ~isempty(policy.refIdx)
        % Mean:
        ma = M(D1+1:D2);                                                    % computed actions @t
        deltaRef_pos = [ma(policy.refIdx,1); 0];                                 % action = change in postion reference = ref_t+1 - ref_t
        
        % absolute references:
        curRefPos = prevRefPos + deltaRef_pos;
        curRefVel = deltaRef_pos/plant.dt; % reference velocity = change in pos ref / sampling time        
        curRef = [curRefPos; curRefVel; 0];
        
        % relative differences
        deltaRef_vel = curRefVel - prevRefVel;        
        deltaRef = [deltaRef_pos;                                           
            deltaRef_vel;
            zeros(mod(length(plant.dyno),2),1)];    % (zeros for force indices)
        
        % complete reference:
        finalRef = curRef;
        finalRef(difi) = deltaRef(difi);
        
        Mnext(refi) = Mnext(refi) + finalRef(refi);                                          % Recentered next state
        
        prevRefVel = curRefVel;     % bookkkeeping
        prevRefPos = curRefPos;
    end
end

% Variance:
% Sa = S(D1+1:D2,D1+1:D2);
% SdeltaRef_pos = Sa(policy.refIdx, policy.refIdx);
% Snext = P*S*P' + SdeltaRef_pos; 

Snext = P*S*P'; 
Snext = (Snext+Snext')/2;


PP = kron(P,P);
dMdm =  P*Mdm; dMds =  P*Mds; dMdp =  P*Mdp;
dSdm = PP*Sdm; dSds = PP*Sds; dSdp = PP*Sdp;

X = reshape(1:D0*D0,[D0 D0]); XT = X';                                      % symmetrize dS
dSdm = (dSdm + dSdm(XT(:),:))/2; dMds = (dMds + dMds(:,XT(:)))/2;
dSds = (dSds + dSds(XT(:),:))/2; dSds = (dSds + dSds(:,XT(:)))/2;
dSdp = (dSdp + dSdp(XT(:),:))/2;


% A1) Separate multiple dynamics models ---------------------------------------
function [dyn i k] = sliceModel(dynmodel,n,ii,D1,D2,D3) % separate sub-dynamics
if isfield(dynmodel,'sub')
    dyn = dynmodel.sub{n}; do = dyn.dyno; D = length(ii)+D1-D2;
    if isfield(dyn,'dyni'), di=dyn.dyni; else di=[]; end
    if isfield(dyn,'dynu'), du=dyn.dynu; else du=[]; end
    if isfield(dyn,'dynj'), dj=dyn.dynj; else dj=[]; end
    i = [ii(di) D1+du D2+dj]; k = D2+do;
    dyn.inputs = [dynmodel.inputs(:,[di D+du]) dynmodel.target(:,dj)];   % inputs
    dyn.target = dynmodel.target(:,do);                                 % targets
else
    dyn = dynmodel; k = D2+1:D3; i = ii;
end

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
