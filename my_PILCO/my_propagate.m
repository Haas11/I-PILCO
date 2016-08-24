%% propagate.m
% *Summary:* Propagate the state distribution one time step forward.
%
%  [Mnext, Snext] = propagate(m, s, plant, dynmodel, policy)
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
%   Mnext             mean of the successor state at time t+1            [E x 1]
%   Snext             covariance of the successor state at time t+1      [E x E]
%
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

function [Mnext, Snext, Mcon, Scon] = my_propagate(m, s, plant, dynmodel, policy, varargin)
%% Code
global diffChecks diffTol gpCheck REF_PRIOR TIME_INPUT
persistent prevRefVel prevRefPos

if nargin == 6
    t = varargin{1};
else
    t = 0;
end

if TIME_INPUT
    m = [m; t];
    s = [s, zeros(length(s),1)];
    s = [s; zeros(1,length(s))];
end

% extract important indices from structures
angi = plant.angi;  % angular indices
poli = plant.poli;  % policy indices
dyni = plant.dyni;  % dynamics-model indices
difi = plant.difi;  % state indices where the model was trained on differences
refi = plant.refi;  % indices for which to encode a prior reference
if isfield(dynmodel,'ref')
    ref = dynmodel.ref;
end

D0 = length(m);                                        % size of the input mean
D1 = D0 + 2*length(angi);          % length after mapping all angles to sin/cos
D2 = D1 + length(policy.maxU);          % length after computing control signal
D3 = D2 + D0;                                         % length after predicting
M = zeros(D3,1); M(1:D0) = m; S = zeros(D3); S(1:D0,1:D0) = s;   % init M and S

% 1) Augment state distribution with trigonometric functions ------------------
i = 1:D0; j = 1:D0; k = D0+1:D1;
[M(k), S(k,k) C] = gTrig(M(i), S(i,i), angi);
q = S(j,i)*C; S(j,k) = q; S(k,j) = q';

sn2 = exp(2*dynmodel.hyp(end,:)); sn2(difi) = sn2(difi)/2;
mm=zeros(D1,1); mm(i)=M(i); ss(i,i)=S(i,i)+diag(sn2);
[mm(k), ss(k,k) C] = gTrig(mm(i), ss(i,i), angi);     % noisy state measurement
q = ss(j,i)*C; ss(j,k) = q; ss(k,j) = q';

% 2) Compute distribution of the control signal -------------------------------
i = poli; j = 1:D1; k = D1+1:D2;
[M(k) S(k,k) C] = policy.fcn(policy, mm(i), ss(i,i));
q = S(j,i)*C; S(j,k) = q; S(k,j) = q';

Mcon = M(k);
Scon = S(k,k);
% [mTest, sTest, cTest]  = policy.fcn(policy, mm(i), zeros(length(i)));
% disp(['control mean for uncertain state: ' num2str(M(k)')])
% disp(['control mean for certain state: ' num2str(mTest')]);
% fprintf('\n');
% disp('control var for uncertain state: '); disp(S(k,k)');
% disp('control var for certain state: '); disp(sTest');
% fprintf('\n\n');


% 3) Compute dynamics-GP prediction              ------------------------------
ii = [dyni D1+1:D2]; j = 1:D2;
if isfield(dynmodel,'sub'), Nf = length(dynmodel.sub); else Nf = 1; end
for n=1:Nf                               % potentially multiple dynamics models
    [dyn i k] = sliceModel(dynmodel,n,ii,D1,D2,D3); j = setdiff(j,k);
    [M(k), S(k,k), C] = dyn.fcn(dyn, M(i), S(i,i));
    
    if diffChecks && gpCheck
        gradTypes = {'dMdm','dMds','dSdm','dSds','dVdm','dVds'};
        for type=1:6
            [dd,~,~] = gpT(gradTypes{type}, dynmodel, M(i), S(i,i));    % dynamical model gradients
            failedIdx = find(dd>diffTol)';
            if ~isempty(failedIdx)
                fprintf('\GP derivative check for %4s failed!\n',...
                    gradTypes{type});
                fprintf('Failed output dimensions:'); disp(failedIdx);
                fprintf('Timestep = %1i \n', t);
                fprintf('Differences are:'); disp(dd(failedIdx)');
            end
        end
    end
    
    q = S(j,i)*C; S(j,k) = q; S(k,j) = q';
    
    j = [j k];                                   % update 'previous' state vector
end

% 4) Compute distribution of the next state -----------------------------------
P = [zeros(D0,D2) eye(D0)];
P(difi,difi) = eye(length(difi));

Mnext = P*M;

% 5) Encode prior reference mean ---------------------------------------------
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
        
        % Variance:
        % Sa = S(D1+1:D2,D1+1:D2);
        % SdeltaRef_pos = Sa(policy.refIdx, policy.refIdx);
        %
        % posIdx = 1:1:length(policy.refIdx);
        % velIdx = posIdx + length(policy.refIdx);
        %
        % SdeltaRef(posIdx,posIdx) = Sa(policy.refIdx, policy.refIdx);
        % SdeltaRef(velIdx,velIdx) = SdeltaRef_pos.*(1/plant.dt)^2 ;
        % SdeltaRef(D0-(D0-length(posIdx*2)),D0-(D0-length(posIdx*2))) = SdeltaRef_pos;
        %
        % Snext = P*S*P';
    end
end

Snext = P*S*P';
Snext = (Snext+Snext')/2;


function [dyn i k] = sliceModel(dynmodel,n,ii,D1,D2,D3) % separate sub-dynamics
% A1) Separate multiple dynamics models ---------------------------------------
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