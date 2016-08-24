%% my_lossAdd.m
% *Summary:* Utility function to add a number of loss functions together, each of which
% can be using a different loss function and operating on a different part of
% the state.
%
%       function [L, dLdm, dLds, S, dSdm, dSds, C, dCdm, dCds] = my_lossAdd(cost, m, s)
%
% *Input arguments:*
%
%   cost            cost struct
%      .fcn         @lossAdd - called to arrive here
%      .sub{n}      cell array of sub-loss functions to add together
%         .fcn      handle to sub function
%         .losi     indices of variables to be passed to loss function
%         .< >      all fields in sub will be passed onto the sub function
%      .expl        (optional) if present cost.expl*sqrt(S) is added to the loss
%
%    m         mean of input distribution                              [D x 1]
%    s         covariance matrix of input distribution                 [D x D]
%a
% *Output arguments:*
%
%   L               expected loss                                  [1   x   1]
%   dLdm            derivative of L wrt input mean                 [1   x   D]
%   dLds            derivative of L wrt input covariance           [1   x   D^2]
%   S               variance of loss                               [1   x   1]
%   dSdm            derivative of S wrt input mean                 [1   x   D]
%   dSds            derivative of S wrt input covariance           [1   x   D^2]
%   C               inv(S) times input-output covariance           [D   x   1]
%   dCdm            derivative of C wrt input mean                 [D   x   D]
%   dCds            derivative of C wrt input covariance           [D   x   D^2]
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-05

function [L, dLdm, dLds, S, dSdm, dSds, C, dCdm, dCds] = my_lossAdd2(cost, m, s, varargin)
%% Code

% Dimensions and Initializations
Nlos = length(cost.sub);                D = length(m);
L = 0;              S = 0;              C = zeros(D,1);
dLdm = zeros(1,D);  dSdm = zeros(1,D);  dLds = zeros(1,D^2);
dCdm = zeros(D);    dCds = zeros(D,D^2);%dSds = zeros(D);

dSds = zeros(1,D^2);
dSds_matrix = zeros(D,D);
dLds_matrix = zeros(D,D);

for n = 1:Nlos                              % Loop over each of the sub-functions
    costi = cost.sub{n}; i = costi.losi;    % slice
    
    % Call the sub loss function
    if nargout < 3
        [Li, ~, ~] = costi.fcn(costi, m(i), s(i,i));
        
        L = L + Li;
        
    else
        [Li, Ldm, Lds, Si, Sdm, Sds, Ci, Cdm, Cds] = costi.fcn(costi, m(i), s(i,i));
        
        L = L + Li;                                 % loss mean
        S = S + Si + Ci'*s(i,:)*C + C'*s(:,i)*Ci;   % loss variance         V(a+b) = V(a)+V(b)+C(a,b)+C(b,a)    (desired in calcCost)
        
        dLdm(i) = dLdm(i) + Ldm;                    % derivative of mean loss wrt mean state
        
        dLds_matrix(i,i) = dLds_matrix(i,i) + reshape(Lds,length(i),length(i));
        dLds = reshape(dLds_matrix,1,D^2);          % derivative of mean loss wrt variance state
        
        Cis = Ci'*(s(i,:) + s(:,i)');
        Cs = C'*(s(:,i) + s(i,:)');
        
        dSdm(i) = dSdm(i) + Sdm + Cs*Cdm;
        dSdm    = dSdm + Cis*dCdm;                  % derivative of variance loss wrt mean state        --> used for exploration
        
        dSds_matrix(i,i) = dSds_matrix(i,i) + reshape(Sds,length(i),length(i))...
            + reshape(Cs*Cds,length(i),length(i));
        dSds_matrix = dSds_matrix + reshape(Cis*dCds,D,D);
        dSds_matrix(i,:) = dSds_matrix(i,:) + Ci*C';
        dSds_matrix(:,i) = dSds_matrix(:,i) + C*Ci';
        dSds = reshape(dSds_matrix,1,D^2);          % derivative of variance los wrt variance state     --> used for exploration
        
        % Input - Output covariance update
        C(i) = C(i) + Ci;                       % must be after S and its derivatives
        ii = sub2ind2(D,i,i);
        dCdm(i,i) = dCdm(i,i) + Cdm;            % must be after dSdm & dSds
        dCds(i,ii) = dCds(i,ii) + Cds;
    end
end

% =========================================================================
% Exploration if required
if isfield(cost,'expl') && cost.expl ~= 0 && abs(S)>1e-12% && nargout > 3   % during calcCost
    L = L + cost.expl*sqrt(S);
    dLdm = dLdm + cost.expl*0.5/sqrt(S)*dSdm;
    dLds = dLds + cost.expl*0.5/sqrt(S)*dSds;       % does the energy penalty affect the loss variance?
end

% =========================================================================
% Energy penalty if required
if isfield(cost,'ep') && cost.ep ~= 0
    % Control mean:
    ma = varargin{1};                                  % control mean          [1 x nU]
    sa = varargin{2};                                  % control variance      [nU x nU]
    iT = cost.iT;
    
    La = cost.ep*(trace(sa*iT) + ma'*iT*ma);            % mean squared energy penalty
    L = L + La;
    
    % energy penalty w/ derivatives
    if nargout > 1 && nargin > 5
        dmadm = varargin{3};                            % deriv of control mean w.r.t. policy input state mean      [nU x lpoli]
        dmads = varargin{4};                            % deriv of control mean w.r.t. policy input state variance  [nU x lpoli^2]
        dsadm = varargin{5};                            % deriv of control variance w.r.t. policy input mean        [nU^2 x lpoli]
        dsads = varargin{6};                            % deriv of control variance w.r.t. policy input variance    [nU^2 x lpoli^2]
        plant = varargin{7};        
        dyno  = plant.dyno;     poli = plant.poli; 
        ldyno = length(dyno);   nU = length(ma);                
        X = reshape(1:D*D,[D D]);
        i = poli; I=0*X; I(i,i)=1; ii=X(I==1)';         % ii = column numbers corresponding to relevant policy input states
            
        dMadm = zeros(nU,ldyno);
        dMadm(:,poli) = dmadm;
        
        dSadm = zeros(nU^2, ldyno);
        dSadm(:,poli) = dsadm;                                              % nUi*lpoli nonzero indices
        
        % partial derivative of energy cost w.r.t. state means
        dLadm = 2*ma'*iT*dMadm + reshape(iT,1,nU^2)*dSadm;                  % [1 x ldyno] = [1 x nU][nU x nU][nU x ldyno] + [nU x nU][nU^2 x ldyno]
        dLdm = dLdm + cost.ep*(dLadm);                                      % [1 x ldyno]
        
        % partial derivative of energy cost w.r.t. state variance
        dMads = zeros(nU,ldyno^2);
        dMads(:,ii) = dmads;
        
        dSads = zeros(nU^2,ldyno^2);
        dSads(:,ii) = dsads;
        
        dLads = 2*ma'*iT*dMads + reshape(iT,1,nU^2)*dSads;                                  %[1 x ldyno^2] = [1 x nU][nU x nU][nU x ldyno^2] + [1 x nU^2][nU^2 x ldyno^2]
        dLds = dLds + cost.ep*dLads;        
    end
end
% =========================================================================

function idx = sub2ind2(D,i,j)
% D = #rows, i = row subscript, j = column subscript
i = i(:); j = j(:)';
idx =  reshape(bsxfun(@plus,D*(j-1),i),1,[]);

% function [S, Mdm, Mds, Sdm, Sds] = ...
%                  fillIn(S,C,mdm,sdm,Cdm,mds,sds,Cds,Mdm,Sdm,Mds,Sds,i,k,D)
% X = reshape(1:D*D,[D D]); XT = X';                    % vectorized indices
% I=0*X; I(i,i)=1; ii=X(I==1)'; I=0*X; I(k,k)=1; kk=X(I==1)';
% I=0*X; I(i,k)=1; ik=X(I==1)'; ki=XT(I==1)';
%
% Mdm(k,:)  = mdm*Mdm(i,:) + mds*Sdm(ii,:);                      % chainrule
% Mds(k,:)  = mdm*Mds(i,:) + mds*Sds(ii,:);
% Sdm(kk,:) = sdm*Mdm(i,:) + sds*Sdm(ii,:);
% Sds(kk,:) = sdm*Mds(i,:) + sds*Sds(ii,:);
% dCdm      = Cdm*Mdm(i,:) + Cds*Sdm(ii,:);
% dCds      = Cdm*Mds(i,:) + Cds*Sds(ii,:);
%
% S(i,k) = S(i,i)*C; S(k,i) = S(i,k)';                        % off-diagonal
% SS = kron(eye(length(k)),S(i,i)); CC = kron(C',eye(length(i)));
% Sdm(ik,:) = SS*dCdm + CC*Sdm(ii,:); Sdm(ki,:) = Sdm(ik,:);
% Sds(ik,:) = SS*dCds + CC*Sds(ii,:); Sds(ki,:) = Sds(ik,:);
