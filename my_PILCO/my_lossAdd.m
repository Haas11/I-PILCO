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

function [L, dLdm, dLds, S, dSdm, dSds, C, dCdm, dCds] = my_lossAdd(cost, m, s, varargin)
%% Code

% Dimensions and Initializations
Nlos = length(cost.sub);                D = length(m);
L = 0;              S = 0;              C = zeros(D,1);
dLdm = zeros(1,D);  dSdm = zeros(1,D);  dLds = zeros(1,D^2);
dCdm = zeros(D);    dCds = zeros(D,D^2);%dSds = zeros(D);

dSds = zeros(1,D^2);
dSds_matrix = zeros(D,D);
dLds_matrix = zeros(D,D);
for n = 1:Nlos                            % Loop over each of the sub-functions
    costi = cost.sub{n}; i = costi.losi;    % slice
    
    % Call the sub loss function
    if nargout < 4                      % Just the expected loss & derivs           (happens in my_value.m)
        [Li, Ldm, Lds] = costi.fcn(costi, m(i), s(i,i));        
        L = L + Li;
        dLdm(i) = dLdm(i) + Ldm;
        
        dLds_matrix(i,i) = dLds_matrix(i,i) + reshape(Lds,length(i),length(i));
        dLds = reshape(dLds_matrix,1,D^2);                
    else                                % Also loss variance                        (happens in calcCost)
        [Li, Ldm, Lds, Si, Sdm, Sds, Ci, Cdm, Cds] = costi.fcn(costi, m(i), s(i,i));
        
        L = L + Li;
        S = S + Si + Ci'*s(i,:)*C + C'*s(:,i)*Ci;   % V(a+b) = V(a)+V(b)+C(a,b)+C(b,a)                
        
        if nargout > 4    % derivatives for exploration term             TODO: Does this ever happen? only for added losses (scenario loss functions only have 4 outputs max)
            dLdm(i) = dLdm(i) + Ldm;
            
            dLds_matrix(i,i) = dLds_matrix(i,i) + reshape(Lds,length(i),length(i));
            dLds = reshape(dLds_matrix,1,D^2);
                        
            Cis = Ci'*(s(i,:) + s(:,i)');
            Cs = C'*(s(:,i) + s(i,:)');
            
            dSdm(i) = dSdm(i) + Sdm + Cs*Cdm; 
            dSdm    = dSdm + Cis*dCdm;
            
            dSds_matrix(i,i) = dSds_matrix(i,i) + reshape(Sds,length(i),length(i))...
                + reshape(Cs*Cds,length(i),length(i));
            dSds_matrix = dSds_matrix + reshape(Cis*dCds,D,D);
            dSds_matrix(i,:) = dSds_matrix(i,:) + Ci*C';
            dSds_matrix(:,i) = dSds_matrix(:,i) + C*Ci';
            dSds = reshape(dSds_matrix,1,D^2);
        end
        
        % Input - Output covariance update
        C(i) = C(i) + Ci;                       % must be after S and its derivatives
        
        ii = sub2ind2(D,i,i);
        dCdm(i,i) = dCdm(i,i) + Cdm;             % must be after dSdm & dSds
        dCds(i,ii) = dCds(i,ii) + Cds;
    end
end

% =========================================================================
% Exploration if required
if isfield(cost,'expl') && cost.expl ~= 0 && abs(S)>1e-12 && nargout > 3 % why this nargout? TODO: otherwise double? exploration also happens in lossSat_xxx-->
    L = L + cost.expl*sqrt(S);
    dLdm = dLdm + cost.expl*0.5/sqrt(S)*dSdm;
    dLds = dLds + cost.expl*0.5/sqrt(S)*dSds;
end
% =========================================================================
% Energy penalty if required
if isfield(cost,'ep') && cost.ep ~= 0
    
    % Control mean:
    ma = varargin{1};                                   % control mean          [1 x nU]
    policy = varargin{2};
    
    if isfield(cost,'refPen') && cost.refPen
        % also penalize reference changes
        idx = [policy.impIdx, policy.refIdx];
        translVec = [ones(size(policy.impIdx)).*2, ones(size(policy.refIdx))];
    else
        % only penalize stiffness magnitude
        idx = policy.impIdx;
        translVec = ones(size(policy.impIdx)).*2;
    end
        
    Ma = ma(idx);                                       % select impedance part
    Ma_norm = Ma./(policy.maxU(idx).*translVec)';		% normalized actions    [nUi x 1]
    
    if cost.epType == 2 || (isfield(cost,'refPen') && cost.refPen)
        La = cost.ep*(1/2)*(Ma_norm'*Ma_norm);                    % mean squared energy penalty [1 x 1]
    elseif cost.epType == 1
        La = cost.ep*sum(Ma_norm,1);                        % 1-norm energy penalty 	  [1 x 1]
    else
        error('Invalid type specified for energy penalty. (choose cost.epType=1 for 1-norm, cost.epType=2 for mean squared)');
    end        
    L = L + La;
    
    if nargout > 1 && nargin > 5        % energy penalty w/ derivatives
        dmadm = varargin{3};                            % deriv of control mean w.r.t. policy input state mean [nU x lpoli] 		
        dmads = varargin{4};                            % deriv of control mean w.r.t. policy input state variance [nU x lpoli^2] 	
                
        plant = varargin{5};
        dyno  = plant.dyno;     poli = plant.poli;
        ldyno = length(dyno);   lpoli = length(poli);
        nU = length(policy.maxU);       
        
        if cost.epType == 2 || (isfield(cost,'refPen') && cost.refPen)
            dLadma = Ma_norm';                                          % derivative for mean squared cost	[1 x nUi]
        elseif cost.epType == 1
            dLadma = ones(1,length(idx));                               % derivative of 1 norm 	[1 x nUi]
        end
                
        % deriv control mean w.r.t state mean
        dmadm_con = [dmadm, zeros(nU,ldyno-lpoli)];                         % concat dif with zeros	[nU  x ldyno]  --> index dependent? !!!
        dMadm = dmadm_con(idx,:);                                 % select impedance part [nUi x ldyno]        
        dMadm = bsxfun(@rdivide,dMadm,(policy.maxU(idx)'.*translVec'));    % normalize
        
        dLdm = dLdm + cost.ep*dLadma*dMadm;                                 % [1 x ldyno] 	= [1 x ldyno] 	+ [1 x 1][1 x nUi][nUi x ldyno]
                
        % deriv control mean w.r.t. state variance             
        dynotemp = 1:1:ldyno;
        comp = ismember(dynotemp,poli);                                     % compare inputs 		[1 x ldyno]         
        relevantIdx = find(comp==1);                                        % relevant indices  	[1 x lpoli]             (only consider derivative of states that are inputs to the policy)
        dmads3D = reshape(dmads,[nU lpoli lpoli]);                          % reshape derivative of conrol mean w.r.t. policy input state covariances to 3D matrix
        
        dmads_con = zeros(nU,ldyno,ldyno);                                  % 3D-matrix of derivative of control mean w.r.t. all GP output state covariances  [nU x ldyno x ldyno]
        dmads_con(:,relevantIdx,relevantIdx) = dmads3D(:,:,:);              % fill this matrix at the input dimensions of the policy
        dmads_con = reshape(dmads_con,[nU ldyno^2]);                        % reshape back to 2D 	[nU x ldyno^2]
        
        dMads = dmads_con(idx,:);                                           % select impedance part [nUi x ldyno^2]                
        dMads = bsxfun(@rdivide,dMads,(policy.maxU(idx)'.*translVec'));     % normalize
        
        dLds = dLds + cost.ep*dLadma*dMads;                                 % [1 x ldyno^2] = [1 x ldyno^2] + [1 x 1][1 x nUi][nUi x ldyno^2]]
    end
end
% =========================================================================

function idx = sub2ind2(D,i,j)
% D = #rows, i = row subscript, j = column subscript
i = i(:); j = j(:)';
idx =  reshape(bsxfun(@plus,D*(j-1),i),1,[]);
