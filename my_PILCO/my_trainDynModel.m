%% trainDynModel.m
% *Summary:* Script to learn the dynamics model
% GP Model learning. Targets are defined as absolute values, or as difference for
% indices 'difi'. For indices 'refi' the prior mean is taken as the reference
% trajectory. If noisyInputs is TRUE,
% NIGP toolbox is called to do noisy input Gaussian Process learning
%
% Adapted from:
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: April 2016
% Victor van Spaandonk
%
%% High-Level Steps
% # Extract states and controls from x-matrix
% # Define the training inputs and targets of the GP
% # Train the GP

%% Code
global REF_PRIOR
if j>1, preHyp = dynmodel.hyp; end %#ok<*IJCL>

% 1. Train GP dynamics model
Du = length(policy.maxU);
Da = length(plant.angi); % no. of ctrl and angles

% Inputs:
xaug = [x(:,dyno) x(:,end-Du-2*Da+1:end-Du)];     % x augmented with angles
dynmodel.inputs = [xaug(:,dyni) x(:,end-Du+1:end)];     % use dyni and ctrl

% Targets:
state_target = y(:,dyno);                                           % x_{t+1}     [H x nX]
state_target(:,difi) = state_target(:,difi) - x(:,dyno(difi));      % delta_t   [H x nX]
dynmodel.targets = state_target;                                    % prior mean = 0


if REF_PRIOR    % prior mean = reference
%     if ~isfield(policy,'refIdx') || isempty(policy.refIdx)  % when not learning trajectories
%         
% %         if size(ref_target_repeat,1) < size(x,1)
% %             ref_target_repeat = [ref_target_repeat; ref_target]; % if not already done (debug for rerun after pauze)
% %         end
% %         
% %         if j>1
% %             if size(xx,1) ~= H                  % if rollout was terminated prematurely
% %                 difLength = H - size(xx,1);
% %                 ref_target_repeat = ref_target_repeat(1:end-difLength,:);
% %             end   
% %         end
%         
% %         dynmodel.targets(:,refi) = dynmodel.targets(:,refi) - ref_target_repeat(:,refi);   % targets with reference prior mean  
% 
%         
        dynmodel.targets(:,refi) = dynmodel.targets(:,refi) - r(:,refi);   % targets with reference prior mean  
%         
%     else                                                    % stiffness + trajectory learning
%         ref_target_total = zeros((J+j-1)*H,length(dyno));
%         
%         for cc=1:J+j-1
%             ma = x((cc-1)*H+1:cc*H, end-Du+1:end);                     % actions
%             
% %             aref_pos = ma(:,policy.refIdx);                           % policy reference positions @t+1             
% %             for kk=1:length(policy.refIdx)
% %                 aref_vel(:,kk) = gradient(aref_pos(:,kk),dt_pilco);  % computed reference velocities
% %             end
% %             aref = [aref_pos, aref_vel, zeros(H,1)];                 % reference for rollout No.#cc     ( delta[X Y dX dY Fx]_ref )
% %             ref_target_total((cc-1)*H+1:cc*H,:) = aref;
%             
%             RefVel = zeros(H+1,length(policy.refIdx));
%             
%             deltaRef_pos = ma(:,policy.refIdx);                             % action = change in postion reference = ref_t+1 - ref_t                        
%             RefVel(1:H,:) = deltaRef_pos/plant.dt;                          % reference velocity = change in pos ref / sampling time
%             RefVel(H+1,:) = RefVel(end,:);
%             
%             deltaRef_vel = RefVel(2:H+1,:) - RefVel(1:H,:);
%             
%             deltaRef = [deltaRef_pos, deltaRef_vel, zeros(H,mod(length(plant.dyno),2))];    % (zeros for force indices)
%             ref_target_total((cc-1)*H+1:cc*H,:) = deltaRef;
%         end        
%         dynmodel.targets = state_target - ref_target_total;   % targets with reference prior mean  (difference between difference in state and difference in reference)
%     end   
end

D = size(dynmodel.inputs,2);     % no. of inputs
E = size(dynmodel.targets,2);    % no. of targets
if ~isfield(dynmodel,'hyp')  % if we don't pass in hyper-parameters, define them
    dynmodel.hyp = zeros(D+2,E);
    dynmodel.hyp = repmat([log(std(dynmodel.inputs)) 0 -1]',1,E);   % init hyp length scales
    dynmodel.hyp(D+1,:) = log( std(dynmodel.targets) );                        % signal std dev
    dynmodel.hyp(D+2,:) = log( (std(dynmodel.targets)/5) );                     % noise std dev (SNR = 5)
    for kk=1:D
        if std(dynmodel.inputs(:,kk)) < 1e-4      % avoid numerical instability on init
            warning('Standard deviation smaller than %6.2f, setting parameters to 1',1e-4)
            dynmodel.hyp(kk,:) = ones(1,E);
        end
    end
end

%% Train dynamical model:
if noisyInputs
    fprintf('Training w/ noisy inputs...');
    trainMode = 0;  %  (0) derivative of post. mean, (1) include uncertainty in deriv
    
    Nls = 200; % (+) max number of line searches, (-) max number of fcn evals.
    
    % TODO: Make possible for learning sparse approximation!
    
    if j==1     % first iteration of learning
        [inputNoiseModel, fx] = trainNIGP(dynmodel.inputs,dynmodel.targets,Nls,trainMode,dynmodel.hyp);%,log(inputNoiseSTD'));
    else        % use previous estimate of input noise variance as startpoint
        [inputNoiseModel, fx] = trainNIGP(dynmodel.inputs,dynmodel.targets,Nls,trainMode,dynmodel.hyp,lsipn);
    end
    
    logHyp = inputNoiseModel.seard;         % trained log GP kernel hyperparameters
    lsipn = inputNoiseModel.lsipn;          % trained log input noise standard deviations
    fprintf('Learned input noise std: '); disp(exp(lsipn'));
    slopeVal = inputNoiseModel.df2;         % squared slope values at each training point
    correctionTerms = inputNoiseModel.dipK; % corrective noise terms to be added to the diagonal of K
    
    dynmodel.hyp = logHyp;
    dynmodel.nigp = correctionTerms;

else
    dynmodel = dynmodel.train(dynmodel, plant, trainOpt);  %  train dynamics GP with deterministic inputs
end

% display some hyperparameters
postHyp = dynmodel.hyp;
% noise standard deviations
disp(['Learned output noise std : ' num2str(exp(postHyp(end,:)))]);
% signal-to-noise ratios (values > 500 can cause numerical problems)
disp(['SNRs                     : ' num2str(exp(postHyp(end-1,:)-postHyp(end,:)))]);
disp(num2str(exp(postHyp)));            

% if j>1
%     disp('Change in hyperparameters: '  );
%     disp(num2str(postHyp - preHyp));
% end

