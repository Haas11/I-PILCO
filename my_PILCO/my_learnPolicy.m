%% learnPolicy.m
% *Summary:* Script to perform the policy search
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-06
%
%% High-Level Steps
% # Learn the policy (call optimizer)
% # Predict trajectory from a sampled start state and compute expected cost

%% Code
global gradTypes GPgradTypes diffChecks diffTol gpCheck conCheck propCheck valueCheck satCheck lossCheck checkFailed FIRST%#ok<NUSED>
gradTypes   = {'dMdm','dMds','dMdp','dSdm','dSds','dSdp','dCdm','dCds','dCdp'};
GPgradTypes = {'dMdm','dMds',       'dSdm','dSds',       'dVdm','dVds'};

if diffChecks && valueCheck
    fprintf('\nPerforming value function derivative checks...\n');
    fprintf('===============================================\n\n');
    [dd,~,~] = valueT(policy.p, mu0Sim', S0Sim, dynmodel, policy, plant, cost, 10);
    failedIdx = find(dd>diffTol)';
    if ~isempty(failedIdx)
        fprintf('\nDerivative check for dJdp failed!\n');
        fprintf('Failed output dimensions:'); disp(failedIdx);
        fprintf('Difference is:\t'); disp(dd(failedIdx)');
        checkFailed = 1;
    end
end
gpCheck = 0;

if isfield(policy.p,'hyp')
    fprintf('\nPolicy hypers:');
    disp(exp(policy.p.hyp));
end

% 1. Update the policy
if K==1
    % optimize for single start state
    [policy.p, fX3] = minimize(policy.p, 'my_value', opt, mu0Sim', S0Sim, ...
        dynmodel, policy, plant, cost, H);
    
    % 2. Predict state trajectory from p(x0) and compute cost trajectory
    [M{j}, Sigma{j}, Mcon{j}, Scon{j}] = my_pred(policy, plant, dynmodel, mu0Sim', S0Sim, H);
    [fantasy.mean{j}, fantasy.std{j}] =  my_calcCost(cost, M{j}, Sigma{j}, policy, plant);        
else
    % randomize multiple start states
    mu0MultiSim = mu0Sim';
    for k=1:K
        mu0MultiSim(:,k+1) = mu0MultiSim(:,1)+((2*rand(size(startStateInterval))-1).*startStateInterval);
    end
    
    [policy.p, fX3] = minimize(policy.p, 'my_multiValue', opt, mu0MultiSim, S0Sim, ...
        dynmodel, policy, plant, cost, H);
    
    % 2. Predict state trajectories and cost:
    [fantasy.mean{j}, fantasy.std{j}, M_multi{j}, S_multi{j}, Mcon_multi{j}, Scon_multi{j}] = my_predcost(mu0MultiSim, S0Sim, dynmodel, plant, policy, cost, H);
    
    % nominal trajectories for plotting in figure(4);
    M{j} = M_multi{j}{1};
    Sigma{j} = S_multi{j}{1};
    
    Mcon{j} = Mcon_multi{j}{1};
    Scon{j} = Scon_multi{j}{1};    
end

if compareToFullModel && numel(dynmodel.induce) ~= 0
    dynmodel.full = true;   % force full model predictions (gp0)
    [Mfull{j}, Sfull{j}] = my_pred(policy, plant, dynmodel, mu0Sim', S0Sim, H);     % compute trajectory based on full GP model
    dynmodel.full = false;  % revert to sparse for next iter. (gp1d)
    if j==25
        compareToFullModel = false;     % stop comparison after 25 trials (saves time)
    end
end

% Fantasy World:
fantasy.acumMean(j) = sum(fantasy.mean{j});
fantasy.acumStd(j) = sum(fantasy.std{j});

%% Derivative Checks [DEBUG]
if diffChecks       % perform derivate checks [DEBUG]
    if lossCheck
        lossVec = {'dMdm','dMds','dSds','dCdm','dCds'};
        for h=1:H
            for k=1:5
                [dd,~,~] = lossT(lossVec{k}, cost, M{j}(:,h), Sigma{j}(:,:,h));     % cost gradients    Doesn't work w/ energy penalty!
                failedIdx = find(dd>diffTol)';
                if ~isempty(failedIdx)
                    fprintf('\nController derivative check for %4s failed!\n',...
                        lossVec{kk});
                    fprintf('Timestep = %1i \n', h);
                    fprintf('Difference is:'); disp(dd(failedIdx)');
                    checkFailed = 1;
                end
            end
        end
    end
    if propCheck
        FIRST = true;
        for h=1:H
            for k=1:6
                [dd,~,~] = propagateT(gradTypes{k}, plant, ...
                    dynmodel, policy, M{j}(:,h), Sigma{j}(:,:,h));                    % state propagation gradients
                failedIdx = find(dd>diffTol)';
                if ~isempty(failedIdx)
                    fprintf('\nPropagate derivative check for %4s failed!\n',...
                        gradTypes{k});
                    fprintf('Failed output dimensions:'); disp(failedIdx);
                    fprintf('Timestep = %1i \n', h);
                    fprintf('Differences are:'); disp(dd(failedIdx)');
                    checkFailed = 1;
                end
            end
        end
    end
    if conCheck
        for h=1:H
            for k=1:9
                [dd,~,~] = conT(gradTypes{k}, policy, M{j}(:,h), Sigma{j}(:,:,h));    % controller gradients
                failedIdx = find(dd>diffTol)';
                if ~isempty(failedIdx)
                    fprintf('\nController derivative check for %4s failed!\n',...
                        gradTypes{k});
                    fprintf('Failed output dimensions:'); disp(failedIdx);
                    fprintf('Timestep = %1i \n', h);
                    fprintf('Differences are:'); disp(dd(failedIdx)');
                    checkFailed = 1;
                end
            end
        end
    end
end

%% Verbosity
if plotting.verbosity > 0
    if ~ishandle(3)                     % predicted immediate costs & variances
        figure(3)
    else
        set(0,'CurrentFigure',3);
    end
<<<<<<< HEAD
    clf(3); 
    errorbar((0:H)*plant.dt,fantasy.mean{j},2*fantasy.std{j},'r');
=======
    clf(3);
    errorbar((0:H)*plant.dt,fantasy.mean{j},2*fantasy.std{j});
>>>>>>> d5abe78d437e0f83a378b815e55750f3b2c4c9c6
    xlabel('time step'); ylabel('immediate cost');
    drawnow;
    
    % Plot overall optimization progress
    if plotting.verbosity > 1
        if ~ishandle(2)
            figure(2)
        else
            set(0,'CurrentFigure',2);
        end
        if j>1 %#ok<*IJCL>   plot horizontal line to indicate iteration count
            hold on;
            ylim = get(gca,'YLim');
            line([prevLength prevLength],[ylim(1),max(fX3)],'Color',[0.5 0.5 0.5],'LineStyle',':');
        end
        hold on; plot(1:prevLength+length(fX3),[NaN(prevLength,1); fX3]');
        xlabel('line search iteration'); ylabel('function value')
        title('Accumulated Cost per Policy Update');
        drawnow;
        prevLength = prevLength + length(fX3);
    end
end