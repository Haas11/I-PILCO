%% pendulum_learn.m
% *Summary:* Script to learn a controller for the pendulum swingup
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-27
%
%% High-Level Steps
% # Load parameters
% # Create J initial trajectories by applying random controls
% # Controlled learning (train dynamics model, policy learning, policy
% application)

%% 1. Initialization
clear; close all;
my_settings_pendulum;            % load scenario-specific settings
basename = 'pendulum_';       % filename used for saving data

%% 2. Initial J random rollouts
for jj = 1:J
    [xx, yy, realCost{jj}, latent{jj}] = ...
        bench_rollout(gaussian(mu0, S0), struct('maxU',policy.maxU), H, plant, cost);
    x = [x; xx]; y = [y; yy];       % augment training sets for dynamics model
    
    if mean(abs(latent{1}(end-5:end,dyno(2))-cost.sub{1}.target(2)')) < 0.01;
        insertSuccess{1}(jj) = 2;                % succesful insertion
    end    
    
    if plotting.verbosity > 0;      % visualization of trajectory        
        if plotting.verbosity > 0
            if ~ishandle(6)         % action iterations
                figure(6);
            else
                set(0,'CurrentFigure',6);
            end
            hold on;
            a = xx(:,end-Du+1:end);
            for i=1:Du
                subplot(2,1,i);
                hold on;
                stairs((1:length(a(:,i)))*plant.dt,a(:,i),strcat(colorVec{jj},'--'));
                legend(iterVec{1:jj});
                xlabel('Timestep');     ylabel(actionTitles{i});
            end
            drawnow;
            
            if plotting.verbosity > 2                                
                if ~ishandle(7)         % recorded cost iterations
                    figure(7);
                else
                    set(0,'CurrentFigure',7);
                end
                hold on; grid on; stairs(1:length(realCost{jj}),realCost{jj},strcat(colorVec{jj},'--'));
                legend(iterVec{1:jj});
                title('Recorded Rollout Cost');
                xlabel('Timestep');     ylabel('Immediate Cost');
                drawnow;
            end
        end
        
        if ~ishandle(1); figure(1); else set(0,'CurrentFigure',1); end; clf(1);
        draw_rollout_pendulum;
    end
end

mu0Sim(odei,:) = mu0; S0Sim(odei,odei) = S0;
mu0Sim = mu0Sim(dyno)'; S0Sim = S0Sim(dyno,dyno);

% Compute cost distributions of real world trials:
tempCost = [];
for i=1:J
    tempCost = [tempCost, realCost{i}];
end
trialAcumCost{1} = sum(tempCost,1);
realWorld.mean(1) = mean(trialAcumCost{1},2);
realWorld.std(1) = std(trialAcumCost{1},0,2);   % flag: 0 = n-1, 1=n

if isempty(find(insertSuccess{1}==2,1))   % None Success
    scoreCard(1) = 0;
elseif length(find(insertSuccess{1}==2,1))==J
    scoreCard(1) = 2;                 % All Success
else
    scoreCard(1) = 1;                 % Partial Success
end

jj=J; initTrial = 0;

%% 3. Controlled learning (N iterations)
for j = 1:N
    trainDynModel;   % train (GP) dynamics model
    my_learnPolicy;     % learn policy
    bench_applyController; % apply controller to system
    
    disp(['controlled trial # ' num2str(j)]);
    
    if plotting.verbosity > 0;      % visualization of trajectory
        if ~ishandle(1); figure(1); else set(0,'CurrentFigure',1); end; clf(1);
        draw_rollout_pendulum;
    end
end