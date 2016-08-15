
%% link3_learnImp.m
% *Summary:* Main file for learning planar Peg-in-Hole insertion task.
% Probabilistic model-based RL,
%
%
%
%% High-Level Steps
% # Load parameters
% # Create J initial trajectories by applying random controls
% # Initialise policy with hand-picked targets
% # Controlled learning (train dynamics model, policy learning, policy
% application)

%% Code
% epp = [0.0001 0.001 0.01 0.1 1];
% expll = [-0.1 -0.2 -0.4 -0.6 -0.8];

%% 1. Initialization
clear; clc;
figHandles = findobj('Type','figure');
for i=1:length(figHandles);     % clear figures but retain positions
    clf(figHandles(i));
end
settings_S_iiwa;                  % load scenario-specific settings
ep=num2str(cost.ep);
expl=num2str(cost.expl);
basename = strcat('KUKA_',date,'_','conLin-','_ep-0p',ep(strfind(ep,'.')+1:end),...
    '_expl-0p',expl(strfind(expl,'.')+1:end),'-');      % filename used for saving data

printTuningParams

%% 2. Initial rollouts
fprintf('\nPerforming initial rollouts...\n');
for jj = 1:J
    
    [xx, yy, realCost{jj}, latent{jj}, rr] = ...
        my_iiwaRollout(policy, plant, cost, H, Hdes, a_init);
    realAcumCost(jj) = sum(realCost{jj});
    rr = rr(:,ref_select);      % filter out the relevant reference dimensions
    rrr = rr(2:H+1,:);
    rrr(:,difi) = rr(2:H+1,difi) - rr(1:H,difi); % dimensions that are differences
    x = [x; xx]; y = [y; yy];   r = [r rrr];
    
    % determine if rollout successful:
    if size(xx,1) > H-5
        insertSuccess{1}(jj) = 1;  %#ok<*SAGROW> not aborted
        if mean(abs(latent{1}(end-5:end,dyno(1))-(xhole(1)+0.05))) < 0.01;
            insertSuccess{1}(jj) = 2;                % succesful insertion
        end
    end
    
    if insertSuccess{1}(jj) ~= 0;
        REF_DIFF = rrr;
        dynmodel.ref = rrr;
    else
        error('Initial rollout aborted. Complete reference unavailable.');
    end
    
    if plotting.verbosity > 0
        if ~ishandle(6)         % action iterations
            figure(6);
        else
            set(0,'CurrentFigure',6);
        end
        hold on;
        a = xx(:,end-Du+1:end);
        for i=1:Du
            subplot(ceil(Du/sqrt(Du)),ceil(sqrt(Du)),i);
            hold on;
            stairs(1:length(a(:,i)),a(:,i),strcat(colorVec{jj},'--'));
            legend(iterVec{1:jj});
            xlabel('Timestep');     ylabel(actionTitles{i});
        end
        drawnow;
                       
        if plotting.verbosity > 2
            q_sim = latent{jj}(:,1:robot.n);
            if ~ishandle(5)         % robot animation
                figure(5);
            else
                set(0,'CurrentFigure',5);
            end
            clf(5);
            patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
            robot.plot(q_sim);
            
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
end

tempCost = [];
for i=1:J
    tempCost = [tempCost, realCost{i}];
end
trialAcumCost{1} = sum(tempCost,1);
realWorld.mean(1) = mean(trialAcumCost{1},2);
realWorld.std(1) = std(trialAcumCost{1},0,2);   % flag: 0 = n-1, 1=n

if isempty(find(insertSuccess{1}==2,2))   % None Success
    scoreCard(1) = 0;
elseif length(find(insertSuccess{1}==2,2))==J
    scoreCard(1) = 2;                 % All Success
else
    scoreCard(1) = 1;                 % Partial Success
end

%% 3. Controlled learning (N iterations)
fprintf('\nPILCO Learning started\n--------------------------------\n');
jj=J;
for j = 1:N
    my_trainDynModel;       % train (GP) dynamics model
    my_learnPolicy;         % update policy
    my_applyController;     % apply controller to system
    
    disp(['\nControlled trial # ' num2str(j)]);
    disp('Insertion successes: '); disp(insertSuccess);
end
