%% applyController.m
% *Summary:* Script to apply the learned controller to a (simulated) system
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Edited by Victor van Spaandonk.
%
% Last modified: 2015-03-01
%
%% High-Level Steps
% # Generate a single trajectory rollout by applying the controller
% # Generate many rollouts for testing the performance of the controller
% # Save the data

%% 1. Generate trajectory rollout given the current policy
% TODO: CONSIDER PASSING VARIABLE START STATE BACK AND USING THAT AS IC FOR
% PREDICTED STATE TRAJECTORY IN PREDCOST.M OR PRED.M

if KUKA
    [xx, yy, realCost{j+J}, latent{j+J}, rr] = ...
        my_iiwaRollout(policy, plant, cost, H, Hdes);
else
    [xx, yy, realCost{j+J}, latent{j+J}, rr] = ...
        my_rollout(mu0, policy, H, plant, robot); %#ok<*IJCL>
end
realAcumCost(j+J) = sum(realCost{j+J});                 % Accumulated cost
rr = rr(:,ref_select);  rrr = rr(2:end,:);
rrr(:,difi) = rr(2:end,difi) - rr(1:end-1,difi);    % dimensions that are differences
r = [r; rrr];   x = [x; xx];    y = [y; yy]; %#ok<*AGROW>

% determine if rollout aborted, failed or successful:
lengthDiff = H - size(realCost{j+J},1);
if lengthDiff == 0
    insertSuccess{j+1}(1) = 1;                         %#ok<*SAGROW> not aborted
    if mean(abs(latent{j+J}(end-5:end,dyno(1))-0.75)) < 0.01;
        insertSuccess{j+1}(1) = 2;                     % successful insertion
    end
else
    insertSuccess{j+1}(1) = 0;
end

tempCost = [realCost{J+j}; ones(lengthDiff,1)];

% 2. Make many rollouts to test the controller quality / robustness?
testLati = cell(1,Ntest);
testCosti = cell(1,Ntest);
for i=1:Ntest
    if KUKA
        [~, ~, testCosti{i}, testLati{i}, ~] = ...
            my_iiwaRollout(policy, plant, cost, H, Hdes);
    else
        [~, ~, testCosti{i}, testLati{i}, ~] = ...
            my_rollout(mu0, policy, H, plant, robot); %#ok<*IJCL>
    end
    
    lengthDiff = H-size(testCosti{i},1);
    if lengthDiff==0                           % rollout was completed
        insertSuccess{j+1}(i+1) = 1;
        if mean(abs(testLati{i}(end-5:end,dyno(1))-0.75)) < 0.01;
            insertSuccess{j+1}(i+1) = 2;                     % successful insertion
        end
        
        tempCost = [tempCost, [testCosti{i}; ones(lengthDiff,1)]];          % concatenate with ones if rollout aborted
    end
end
testLat{j} = testLati;
testCost{j} = testCosti;

% Real World:
trialAcumCost{j+1} = sum(tempCost,1);
realWorld.mean(j+1) = mean(trialAcumCost{j+1},2);
realWorld.std(j+1) = std(trialAcumCost{j+1},0,2);   % flag: 0 = n-1, 1=n
if isempty(find(insertSuccess{j+1}==2,2))   % None Success
    scoreCard(j+1) = 0;
elseif all(insertSuccess{j+1}==2)
    scoreCard(j+1) = 2;                 % All Success
else
    scoreCard(j+1) = 1;                 % Partial Success
end

%% verbosity
if plotting.verbosity > 0
    % Cost During Last Trials
    if ~ishandle(3)
        figure(3);
    else
        set(0,'CurrentFigure',3);
    end
    hold on;
    stairs((1:length(realCost{J+j}))*plant.dt,realCost{J+j},'r'); hold on;
    for ii=1:Ntest
        stairs((1:length(testCosti{ii}))*plant.dt,testCosti{ii},'g');
        hold on;
    end
    title('Predicted (uncertain) & Rollout (deterministic) Immediate Cost');
    xlabel('Time [s]');     ylabel('Immediate Cost');
    legend('predicted','test','verifications');
    drawnow;
    
    % Actions over all iterations:
    if ~ishandle(6)
        figure(6);
    else
        set(0,'CurrentFigure',6);
    end
    hold on;
    a = xx(:,end-Du+1:end);
    for i=1:Du
        if Du==2
            subplot(2,1,i);
        else
            subplot(2,2,i);
        end
        hold on;
        stairs((1:length(a(:,i)))*plant.dt,a(:,i),colorVec{J+j});
        legend(iterVec{1:J+j});
        xlabel('Time [s]');     ylabel(actionTitles{i});
        axis tight
    end
    
    run plotCost.m
    
    if plotting.verbosity > 1
        run plotModel.m
        
        % Actions during latest rollout:
        if ~ishandle(11)
            figure(11);
        else
            set(0,'CurrentFigure',11)
        end
        clf(11);
        for i=1:Du       % plot the rollouts on top of predicted error bars
            if Du==2
                subplot(2,1,i); hold on;
            elseif Du==4
                subplot(2,2,i); hold on;
            end
            errorbar( (1:length(Mcon{j}(i,:)))*plant.dt, Mcon{j}(i,:), ...
                2*sqrt(squeeze(Scon{j}(i,i,:))),'r');  hold on;
            stairs((1:length(a(:,i)))*plant.dt,a(:,i),'b');
            if i==Du; xlabel('\fontsize{14}Time [s]'); end;
            ylabel(strcat('\fontsize{14}',actionTitles{i}),'interpreter','Tex');
            if i==1; title('Predicted and Recorded Policy Outputs'); end
            axis tight;     grid on;
        end
        drawnow;
        
        if plotting.verbosity > 2
            % robot animation
            q_sim = latent{j+J}(:,1:robot.n);
            if ~ishandle(5)
                figure(5);
            else
                set(0,'CurrentFigure',5);
            end
            clf(5);
            patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
            robot.plot(q_sim);
            
            if ~ishandle(7)         % cost iterations
                figure(7);
            else
                set(0,'CurrentFigure',7);
            end
            hold on; grid on;
            stairs(1:length(realCost{J+j}),realCost{J+j},colorVec{J+j});
            legend(iterVec{1:J+j});
            title('Recorded Cost in Rollout');
            xlabel('Time [s]');     ylabel('Immediate Cost');
            drawnow;
        end
    end
end

%% 3. Save data
filename = [basename num2str(j) '_H' num2str(H)];
save(filename);
