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
if isfield(plant,'constraint')
    HH = maxH;
else
    HH = H;
end

if KUKA
    [xx, yy, realCost{j+J}, latent{j+J}, rr] = ...
        my_iiwaRollout(policy, plant, cost, HH, Hdes);  % Experiment
    
else
    [xx, yy, realCost{j+J}, latent{j+J}, rr] = ...
        my_rollout(mu0, policy, HH, plant, robot);      % Simulation
end
rr = rr(:,ref_select);      % filter out the relevant reference dimensions
rrr = rr(2:end,:);
rrr(:,difi) = rr(2:end,difi) - rr(1:end-1,difi); % dimensions that are differences
r = [r; rrr];      %#ok<*AGROW> % augment training sets for dynamics model
x = [x; xx]; y = [y; yy];

realAcumCost(j+J) = sum(realCost{j+J});                 % Accumulated cost
realAcumCost2{j+J}(1) = sum(realCost{j+J});                 % Accumulated cost

lengthDiff = H - size(realCost{j+J},1);
if lengthDiff == 0
    insertSuccess{j+1}(1) = 1;                         %#ok<*SAGROW> not aborted
    if abs(mean(latent{j+J}(end-5:end,dyno(1)) - (xhole(1)+0.05))) < 0.01;
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
        [xx, yy, testCosti{i}, testLati{i}, rr] = ...
            my_iiwaRollout(policy, plant, cost, HH, Hdes);  % Experiment
        
    else
        [xx, yy, realCost{i}, latent{i}, rr] = ...
            my_rollout(mu0, policy, HH, plant, robot);      % Simulation
    end
    lengthDiff = H-size(testCosti{i},1);
    if lengthDiff==0                           % rollout was completed
        insertSuccess{j+1}(i+1) = 1;
        if abs(mean(testLati{i}(end-5:end,dyno(1)) - (xhole(1)+0.05))) < 0.01;
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
elseif length(find(insertSuccess{j+1}==2,2))==Ntest+1
    scoreCard(j+1) = 2;                 % All Success
else
    scoreCard(j+1) = 1;                 % Partial Success
end

%% verbosity
if plotting.verbosity > 0
    if ~ishandle(3)         % Cost plot
        figure(3);
    else
        set(0,'CurrentFigure',3);
    end
    hold on;
    stairs(1:length(realCost{J+j}),realCost{J+j},'r'); hold on;
    for ii=1:Ntest
        stairs(1:length(testCost{j}{ii}),testCost{j}{ii},'g');
        hold on;
    end
    title('Predicted (uncertain) & Rollout (deterministic) Immediate Cost');
    xlabel('Time [s]');     ylabel('Immediate Cost');
    legend('predicted','test','verifications');
    drawnow;
    
    if ~ishandle(6)         % Action plot
        figure(6);
    else
        set(0,'CurrentFigure',6);
    end
    hold on;
    a = xx(:,end-Du+1:end);
    for i=1:Du
        subplot(ceil(Du/sqrt(Du)),ceil(sqrt(Du)),i);
        hold on;
        stairs(1:length(a(:,i)),a(:,i),colorVec{J+j});
        legend(iterVec{1:J+j});
        xlabel('Timestep');     ylabel(actionTitles{i});
    end
    
    run plotCost.m
    
    if plotting.verbosity > 1
        if ~ishandle(4)
            figure(4);
        else
            set(0,'CurrentFigure',4)
        end
        clf(4);
        ldyno = length(dyno);
        for i=1:ldyno       % plot the rollouts on top of predicted error bars
            subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
            
            % Full model:
            if compareToFullModel && ~isempty(Mfull{j})
                errorbar(0:length(Mfull{j}(i,:))-1, Mfull{j}(i,:), ...
                    2*sqrt(squeeze(Sfull{j}(i,i,:))), 'y');
            end
            
            % Full or Sparse model:
            errorbar( 0:length(M{j}(i,:))-1, M{j}(i,:), ...
                2*sqrt(squeeze(Sigma{j}(i,i,:))) );
            
            % Test trial:
            for ii=1:Ntest
                stairs( 0:size(testLat{j}{ii}(:,dyno(i)),1)-1, testLat{j}{ii}(:,indices(dyno(i))), 'g' );        % recorded latent states in multiple robustness test-rollouts
            end
            
            % Model Trial:
            stairs( 0:size(latent{j+J}(:,dyno(i)),1)-1, latent{j+J}(:,indices(dyno(i))),'r');      % recorded latent states in apply_controller roll-out
            
            % Inducing Inputs:
            if i <= length(dyni) && numel(dynmodel.induce) ~= 0
                plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
            end
            
            title(dynoTitles{i});
            if i==1
                if compareToFullModel && numel(dynmodel.induce) ~= 0
                    if Ntest > 0
                        legend('Full Model','Sparse Model','Test','Trial','Inducing Inputs','Location','Best');
                    else
                        legend('Full Model','Sparse Model','Trial','Inducing Inputs','Location','Best');
                    end
                else
                    if Ntest > 0
                        legend('Full Model','Test','Trial','Location','Best');
                    else
                        legend('Full Model','Trial','Location','Best');
                    end
                end
            end
            axis tight
        end
        drawnow;
        
        if plotting.verbosity > 2
            q_sim = latent{j+J}(:,1:robot.n);
            if ~ishandle(5)         % robot animation
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
