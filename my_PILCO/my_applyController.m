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

[xx, yy, realCost{j+J}, latent{j+J}, rr] = ...
    my_rollout(mu0, policy, HH, plant, robot); %#ok<*IJCL>
x = [x; xx]; y = [y; yy];
robs = [mu0(1,dyno); rr(:,ref_select)];  rrr = rr(:,ref_select);
rrr(:,difi) = rr(:,ref_select(difi)) - robs(1:size(rr,1),difi);
r = [r; rrr];      %#ok<*AGROW> % augment training sets for dynamics model
realAcumCost(j+J) = sum(realCost{j+J});                 % Accumulated cost
realAcumCost2{j+J}(1) = sum(realCost{j+J});                 % Accumulated cost

% determine if rollout aborted, failed or successful:

if size(latent{j+J},1) > H-5
    insertSuccess(j+J) = 1;                         %#ok<*SAGROW> not aborted
    if abs(mean(latent{j+J}(end-5:end,dyno(1)) - (xhole(1)+0.05))) < 0.01;
        insertSuccess(j+J) = 2;                     % successful insertion
    end
end


% 2. Make many rollouts to test the controller quality / robustness?
testLat = cell(1,Ntest);
testCost = cell(1,Ntest);
for i=1:Ntest
    [~,~,testCost{i},testLat{i}] = my_rollout(mu0, policy, HH, plant, robot);
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
        stairs(1:length(testCost{ii}),testCost{ii},'g');
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
    clf;
    a{j+J} = xx(:,end-Du+1:end);
    for k=1:j+J
        for i=1:Du
            subplot(ceil(Du/sqrt(Du)),ceil(sqrt(Du)),i);
            hold on;
            stairs(1:length(a{k}(:,i)),a{k}(:,i),colorVec{k});
            xlabel('Timestep');     ylabel(actionTitles{i});
        end
        legend(iterVec{1:J+j});
    end
    if numel(dynmodel.induce) ~= 0
        for i=1:Du            
            plot(zeros(nii,1),dynmodel.induce(:,end-Du+1:end),'kx');
        end
    end
    drawnow;
    
    if ~ishandle(10)         % cost iterations
        figure(10);
    else
        set(0,'CurrentFigure',10);
    end
    clf(10);
    hold on; grid on;
    failMedSuc{1} = find(insertSuccess(1:j+J)==0);
    failMedSuc{2} = find(insertSuccess(1:j+J)==1);
    failMedSuc{3} = find(insertSuccess(1:j+J)==2);
    color = {'rx','b+','go'};
    
    % plot recorded accumulated cost:
    for k=1:3
        for i=1:length(failMedSuc{k})
            plot(failMedSuc{k},realAcumCost(failMedSuc{k}),color{k},'MarkerSize',10,'LineWidth',1.5);
            hold on
        end
        hb(k) = plot(0,0,color{k}, 'visible', 'off','MarkerSize',10,'LineWidth',1.5);   % dummy plot for legend
    end
    plot(realAcumCost(1:j+J),'r--','LineWidth',0.1);
    
    % plot predicted accumulated costs:
    for i=1:j, errorbar(i+J,sum(fantasy.mean{i}), 2*sum(fantasy.std{i}), 'k*');    end
    hb(4) = plot(0,0,'k*','visible','off','MarkerSize',10,'LineWidth', 1.5);
    legend(hb,'Aborted','Failed','Successful','Predicted');
    title('Accumulated Rollout Cost');   xlabel('Learning iteration');   ylabel('Total Cost');
    ax = gca; ax.XTick = 1:1:J+j;
    
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
            
            
            % Full/Sparse model:
                errorbar(0:length(M{j}(i,:))-1, M{j}(i,:), ...
                    2*sqrt(squeeze(Sigma{j}(i,i,:))),'b');
                           
            % Full model:
            if compareToFullModel && numel(dynmodel.induce)~=0
                errorbar(0:length(Mfull{j}(i,:))-1, Mfull{j}(i,:), ...
                    2*sqrt(squeeze(Sfull{j}(i,i,:))), 'y');
            end
            
            % Test trials:
            for ii=1:Ntest
                stairs(0:size(testLat{ii}(:,dyno(i)),1)-1, testLat{ii}(:,indices(dyno(i))), 'g' );        % recorded latent states in multiple robustness test-rollouts
            end
            
            % Real trial:
            stairs(0:size(latent{j+J}(:,dyno(i)),1)-1, latent{j+J}(:,indices(dyno(i))),'r');      % recorded latent states in apply_controller roll-out
            
            % Inducing inputs Locations:
            if i <= length(dyni) && numel(dynmodel.induce) ~= 0
                plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
            end
            
            title(dynoTitles{i});
            if i==1
                if numel(dynmodel.induce)==0
                    legend('Full Model','Trial','Sparse Inputs','Location','Best');
                else
                    legend('Sparse Model','Full Model','Trial','Sparse Inputs','Location','Best');
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

%% Plot Policy Progression
if j==N
    figure;
    for i=1:N+J
        plot(latent{i}(:,end-2:end));
        title(num2str(i));
        legend('x','y','r');
        pause;
    end
end