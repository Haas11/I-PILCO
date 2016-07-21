
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
settings_S;                  % load scenario-specific settings
ep=num2str(cost.ep);
expl=num2str(cost.expl);
basename = strcat('_',date,'_','conLin-','_ep-0p',ep(strfind(ep,'.')+1:end),...
    '_expl-0p',expl(strfind(expl,'.')+1:end),'-');      % filename used for saving data

% numerically test my_gSat for proper means, variances and gradients
if diffChecks && satCheck,    gSinSatT(@my_gSat); end

printTuningParams

%% 2. Initial rollouts
fprintf('\nPerforming initial rollouts...\n');
initRollout = 1;    %#ok<*NASGU>     Random walk (OU process)
constMean = 1;      % Gaussian inputs w/ constant mean
for jj = 1:J
    [xx, yy, realCost{jj}, latent{jj}, rr] = ...
        my_iiwaRollout(mu0, policy, H, plant, robot, a_init);   
    realAcumCost(jj) = sum(realCost{jj});   
    x = [x; xx]; y = [y; yy]; 
    robs = [mu0(1,dyno); rr(:,ref_select)];  rrr = rr(:,ref_select);
    rrr(:,difi) = rr(:,ref_select(difi)) - robs(1:length(rr),difi);
    r = [r; rrr];      %#ok<*AGROW> % augment training sets for dynamics model
    
    % determine if rollout successful:
    if peg
        if size(xx,1) == H
            insertSuccess(jj) = 1;  %#ok<*SAGROW> not aborted
            if mean(abs(latent{1}(end-5:end,dyno(1))-(xhole(1)+0.05))) < 0.01;
                insertSuccess(jj) = 2;                % succesful insertion
            end
        end
    end
   
    if insertSuccess(jj) ~= 0;
        REF_DIFF = rrr;
        dynmodel.ref = rrr;
    else
        warning('Initial rollout aborted. Complete reference unavailable.');
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
        
        if ~ishandle(10)         % cost iterations
            figure(10);
        else
            set(0,'CurrentFigure',10);
        end
        clf(10);
        hold on; grid on;
        failMedSuc{1} = find(insertSuccess(1:jj)==0);
        failMedSuc{2} = find(insertSuccess(1:jj)==1);
        failMedSuc{3} = find(insertSuccess(1:jj)==2);
        color = {'rx','bo','g+'};
        for k=1:3
            for i=1:length(failMedSuc{k})
                plot(failMedSuc{k},realAcumCost(failMedSuc{k}),color{k},'MarkerSize',10,'LineWidth',1.5);
                hold on
            end
            hb(k) = plot(0,0,color{k}, 'visible', 'off','MarkerSize',10,'LineWidth',1.5);   % dummy plot for legend
        end
        hb(4) = plot(0,0,'k*','visible','off','MarkerSize',10,'LineWidth', 1.5);
        legend(hb,'Aborted','Failed','Successful','Predicted');
        title('Accumulated Rollout Cost');   xlabel('Learning iteration');   ylabel('Total Cost');
        ax = gca; ax.XTick = 1:1:jj;
        
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

%% 3. Controlled learning (N iterations)
fprintf('\nPILCO Learning started\n--------------------------------\n');
jj=J;
initRollout = 0;
constMean = 0;
for j = 4:N
    my_trainDynModel;       % train (GP) dynamics model
    my_learnPolicy;         % update policy
    my_applyController;     % apply controller to system
    
    disp(['\nControlled trial # ' num2str(j)]);
    disp('Insertion successes: '); disp(insertSuccess);
end
%%
figure;
for i=1:N+J
    plot(latent{i}(:,end-2:end));
    title(num2str(i));
    legend('x','y','r');
    pause;
end