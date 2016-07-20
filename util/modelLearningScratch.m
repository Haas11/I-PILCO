%% 1. Initialization
clear all; clc; close all;%#ok<CLALL>

%% Inits
covfunc = {'covSum', {'covSEard', 'covNoise'}};        % specify ARD covariance
Mgpr = zeros(H,length(dyno));               Mgp0 = zeros(H,length(dyno));
Mgp0_prob = zeros(H,length(dyno));          Sgp0_prob = zeros(length(dyno),length(dyno),H);
Mgp0_prob2 = zeros(H,length(dyno));         Sgp0_prob2 = zeros(length(dyno),length(dyno),H);
Sgpr = zeros(length(dyno),length(dyno),H);  Sgp0 = zeros(length(dyno),length(dyno),H);

Mgpr2(1,:) = dynmodel.inputs(1,1:length(dyno));
Mgp02(1,:) = dynmodel.inputs(1,1:length(dyno));
Mgp02_prob(1,:) = dynmodel.inputs(1,1:length(dyno));
Mgp02_prob2(1,:) = dynmodel.inputs(1,1:length(dyno));

%% Rollout
for jj=1:J
    random = 1;     %#ok<NASGU>
    [xx, yy, ~, ~] = my_rollout(mu0, policy, H, plant, robot);
    x = [x; xx]; y = [y; yy];
    random = 0;
end

trajectory = x(:,dyno);

%% Train
j=1;
my_trainDynModel;

%% learn policy
plant.prop = @propagated;
[policy.p, ~] = minimize(policy.p, 'my_value', opt, dynmodel.inputs(1,1:length(dyno))', S0Sim, ...
    dynmodel, policy, plant, cost, H);

%% 1) plot full probabilistic trajectories
j=1;
plant.prop = @my_propagated;    % also outputs control trajectories
[M{j}, Sigma{j}] ...
    = my_pred(policy, plant, dynmodel, dynmodel.inputs(1,1:length(dyno))', S0Sim, H);

figure;                     % state trajectories
ldyno = length(plant.dyno);
for i=1:ldyno
    subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
    errorbar( 0:length(M{j}(i,:))-1, M{j}(i,:), ...
        2*sqrt(squeeze(Sigma{j}(i,i,:))) );
    
    %     stairs( 0:size(latent(:,dyno(i)),1)-1, latent(:,dyno(i)),'g');      % recorded latent states in apply_controller roll-out
%     plot(trajectory(1:H,i),'r+');
%     hold on
%     stairs(trajectory(H+1:2*H,i),'g+');
%     hold on
%     stairs(trajectory(2*H+1:3*H,i),'k+');
    title(dynoTitles{i});
    axis tight
    hold on;
end
drawnow;

figure;                     % control trajectories
nU = length(policy.maxU);
for i=1:nU
    subplot(1,2,i); hold on;
    errorbar( 0:length(Mcon{j}(i,:))-1, Mcon{j}(i,:), ...
        2*sqrt(squeeze(Scon{j}(i,i,:))));
    xlabel('Timestep'); if i==1, ylabel('Stiffness [N/m]'); end
    hold on
    stairs(a(1:H,i),'r');
%     plot(dynmodel.inputs(1:H,5+i),'r+');
%     hold on
%     plot(dynmodel.inputs(101:2*H,5+i),'g+');
%     hold on
%     plot(dynmodel.inputs(201:3*H,5+i),'g+');
    legend('pred','rollout');
    axis tight
    hold on;
    title(actionTitles{i});
end
drawnow;

figure;
subplot(5,1,1:2); 
for i=1:2    
    errorbar( 1:length(M{j}(i,:)), M{j}(i,:), ...
        2*sqrt(squeeze(Sigma{j}(i,i,:))) );
    hold on;
    stairs(xx(:,dyno(i)),'r');
    xaxis([0 100])
end
grid on;
subplot(5,1,3);
errorbar( 1:length(M{j}(end,:)), M{j}(end,:), ...
    2*sqrt(squeeze(Sigma{j}(end,end,:))) );
hold on; grid on;
stairs(xx(:,dyno(end)),'r');
xaxis([0 100])
ylabel('Ext. Force [N]')
for i=1:2
    subplot(5,1,3+i);
    errorbar( 1:length(Mcon{j}(i,:)), Mcon{j}(i,:), ...
        2*sqrt(squeeze(Scon{j}(i,i,:))) );    
    hold on
    stairs(a(1:H,i),'r');
    hold on;
    xaxis([0 100])
    grid on
end
xlabel('Timestep');

%% 2) GP0 predictions with deterministic inputs centered at training locations
for i=1:H
    mgp0 = dynmodel.inputs(i,:)';
    [mgp0_out, sgp0_out, ~] ...
        = gp0(dynmodel, mgp0, zeros(length(mgp0)));     % predictions at deterministic inputs
    Mgp0(i,:)   = mgp0_out';
    Sgp0(:,:,i) = sgp0_out;
    Mgp02(i+1,:)= mgp0_out' + Mgp02(i,:);
end

figure;
ldyno = length(plant.dyno);
for i=1:ldyno
    subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i);
    
    errorbar( 1:length(Mgp02(:,i))-1, Mgp02(1:100,i)', 2*sqrt(squeeze(Sgp0(i,i,:))) );
    hold on;
    plot(trajectory(1:100,i),'r+');
    
    title(dynoTitles{i});
    legend('gp0','rollout');
    axis tight
    hold on;
end
drawnow;

%% 3) plot output of GP0 uncertain state inputs at training locations

% 3) GP0 predictions with uncertain inputs centered at training locations  (constant action variance)
sgp0_in = diag(ones(length(dyni)+length(policy.maxU),1)*0);
for i=1:H
    mgp0_in = dynmodel.inputs(i,:)';    % all test inputs at training locations;
    [mgp0_out, sgp0_out, ~] ...
        = gp0(dynmodel, mgp0_in, sgp0_in);     % predictions at deterministic inputs
    Mgp0_prob(i,:)  = mgp0_out';
    Mgp02_prob(i+1,:)= mgp0_out' + Mgp02_prob(i,:);
    
    Sgp0_prob(:,:,i) = sgp0_out;
    sgp0_in = diag([diag(sgp0_out(1:length(dyni),1:length(dyni)))',1,1]);
end


figure;
ldyno = length(plant.dyno);
for i=1:ldyno
    subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i);
    
    errorbar( 1:length(Mgp02_prob(:,i))-1, Mgp02_prob(1:100,i)', 2*sqrt(squeeze(Sgp0_prob(i,i,:))) );
    hold on;
    plot(trajectory(1:100,i),'r+');
    hold on
%     plot(traj(101:200,i),'g+');
%     hold on
%     plot(traj(201:300,i),'k+');
%     
    title(dynoTitles{i});
    legend('gp0','rollout');
    axis tight
    hold on;
end
drawnow;

%% 4) plot output of GP0 uncertain propagated states & centered deterministic actions

% 4) GP0 predictions with uncertain propagated states & deterministic
% actions at training locations
sgp0_in = diag(ones(length(dyno)+length(policy.maxU),1)*0);
mgp0_in = dynmodel.inputs(1,:)';
for i=1:H
    [mgp0_out, sgp0_out, ~] ...
        = gp0(dynmodel, mgp0_in, sgp0_in);     % prediction at uncertain state & deterministic action at training location
    Mgp0_prob2(i,:)  = mgp0_out';
    
    Mgp02_prob2(i+1,:)= mgp0_out' + Mgp02_prob2(i,:);
    mgp0_in = [Mgp02_prob2(i+1,:)'; dynmodel.inputs(i+1,end-1:end)'];    % propagated state input;
    
    Sgp0_prob2(:,:,i) = sgp0_out;
    sgp0_in(1:length(dyno),1:length(dyno)) = sgp0_out;
end
figure;
ldyno = length(plant.dyno);
for i=1:ldyno
    subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i);
    
    errorbar( 1:length(Mgp02_prob2(:,i))-1, Mgp02_prob2(1:100,i)', 2*sqrt(squeeze(Sgp0_prob2(i,i,:))) );
    hold on;
    plot(trajectory(1:100,i),'r+');
    hold on
    plot(trajectory(101:200,i),'k*');
    hold on
    plot(trajectory(201:300,i),'go');
    
    title(dynoTitles{i});
    legend('gp0','rollout 1','rollout 2', 'rollout 3');
    axis tight
    hold on;
end
drawnow;

%% 5) Propagate with informative prior mean

for t = 1:H
    mgp0 = dynmodel.inputs(t,1:5)';
    [m, ~] = plant.prop(mgp0, zeros(length(mgp0)), plant, dynmodel, policy, t);
        
    Mgp02(t+1,:) = m';
end

figure;
ldyno = length(plant.dyno);
for i=1:ldyno
    subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i);
    
    errorbar( 1:length(Mgp02(:,i))-1, Mgp02(1:100,i)', 2*sqrt(squeeze(Sgp0(i,i,:))) );
    hold on;
    plot(trajectory(1:100,i),'r+');
    
    title(dynoTitles{i});
    legend('gp0','rollout');
    axis tight
    hold on;
end
drawnow;


