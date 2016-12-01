%% settings_pendulum.m
% *Summary:* Script set up the pendulum scenario
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-24
%
%% High-Level Steps
% # Define state and important indices
% # Set up scenario
% # Set up the plant structure
% # Set up the policy structure
% # Set up the cost structure
% # Set up the GP dynamics model structure
% # Parameters for policy optimization
% # Plotting verbosity
% # Some array initializations

%% Code

warning('off','all'); format short; format compact; 

% % include some paths
% try
%   rd = '../../';
%   addpath([rd 'base'],[rd 'util'],[rd 'gp'],[rd 'control'],[rd 'loss']);
% catch
% end

rand('seed',5); randn('seed',13); 
stateNames = {'dtheta','theta','sin(theta)','cos(theta)','\sigma_f','\sigma_w'};

global diffChecks REF_PRIOR
REF_PRIOR = 0;
global diffTol conCheck gpCheck propCheck valueCheck satCheck lossCheck checkFailed
diffChecks = 0;  diffTol   = 1e-3;  checkFailed = 0;
conCheck   = 0;  gpCheck   = 0;     propCheck   = 0;    % GP check = very heavy!
valueCheck = 1;  satCheck  = 0;     lossCheck   = 0;
if diffChecks
    fprintf('==================================\n'); %#ok<UNRCH>
    fprintf('DERIVATIVE CHECKS ENABLED! [DEBUG]\n');
    fprintf('Ctrl: %1i, GP: %1i, Prop: %i, Value: %1i, Sat: %1i, Loss: %1i',...
        conCheck, gpCheck, propCheck, valueCheck, satCheck, lossCheck);
    fprintf('\n==================================\n\n');
end
compareToFullModel = true;
Ntest = 1;
stateLength = 2;
%% 1. Define state and important indices

% 1a. Full state representation (including all augmentations)
%  1  dtheta        angular velocity of inner pendulum
%  2  theta2        angle inner pendulum
%  3  sin(theta)    complex representation ...
%  4  cos(theta)    ... of angle of pendulum
%  5  u             torque applied to the pendulum

% 1b. Important indices
% odei  indicies for the ode solver
% augi  indicies for variables augmented to the ode variables
% dyno  indicies for the output from the dynamics model and indicies to loss
% angi  indicies for variables treated as angles (using sin/cos representation)
% dyni  indicies for inputs to the dynamics model
% poli  indicies for variables that serve as inputs to the policy
% difi  indicies for training targets that are differences (rather than values)

odei = [1 2];               % varibles for the ODE solver
augi = [];                  % variables to be augmented
dyno = [1 2];               % variables to be predicted (and known to loss)
angi = [];                 % angle variables
dyni = [1 2];             % variables that serve as inputs to the dynamics GP
poli = [1 2];             % variables that serve as inputs to the policy [1 3 4]
difi = [1 2];               % variables that are learned via differences

dynoTitles = stateNames(dyno);
actionTitles = {'Kp_x  [N/m]', 'Kp_y  [N/m]'};%, 'Kp_{rot} [Nm/rad]'};
hyperTitles = [dynoTitles, actionTitles, {'\sigma_f','\sigma_w'}];

%% 2. Set up the scenario
dt = 0.1;                      % [s] sampling time
T = 4;                         % [s] prediction time
H = ceil(T/dt);                % prediction steps (optimization horizon)
mu0 = [0 0]';                  % initial state mean
S0 = 0.01*eye(2);              % initial state variance
N = 10;                        % number of policy optimizations
J = 1;                         % no. of inital training rollouts (of length H)
K = 1;                         % number of initial states for which we optimize
nc = 20;                       % size of controller training set

%% 3. Set up the plant structure
plant.dynamics = @dynamics_pendulum;    % dynamics ODE function
plant.noise = diag([0.1^2 0.01^2]);     % measurement noise
plant.dt = dt;
plant.ctrl = @zoh;                  % controler is zero-order-hold
plant.odei = odei;                  % indices of the varibles for the ODE solver
plant.augi = augi;                  % indices of augmented variables
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.refi = [];
plant.prop = @my_propagated;        

%% 4. Set up the policy structure
policy.fcn = @(policy,m,s)conCat(@congp,@gSat,policy,m,s);% controller 
                                                          % representation
policy.maxU = 2.5;                                        % max. amplitude of 
Du = length(policy.maxU);                                 % torque
policy.impIdx = 1; policy.refIdx = [];
translVec = [ones(size(policy.impIdx)), ones(size(policy.refIdx))];

[mm, ss, cc] = gTrig(mu0, S0, plant.angi);                  % represent angles 
mm = [mu0; mm]; cc = S0*cc; ss = [S0 cc; cc' ss];         % in complex plane      
policy.p.inputs = gaussian(mm(poli), ss(poli,poli), nc)'; % init. location of 
                                                          % basis functions
policy.p.targets = 0.1*randn(nc, length(policy.maxU));    % init. policy targets 
                                                          % (close to zero)
policy.p.hyp = log([1 0.7 1 0.01]');                  % initialize 
                                                          % hyper-parameters
%% 5. Set up the cost structure
cost.fcn   = @my_lossAdd2;                      % cost function
cost.gamma = 1;                                 % discount factor  =1 for finite horizon
cost.expl  = 0;                                 % exploration weight (UCB) smoothes the value function out and simplifies the optimization problem.

% Energy penalty parameters:
cost.ep  = 0.0;                                % weight
idx = policy.impIdx;
nonIdx = find(~ismember(1:Du,idx));
normalizer = (policy.maxU*diag(translVec)).^2;
quadraticWidth  = diag(normalizer);             % normalization matrix for quadratic ep
iT = inv(quadraticWidth);
iT(nonIdx,nonIdx)= 0;    
cost.iT = iT;

cost.sub{1}.fcn     = @bench_loss_pendulum;
cost.sub{1}.p       = 1;
cost.sub{1}.losi    = [1 2];                    % indicies for saturating cost states
cost.sub{1}.target  = [0 pi]';                  % target state
cost.sub{1}.width   = 0.5;
cost.sub{1}.angle   = plant.angi;

%% 6. Set up the GP dynamics model structure
dynmodel.fcn = @gp1d;                % function for GP predictions
dynmodel.train = @train;             % function to train dynamics model
dynmodel.induce = zeros(300,0,1);    % shared inducing inputs (sparse GP)
trainOpt = [300 500];                % defines the max. number of line searches
                                     % when training the GP dynamics models
                                     % trainOpt(1): full GP,
                                     % trainOpt(2): sparse GP (FITC)
                                     
%% 7. Parameters for policy optimization
opt.length = 50;                         % max. number of line searches
opt.MFEPLS = 30;                         % max. number of function evaluations
                                         % per line search
opt.verbosity = 3;                       % verbosity: specifies how much 
                                         % information is displayed during
                                         % policy learning. Options: 0-3
opt.method = 'BFGS';                     % optimization algorithm. Options:
                                         % 'BFGS' (default), 'LBFGS', 'CG'


%% 8. Plotting verbosity
plotting.verbosity = 2;            % 0: no plots
                                   % 1: some plots
                                   % 2: all plots

%% 9. Some initializations
prevLength = 0;
x = []; y = []; r = [];
fantasy.mean = cell(1,N+J); fantasy.std = cell(1,N+J);
realCost = cell(1,N+J);  latent = cell(1,N+J);  realAcumCost = zeros(1,N+J);
testLat = cell(N,1);    testCost = cell(N,1);
M = cell(N,1);  Sigma = cell(N,1); Mcon = cell(N,1); Scon = cell(N,1);
Mfull = cell(N,1); Sfull = cell(N,1);
insertSuccess = cell(1,N+1);    scoreCard = zeros(1,N+1);
reference = zeros(H+1,length(dyno));

colorVec = {'r','b','g','k','m','c','y','r-.','b-.','g-.','k-.','m-.','c-.','y-.',...
    'r:','b:','g:','k:','m:','c:','y:','r--','b--','g--','k--','m--','c--','y--'};
colorVec = [colorVec,colorVec];
iterVec = cell(1,J+N);
iterVec{1} = 'iter. 1';
for i=2:J+N
    iterVec{i} = num2str(i);
end