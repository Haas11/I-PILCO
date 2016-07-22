%% settings_Link3_learnImp.m
% *Summary:* Script to set up impedance parameter learning for a 3-link
% robot.
%
%

%% Indices & State definitions:
% command window:   Printed function values                     [opt  >= ?]
% figure(1):        Optimization progress                       [opt  >= 3]
% figure(2):        accumulated cost per iteration              [plot >= 1]
% figure(3):        Predicted & recorded immediate cost         [plot >= 1]
% figure(4):        Predicted & recorded state trajectories     [plot >= 1]
% figure(5):        Robot trajectory animation                  [plot >= 3]
% figure(6):        Rollout actions actions per iteration       [plot >= 1]
% figure(7):        Immediate rollout cost per iteration        [plot >= 1]
% figure(8):        Policy initialization result                [plot >= 2]
% figure(9):        Interaction forces                          [plot >= 2]
% figure(10):       Accumulated rollout costs                   [plot >= 1]

% Full state representation
% - angles              [jointspace]
% - angular velocities  [jointspace]
% - translation         [worldframe]
% - rotation            [worldframe, euler]
% - interaction forces  [worldframe]
% - complex angles      [jointspace]
% - impedance parameters (actions)

%  1  X              x-position end-effector in world frame
%  2  Y              y-position end-effector in world frame
%  3  Z              z-position end-effector in world frame
%  4  Alpha
%  5  Gamma
%  6  Phi            angle around z

% 7  dX             x-velocity end-effector in world frame
% 8  dY             y-velocity end-effector in world frame
% 9  dZ             z-velocity end-effector in world frame
% 10  dAlpha         Alpha-angular velocity  in world frame
% 11  dGamma         Gamma-angular velocity  in world frame
% 12  dPhi           Phi-angular velocity    in world frame

% 13  Fx             Cartesian x-force in world frame
% 14  Fy             Cartesian y-force in world frame
% 15  Fz             Cartesian z-force in world frame
% 16  Fa             Cartesian alpha torque in world frame
% 17  Fg             Cartesian gamma torque in world frame
% 18  Fp             Cartesian phi   torque in world frame

% 19  Kp_x           Cartesian stiffness in x
% 20  Kp_y           Cartesian stiffness in y

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

%% Inits
warning('on','all');
rng(0,'twister');
format short; format compact;
stateNames = {'xe_x','xe_y','xe_z','alpha','gamma','phi','dxe_x','dxe_y','dxe_z',...
    'dalpha','dgamma','dphi','Fx','Fy','Fz','Tx','Ty','Tz'};
global REF_PRIOR 

opt.verbosity = 3;                      % optimization verbosity      [0-3]
plotting.verbosity = 2;                 % plotting verbosity          [0-3]
%% 1. Define state indices
KUKA = true;

stateLength = 24;
indices = 1:1:stateLength;
n=7;
odei = 1:1:2*n+6;
augi    = [];                                           % augi  indicies for variables augmented to the ode variables
angi    = [];                                           % angi  indicies for variables treated as angles (using sin/cos representation) (subset of indices)
dyno    = [3];                          % dyno  indicies for the output from the dynamics model and indicies to loss    (subset of indices)
dyni    = [1];                              % dyni  indicies for inputs to the dynamics model                               (subset of dyno)
difi    = [1];                              % difi  indicies for training targets that are differences                      (subset of dyno)
poli    = [1];                                    % poli  indicies for variables that serve as inputs to the policy               (subset of dyno)
refi    = [1];                                 % indices for which to encode a reference as  prior mean
REF_PRIOR  = 1;
ref_select = [3];                    % indices of reference corresponding to dyno    [xe dxe F]

dynoTitles = stateNames(indices(dyno));
actionTitles = {'Kp_x  [N/m]', 'Kp_{y/z}  [N/m]', 'Kp_{rot}  [Nm/rad]'};%, 'Kp_{rot} [Nm/rad]'};

%% 2. Set up the scenario
dynPert = 0.10;          % [%] Perturbation of dynamics during simulation 
dt_sim = 0.01;
dt = dt_sim;             % [s] controller sampling time
dt_pilco = 0.1;          % [s] PILCO sampling rate 
fprintf('\nInitializing robot model');
run init_3Lbot.m

T = 10.0;                % [s] Rollout time
t_pilco = (0:dt_pilco:T)';

genKUKATrajectory;

mu0Sim = 

H = ceil(T/dt_pilco);              % no. of timesteps per rollout
N = 30;                            % no. of controller optimizations
Ntest = 1;                         % no. of roll outs to test controller quality
J = 1;                             % no. of initial training rollouts
K = 1;                             % no. of initial states for which we optimize
colorVec = {'r','b','g','k','m','c','y','r-.','b-.','g-.','k-.','m-.','c-.','y-.',...
    'r:','b:','g:','k:','m:','c:','y:','r--','b--','g--','k--','m--','c--','y--'};
colorVec = [colorVec,colorVec];
iterVec = cell(1,J+N);
iterVec{1} = 'iter. 1';
for i=2:J+N
    iterVec{i} = num2str(i);
end

%% 3. Set up the plant structure
outputNoiseSTD = ones(1,length(odei))*0.005.^2;                          % noise added to odei indicies in simulation
outputNoiseSTD(1,end-5:end) = 0.05^2;
initRollOutNoise = 1e-4;

plant.noise = diag(outputNoiseSTD);
plant.dt = dt_pilco;
plant.ctrl = @zoh;                  % controler is zero-order-hold
plant.odei = odei;                  % indices to the varibles for the ode solver
plant.augi = augi;                  % indices of augmented variables
plant.angi = angi;
plant.poli = poli;
plant.dyno = dyno;
plant.dyni = dyni;
plant.difi = difi;
plant.refi = refi;
plant.prop = @my_propagated;   % handle to function that propagates state over time
plant.indices = indices;

%% 4. Set up the policy structure
% GP Controller:
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  % linear saturating controller
nc = 10;
policy.maxU  = [150 150]./2; policy.minU  = [10 10];
policy.impIdx = [1 2]; policy.refIdx = [];
Du = length(policy.maxU);

initMean = [1000 1000];
initVar = [100 100];              % initial target variance

targets = initMean;                % targets for initializing hyperparameters
seedMatrix = 1:1:J*Du;
seedMatrix = reshape(seedMatrix,J,[]);

a_init = gaussian(initMean,diag(initVar),H)';

scale_factor = (policy.maxU.*(3/2)); scaled_targets = targets./scale_factor;
scaled_var = initVar./scale_factor;

policy.p.inputs  = gaussian(mu0Sim(poli), diag(ones(1,length(poli))*0.1), nc)';                % policy pseudo inputs   [ N  x  d ]
policy.p.targets = ...
    repmat(scaled_targets, nc, 1) + repmat(scaled_var,nc,1).*randn(nc, length(policy.maxU));   % policy pseudo targets  [ N  x  D ]
policy.p.hyp = ...                                                                             % GP-log hyperparameters [(d+2) x  D ]
    repmat(log([ones(1,length(poli))*1, 1, 0.01]'), 1, length(policy.maxU));

%% 5. Set up the cost structure

% cost.fcn = @my_lossSat;                     % cost function
% cost.gamma = 1;                             % discount factor
% cost.width = [0.02 0.05];                           % cost function width
% cost.angle = plant.angi;                    % index of angle (for cost function)
% cost.target  = ([xhole(1:2) 0 0] + [0.05 0 0 0])';           % target state
% cost.losi = [1 2];        % relevant indices

% Compounded loss function:
cost.fcn   = @my_lossAdd;                     % cost function
cost.gamma = 1;                               % discount factor  =1 for finite horizon
cost.expl  = -0.25;                           % exploration parameter (UCB) smoothes the value function out and simplifies the optimization problem.
cost.ep    = 0.001;                           % energy penalty

cost.sub{1}.fcn     = @lossSat_2dPIH;
cost.sub{1}.losi    = [1 2];                            % indicies for saturating cost states
cost.sub{1}.target  = ([xhole(1:2)] + [0.05 0])';           % target state
cost.sub{1}.width   = 0.05;
cost.sub{1}.angle   = plant.angi;

cost.sub{2}.fcn     = @lossSat_2dPIH;
cost.sub{2}.losi 	= 5;                            % indicies for force
cost.sub{2}.target  = 0;                            % target state
cost.sub{2}.width   = 10;                           % Weight matrix
cost.sub{2}.angle   = plant.angi;                   % index of angle (for cost function)

%% 6. Set up the GP dynamics model structure
dynmodel.fcn    = @my_gp1d;                    % function for GP predictions
dynmodel.train  = @my_train;                % function to train dynamics model
nii             = 300;                      % no. of inducing inputs
dynmodel.induce = zeros(nii,0,1);           % shared/individual inducing inputs per target dim (sparse GP)
noisyInputs     = false;                    % if true -> train/regress w/ assumed input noise hyperparams
inputNoiseSTD   = [ones(1,length(dyno))*0.01^2, ones(1,length(policy.maxU))*1e-10.^2];      % starting estimate for the noisy input GP training
dynmodel.parallel = false;                  % train individual target dimensions in parellel
dynmodel.full   = true;
compareToFullModel = true;                  % Computes the state trajectory of the full model for comparison to sparse approximation
trainOpt        = [200 300];                % max. number of line searches [full, sparse]

%% 7. Parameters for policy optimization
opt.fh = 1;
opt.method = 'BFGS';                    % 'BFGS' (default), 'LBFGS' (x>1000), 'CG'
opt.length = 75;                        % (+): max. number of line searches
opt.MFEPLS = 25;                        % max. number of function evaluations per linesearch
% opt.MSR = 100;                        % max. slope ratio (default=100)

%% Final inits
prevLength = 0;
x = []; y = []; r = [];
fantasy.mean = cell(1,N+J); fantasy.std = cell(1,N+J);
realCost = cell(1,N+J);  latent = cell(1,N+J);  realAcumCost = zeros(1,N+J);
M = cell(N,1);  Sigma = cell(N,1);
Mfull = cell(N,1); Sfull = cell(N,1);

insertSuccess = zeros(1,N+J); reference = zeros(H+1,length(dyno));