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

%  1  theta1         angle first link
%  2  theta2         angle second link
%  3  theta3         angle third link

%  4  dtheta1        angular velocity of 1st pendulum
%  5  dtheta2        angular velocity of 2nd pendulum
%  6  dtheta3        angular velocity of 2nd pendulum

%  7  X              x-position end-effector in world frame
%  8  Y              y-position end-effector in world frame
%  9  Z              z-position end-effector in world frame
% 10  Alpha
% 11  Gamma
% 12  Phi            angle around z

% 13  dX             x-velocity end-effector in world frame
% 14  dY             y-velocity end-effector in world frame
% 15  dZ             z-velocity end-effector in world frame
% 16  dAlpha         Alpha-angular velocity  in world frame
% 17  dGamma         Gamma-angular velocity  in world frame
% 18  dPhi           Phi-angular velocity    in world frame

% 19  Fx             Cartesian x-force in world frame
% 20  Fy             Cartesian y-force in world frame
% 21  Fz             Cartesian z-force in world frame
% 22  Fa             Cartesian alpha torque in world frame
% 23  Fg             Cartesian gamma torque in world frame
% 24  Fp             Cartesian phi   torque in world frame

% 25 t               time vector

% 25  Kp_x           Cartesian stiffness in x
% 26  Kp_y           Cartesian stiffness in y
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
% 
% cost.ep = epp(count);
% cost.expl = expll(count);

%% Inits
warning('on','all');
rng(0,'twister');
format short; format compact;
global diffChecks REF_PRIOR 
global diffTol conCheck gpCheck propCheck valueCheck satCheck lossCheck checkFailed
diffChecks = 0;  diffTol   = 1e-3;  checkFailed = 0;
conCheck   = 0;  gpCheck   = 0;     propCheck   = 0;    % GP check = very heavy!
valueCheck = 1;  satCheck  = 0;     lossCheck   = 0;
stateNames = {'theta1','theta1','theta1','dtheta1','dtheta1','dtheta1',...
    'xe_x','xe_y','xe_z','alpha','gamma','phi','dxe_x','dxe_y','dxe_z',...
    'dalpha','dgamma','dphi','Fx','Fy','Fz','Tx','Ty','Tz','t'};

if diffChecks
    fprintf('==================================\n'); %#ok<UNRCH>
    fprintf('DERIVATIVE CHECKS ENABLED! [DEBUG]\n');
    fprintf('Ctrl: %1i, GP: %1i, Prop: %i, Value: %1i, Sat: %1i, Loss: %1i',...
        conCheck, gpCheck, propCheck, valueCheck, satCheck, lossCheck);
    fprintf('\n==================================\n\n');
end

opt.verbosity = 3;                      % optimization verbosity      [0-3]
plotting.verbosity = 2;                 % plotting verbosity          [0-3]


%% 1. Define state indices
stateLength = 25;
indices = 1:1:stateLength;
n=3;
odei = 1:1:2*n+6;
augi    = [];                                           % augi  indicies for variables augmented to the ode variables
angi    = [];                                           % angi  indicies for variables treated as angles (using sin/cos representation) (subset of indices)
dyno    = [7 8 12 13 14 18 19];                          % dyno  indicies for the output from the dynamics model and indicies to loss    (subset of indices)
dyni    = [1 2 3 4 5 6];                              % dyni  indicies for inputs to the dynamics model                               (subset of dyno)
difi    = [];                              % difi  indicies for training targets that are differences                      (subset of dyno)
poli    = [1 2 3];                                    % poli  indicies for variables that serve as inputs to the policy               (subset of dyno)
refi    = [1 2 3 4 5 6];                                 % indices for which to encode a reference as  prior mean
REF_PRIOR  = 1;
ref_select = [1 2 6 7 8 12 13];                    % indices of reference corresponding to dyno    [xe dxe F]

dynoTitles = stateNames(indices(dyno));
actionTitles = {'Kp_x  [N/m]', 'Kp_{y/z}  [N/m]', 'Kp_{rot}  [Nm/rad]'};%, 'Kp_{rot} [Nm/rad]'};

% LOWER INSERTION STIFFNESS
% INCREASE EXPLORATION PARAMETER
% DECREASE LINESEARCHES 
% OFFSET INITIAL ORIENTATION
% REMOVE PENALTY ON PHI
%% 2. Set up the scenario
dynPert = 0.10;          % [%] Perturbation of dynamics during simulation w.r.t ID
dt_sim = 0.01;
dt = dt_sim;             % [s] controller sampling time
dt_pilco = 0.1;          % [s] PILCO sampling rate 
fprintf('\nInitializing robot model');
run init_3Lbot.m
% global ROBOT;
% ROBOT = perturbedRobotNF;

T = 10.0;                % [s] Rollout time
t_pilco = (0:dt_pilco:T)';
peg = 1;                 % [bool]  peg insertion trajectory
xhole = [0.5, 0.3, 0];   % center hole location [x, y, phi/z]
xc    = [0.45, 2, 2, 10, 10, 10]';  % [m] environment constraint location
x0    = [0.4 0 0];
T0    = transl(x0);      % start pose end-effector
T11   = transl([0.5 0.02 0]);
[mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf, Rd, Td]...
    = genTrajectory(robot, peg, T0, T11, xhole, xc, T, dt);
t = 0:dt:T;
Hd = timeseries(Td(:,:,1:length(t)),t);
T_e_init = cardatol(tr2rpy(Hd.data(:,:,1)),1,2,3);      %xyz

initPredVar = 1e-4;
mu0Sim = mu0(indices(dyno));
S0Sim = diag(ones(1,length(dyno)).*initPredVar);

fprintf('\nFinal transformation: \nH^0_n(T) = \n\n');
disp(Hf);

% Environment:
Kp_env = [5e3, 1e3, 0, 0, 0, 0];              %[N/m]  stiffness  (x, y, z, rotx, roty, rotz)
Kd_env = [200, 50, 0, 0, 0, 0];              %[Ns/m] damping
fprintf('\nEnvironment stiffness \t= \t %6.2f [N/m]\n', Kp_env(1));
fprintf('Environment damping \t= \t %6.2f [Ns/m]\n',  Kd_env(1));
fprintf('X-location of environment = \t %6.2f [Ns/m]\n',  xc(1));
fprintf('Hole location \t\t = \t')
disp(xhole)

% Additional:
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

% reference = [downsample(xe_des(:,2:end), ceil(dt_pilco/dt)), downsample(dxe_des(:,2:end), ceil(dt_pilco/dt)), zeros(H+1,6)];
% reference = reference(:,ref_select);
% ref_target = reference(2:H+1,:);
% ref_target(:,difi) = reference(2:H+1,difi) - reference(1:H,difi);
% REF_DIFF = ref_target;                                      % global variable used in my_propagate(d).m
% ref_target_repeat = [];%repmat(ref_target, J-1, 1);

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
plant.simconstraint = @constraint_check;
plant.rollout_model = 'AltEulerAngImpedanceControl';
plant.indices = indices;

%% 4. Set up the policy structure
% policy.fcn = @(policy,m,s)my_mixedConCat(@my_conlin,@my_mixedGSat,policy,m,s);  % linear saturating controller
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  % linear saturating controller
policy.maxU  = [250 250]./2; policy.minU  = [10 10];
policy.impIdx = [1 2]; policy.refIdx = [];
Du = length(policy.maxU);

% multTargets(1,:) = [20 40 40];
% multTargets(2,:) = [30 60 60];
% multTargets(3,:) = [40 80 80];
% varTargets = [20 40 40];              % initial target variance

multTargets(1,:) = [100 100];
multTargets(2,:) = [80 80];
multTargets(3,:) = [120 120];
varTargets = [200 200];              % initial target variance

targets = multTargets(2,:);                % targets for initializing hyperparameters
seedMatrix = 1:1:J*Du;
seedMatrix = reshape(seedMatrix,J,[]);

% policy.p.w = rand(Du,length(poli));
% policy.p.b = ones(Du,1).*targets';

nc = 5;
scale_factor = (policy.maxU.*(3/2)); scaled_targets = targets./scale_factor;
scaled_var = varTargets./scale_factor;

policy.p.inputs  = gaussian(mu0Sim(poli), diag(ones(1,length(poli))*0.1), nc)';                % policy pseudo inputs   [ N  x  d ]
policy.p.targets = ...
    repmat(scaled_targets, nc, 1) + repmat(scaled_var,nc,1).*randn(nc, length(policy.maxU));   % policy pseudo targets  [ N  x  D ]
policy.p.hyp = ...                                                                             % GP-log hyperparameters [(d+2) x  D ]
    repmat(log([ones(1,length(poli))*1, 1, 0.01]'), 1, length(policy.maxU));

% Ohrnstein-Uhlenbeek Stochastic Process for initial rollouts:
% ou_opts = sdeset('RandSeed',2);
% K0 = cell(1,J);
% th = [0.3; 0.2; 0.2];   % drift rate parameters          
% sig = [20; 12; 10];   
% startPoint = [20 20 20;
%               40 40 40;
%               60 125 125];
%         mu = [150 150 150;
%               20 75 75;
%               40 50 50];
% for i=1:J
%     K0{i}(:,1) = t_pilco;
%     K0{i}(:,2:1+length(policy.impIdx)) = min(policy.maxU(1)*2,max(policy.minU(1),sde_ou(th(i),mu(i,:),sig(i),t_pilco,startPoint(i,:),ou_opts)));
% end 

%% 5. Set up the cost structure
cost.fcn   = @my_lossAdd;                        % cost function
% cost.fcn = @my_lossSat;                     % cost function
cost.gamma = 1;                                  % discount factor  =1 for finite horizon
cost.expl  = -0.2;                               % exploration parameter (UCB) smoothes the value function out and simplifies the optimization problem.
cost.ep    = 0.001;                              % energy penalty
% cost.gamma = 1;                             % discount factor
% cost.width = [0.02 0.05];                           % cost function width
% cost.angle = plant.angi;                    % index of angle (for cost function)
% cost.target  = ([xhole(1:2) 0 0] + [0.05 0 0 0])';           % target state
% cost.losi = [1 2];        % relevant indices

cost.sub{1}.fcn     = @lossSat_2dPIH;
cost.sub{1}.losi    = [1 2];                            % indicies for saturating cost states
cost.sub{1}.target  = ([xhole(1:2)] + [0.05 0])';           % target state
cost.sub{1}.width   = 0.05;
cost.sub{1}.angle   = plant.angi;

cost.sub{2}.fcn     = @lossSat_2dPIH;
cost.sub{2}.losi 	= 7;                            % indicies for force
cost.sub{2}.target  = 0;                            % target state
cost.sub{2}.width   = 10;                           % Weight matrix
cost.sub{2}.angle   = plant.angi;                   % index of angle (for cost function)
%% 6. Set up the GP dynamics model structure
dynmodel.fcn    = @gp1d;                    % function for GP predictions
dynmodel.train  = @my_train;                % function to train dynamics model
nii             = 300;                      % no. of inducing inputs
dynmodel.induce = zeros(nii,0,1);% shared/individual inducing inputs per target dim (sparse GP)
noisyInputs     = false;                    % if true -> train/regress w/ assumed input noise hyperparams
inputNoiseSTD   = [ones(1,length(dyno))*0.01^2, ones(1,length(policy.maxU))*1e-10.^2];      % starting estimate for the noisy input GP training
dynmodel.parallel = false;                  % train individual target dimensions in parellel
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
insertSuccess = zeros(1,N+J);
reference = zeros(H+1,length(dyno));

% %%
% fantasy.mean{N} = []; fantasy.std{N} = [];
% realCost{N}=[];  latent{N}=[];  realAcumCost(N) = 0;
% M{N} = [];  Sigma{N} = [];
% insertSuccess(N) = 0;