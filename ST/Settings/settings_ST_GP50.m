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
    'x','y','z','\alpha','\beta','\gamma','dx/dt','dy/dt','dz/dt',...
    'd\alpha}','d\beta}','d\gamma}','f_x','f_y','f_z','\tau_x','\tau_y','\tau_z','t','\sigma_f','\sigma_w'};

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
stateLength = 24;
indices = 1:1:stateLength;
n=3;
odei = 1:1:2*n+6;
augi    = [];                                           % augi  indicies for variables augmented to the ode variables
angi    = [];                                           % angi  indicies for variables treated as angles (using sin/cos representation) (subset of indices)
dyno    = [7 8 13 14 19 20];                            % dyno  indicies for the output from the dynamics model and indicies to loss    (subset of indices)
dyni    = [1 2 3 4];                                    % dyni  indicies for inputs to the dynamics model                               (subset of dyno)
difi    = [1 2 3 4];                                    % difi  indicies for training targets that are differences                      (subset of dyno)
poli    = [1 2 3 4];                                    % poli  indicies for variables that serve as inputs to the policy               (subset of dyno)

% force as input just doesnt work
% force as difference, not input, gives good velocities, bad forces
% force as absolute, not as input...    

REF_PRIOR  = 0;                                         % encode prior reference mean?
refi    = [];                                           % indices for which to encode a reference as  prior mean
ref_select = [1 2 7 8 13 14];                           % indices of reference corresponding to dyno    [xe dxe F]

dynoTitles = stateNames(indices(dyno));
actionTitles = {'Kp_x  [N/m]', 'Kp_{y}  [N/m]', 'x_{ref}','y_{ref}'};%, 'Kp_{rot} [Nm/rad]'};

%% 2. Set up the scenario
% General:
T = 10.0;                % [s] Rollout time
N = 25;                            % no. of controller optimizations
Ntest = 3;                         % no. of roll outs to test controller quality
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

% Timing constraints:
dt = 0.005;              % [s] controller sampling time
dt_pilco = 0.1;         % [s] PILCO sampling rate 
t_pilco = (0:dt_pilco:T)';
t = 0:dt:T;
H = ceil(T/dt_pilco);              % no. of timesteps per rollout

% Robot model:
dynPert = 0.15;
run init_3Lbot.m

% Spatial constraints:
mode = 0;    % mode
xhole = [0.7, 0.2, 0];   % center hole location [x, y, phi/z]
xc    = [0.65, 2, 2, 10, 10, 10]';  % [m] environment constraint location

if ~ishandle(5)         % robot animation
    figure(5);
else
    set(0,'CurrentFigure',5);
end
% First rollout:
H0   = transl([0.5 0 0]);      % start pose end-effector (nominal start pose)
H1   = transl([0.7 0 0]);
H2 = transl(0.7, 0.2, 0);
H3 = transl(0.75, 0.2, 0);
[mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf, Rd, Td]...
    = genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T, dt);
deltaXe_des = diff(xe_des(1:length(t),2:end));
aR_init{1} = timeseries(deltaXe_des',t(1:length(deltaXe_des)));
robot.plot(mu0(1:robot.n));

% Second rollout:
H0 = transl([0.3 0.3 0]);       % start pose end-effector
H1 = transl([0.5 0.15 0]);       
H2 = transl(0.6, 0.25, 0);
H3 = transl(0.55, 0.15, 0);
[mu01, ~, xe_des, ~, ~, ~, ~, ~, ~]...
    = genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T, dt);
deltaXe_des = diff(xe_des(1:length(t),2:end));
aR_init{2} = timeseries(deltaXe_des',t(1:length(deltaXe_des)));
robot.plot(mu01(1:robot.n));

% Third rollout:
H0 = transl([0.35 0.1 0]);       % start pose end-effector
H1 = transl([0.5 0.15 0]);       
H2 = transl(0.6, 0.05, 0);
H3 = transl(0.55, 0.15, 0);
[mu02, ~, xe_des, ~, ~, ~, ~, ~, ~]...
    = genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T, dt);
deltaXe_des = diff(xe_des(1:length(t),2:end));
aR_init{3} = timeseries(deltaXe_des',t(1:length(deltaXe_des)));

initialMu0 = [mu0; mu01; mu02];

initPredVar = 0.001^2;                               % initial state variance around mean
startStateInterval = [0 0 0 0 0 0]';             % interval for which start states may vary [x y z]
mu0Sim = mu0(dyno);                                 % initial mean for simulation
S0Sim = diag(ones(1,length(dyno))*initPredVar);     % covariance matrix of initial state distribution during simulation 

if plotting.verbosity > 1
    figure(15);
    subplot(3,1,1)
    plot(xe_des(:,2:end));
    title('Reference Trajectories');
    ylabel('Positions [m]');
    grid on
    axis tight
    subplot(3,1,2)
    plot(dxe_des(:,2:end))
    ylabel('Velocity [m/s]')
    legend('x_e','y_e','z_e','Location','Best');
    axis tight 
    grid on
    subplot(3,1,3)
    plot(ddxe_des(:,2:end))
    ylabel('Acceleration [m/s/s]'); xlabel(strcat('Time steps   (d_t = ',num2str(dt), ')'));
    axis tight
    grid on
end

% Environment:
Kp_env = [4e3, 2e3 1e4, 0, 0, 0];            %[N/m]  stiffness  (x, y, z, rotx, roty, rotz)
Kd_env = [1, 1, 1, 0, 0, 0];              %[Ns/m] damping

% Display Scenario in Console:
fprintf('\nFinal transformation: \nH^0_n(T) = \n\n');
disp(Hf);
fprintf('\nEnvironment stiffness \t= \t %6.2f [N/m]\n', Kp_env(1));
fprintf('Environment damping \t= \t %6.2f [Ns/m]\n',  Kd_env(1));
fprintf('X-location of environment = \t %6.2f [Ns/m]\n',  xc(1));
fprintf('Hole location \t\t = \t')
disp(xhole)


%% 3. Set up the plant structure
outputNoiseSTD = ones(1,length(odei))*0.0017.^2;                          % noise added to odei indicies in simulation
outputNoiseSTD(1,robot.n+1:robot.n*2) = 0.0001.^2;                          % noise added to odei indicies in simulation
outputNoiseSTD(1,end-5:end) = 0.1^2;

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
plant.rollout_model = 'IPILCO_relativeRPY_ST';
plant.indices = indices;
plant.startStateInterval = startStateInterval;

%% 4. Set up the policy structure
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  % linear saturating controller
maxVel = 0.25;         % maximum Cartesian velocity
policy.maxU  = [250/2 250/2, dt*maxVel  dt*maxVel];
policy.minU  = [10    10,   -dt*maxVel -dt*maxVel];
policy.impIdx = [1, 2]; 			% non-negative indices of policy outputs (saturate + translate)
policy.refIdx = [3, 4];             % reference indices (only saturate)
policy.SNR = 100;
Du = length(policy.maxU);
translVec = [ones(size(policy.impIdx)).*2, ones(size(policy.refIdx))];

aK_init = genInitActions(policy, J, 2, actionTitles, t_pilco, H);
% aK_init = genInitActions(policy, J, 3, actionTitles, t_pilco, 6, 0.2);
% aK_init = genInitActions(policy, J, 4, actionTitles, t_pilco, 3, 0.5);

nc = 50;
policy.p.inputs  = gaussian(mu0Sim(poli), diag(ones(1,length(poli))*0.1), nc)';                % policy pseudo inputs   [ N  x  d ]
policy.p.targets = 0.1*randn(nc, length(policy.maxU));                                         % init. policy targets 
policy.p.hyp = ...                                                                             % GP-log hyperparameters [(d+2) x  D ]
    repmat(log([ones(1,length(poli))*1, 1, 1/policy.SNR]'), 1, length(policy.maxU));

%% 5. Set up the cost structure
cost.fcn   = @my_lossAdd2;                       % cost function
cost.gamma = 1;                                 % discount factor  =1 for finite horizon
cost.expl  = -0.25;                             % exploration weight (UCB) smoothes the value function out and simplifies the optimization problem.

% Energy penalty parameters:
cost.ep  = 0.01;                              % weight
idx = policy.impIdx;

nonIdx = find(~ismember(1:Du,idx));
normalizer = (policy.maxU*diag(translVec)).^2;
quadraticWidth  = diag(normalizer);          % normalization matrix for quadratic ep
iT = inv(quadraticWidth);
iT(nonIdx,nonIdx)= 0;    
cost.iT = iT;

cost.sub{1}.fcn     = @lossSat_2dPIH;
cost.sub{1}.losi    = [1];                        % indicies for saturating cost states
cost.sub{1}.target  = (xhole(1) + [0.05])';   % target state
cost.sub{1}.width   = 0.05;
cost.sub{1}.angle   = plant.angi;

cost.sub{2}.fcn     = @lossSat_2dPIH;
cost.sub{2}.losi    = [5 6];                        % indicies for saturating cost states
cost.sub{2}.target  = [0 0];   % target state
cost.sub{2}.width   = 25;
cost.sub{2}.angle   = plant.angi;

%% 6. Set up the GP dynamics model structure
dynmodel.fcn    = @my_gp1d;                    % function for GP predictions
dynmodel.train  = @my_train;                % function to train dynamics model
nii             = 300;                      % no. of inducing inputs
dynmodel.induce = zeros(nii,0,1);% shared/individual inducing inputs per target dim (sparse GP)
noisyInputs     = false;                    % if true -> train/regress w/ assumed input noise hyperparams
inputNoiseSTD   = [ones(1,length(dyno))*0.005^2, ones(1,length(policy.maxU))*1e-10.^2];      % starting estimate for the noisy input GP training
dynmodel.parallel = true;                  % train individual target dimensions in parellel
trainOpt        = [200 300];                % max. number of line searches [full, sparse]
compareToFullModel = true;

%% 7. Parameters for policy optimization
opt.fh = 1;
opt.method = 'BFGS';                    % 'BFGS' (default), 'LBFGS' (x>1000), 'CG'
opt.length = 50;                        % (+): max. number of line searches
opt.MFEPLS = 25;                        % max. number of function evaluations per linesearch
% opt.MSR = 100;                        % max. slope ratio (default=100)

%% Final inits
prevLength = 0;
x = []; y = []; r = [];
fantasy.mean = cell(1,N+J); fantasy.std = cell(1,N+J);
realCost = cell(1,N+J);  latent = cell(1,N+J);  realAcumCost = zeros(1,N+J);
testLat = cell(1,N); testCost = cell(1,N);
M = cell(N,1);  Sigma = cell(N,1);
Mfull = cell(N,1); Sfull = cell(N,1);
insertSuccess = cell(1,N+1);    scoreCard = zeros(1,N+1);
reference = zeros(H+1,length(dyno));