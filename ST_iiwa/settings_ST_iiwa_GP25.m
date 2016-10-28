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

%   State Vector 'x':
%  1  X              x-position end-effector in world frame
%  2  Y              y-position end-effector in world frame
%  3  Z              z-position end-effector in world frame
%  4  w              quaternion rotation
%  5  x              quaternion vector x coordinate
%  6  y              quaternion vector y coordinate
%  7  z              quaternion vector z coordinate

% 8  dX              x-velocity end-effector in world frame
% 9  dY              y-velocity end-effector in world frame
% 10  dZ             z-velocity end-effector in world frame
% 11  dw        
% 12  dx         
% 13  dy
% 14  dz            
 
% 15  Fx             Cartesian x-force in world frame
% 16  Fy             Cartesian y-force in world frame
% 17  Fz             Cartesian z-force in world frame
% 18  Fa             Cartesian alpha torque in world frame
% 19  Fg             Cartesian gamma torque in world frame
% 20  Fp             Cartesian phi   torque in world frame

% 21  t              time vector

% 21  Kp_x           Cartesian stiffness in x
% 22  Kp_y           Cartesian stiffness in y
% 23  Kp_z

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
global diffChecks REF_PRIOR
global diffTol conCheck gpCheck propCheck valueCheck satCheck lossCheck checkFailed
diffChecks = 0;  diffTol   = 1e-3;  checkFailed = 0;
conCheck   = 0;  gpCheck   = 0;     propCheck   = 0;    % GP check = very heavy!
valueCheck = 1;  satCheck  = 0;     lossCheck   = 0;
stateNames = {'x','y','z','x_q','y_q','z_q','w_q',...
    'dx/dt','dy/dt','dz/dt','dx_q/dt','dy_q/dt','dz_q/dt','dw_q/dt','f_x','f_y','f_z','\tau_x','\tau_y','\tau_z','t'};

opt.verbosity = 2;                      % optimization verbosity      [0-3]
plotting.verbosity = 2;                 % plotting verbosity          [0-3]

%% 1. Define state indices
stateLength = 21;
indices = 1:1:stateLength;
n=7;
odei = 1:1:2*n+6;
augi    = [];                                           % augi  indicies for variables augmented to the ode variables
angi    = [];                                           % angi  indicies for variables treated as angles (using sin/cos representation) (subset of indices)
dyno    = [1 2 8 9 15 16];                           % dyno  indicies for the output from the dynamics model and indicies to loss    (subset of indices)
dyni    = [1 2];                                    % dyni  indicies for inputs to the dynamics model                               (subset of dynot)
difi    = [1 2 3 4];                                    % difi  indicies for training targets that are differences                      (subset of dynot)
poli    = [1 2 3 4];                                    % poli  indicies for variables that serve as inputs to the policy               (subset of dynot)

KUKA = true;        % boolean to indicate real world experiments
REF_PRIOR  = 0;
refi = [1 2];                                 % indices for which to encode a reference as  prior mean
ref_select = dyno;                          % indices of reference corresponding to dyno    [xe dxe F]

dynoTitles = stateNames(dyno);
actionTitles = {'Kp_{x}','Kp_{y}','Ref_{x}','Ref_{y}'};

%% 2. Set up the scenario
N = 25;                            % no. of controller optimizations
Ntest = 1;                         % no. of roll outs to test controller quality
J = 2;                             % no. of initial training rollouts
K = 1;                             % no. of initial states for which we optimize

% Timing:
T = 10.0;                   % [s] Duration of a single trial
dt_pilco = 0.15;            % [s] PILCO sampling rate 
dt = dt_pilco;              % [s] controller sampling time during rollouts
t_pilco = (0:dt_pilco:T)';
t = 0:dt:T;
H = ceil(T/dt_pilco);              % no. of timesteps per rollout

% REMAINS FROM SIMULATION (DEFINED TO AVOID ERRORS)
xc    = [0, 0, 0.0925, 10, 10, 10]';  % [m] environment constraint location 
xhole = [0.75, 0.125, 0.095];      % center hole location [x, y, phi/z]      (used for judging success with KUKA)
robot = 7;
initRot = t2r(quat2tform([0 1 0 0]));

mode=3;                   % 0= free motion sim, 1=peg sim, 2=KUKA experiment with 2 way-points, 3= KUKA w/ 3 way-points
H0 = rt2tr(initRot, [0.45 0 0.08]');
H1 = rt2tr(initRot, [0.7 0.2 0.08]');         % 0.3 meters in 20 timesteps at 0.15 s = 3 seconds --> max speed = 0.1 m/s
H2 = rt2tr(initRot, [0.7 0.05 0.08]'); 
H3 = rt2tr(initRot, [0.775 0.16 0.08]'); 
[mu0, S0, xe_des, dxe_des, ddxe_des, ~, ~, Hd]=...
    genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T,  dt);
Hdes = Hd;
aR_init{1} = dxe_des(1:H,2:3)*dt;

mode=3;
H1 = rt2tr(initRot, [0.65 0 0.08]');         % 0.3 meters in 20 timesteps at 0.15 s = 3 seconds --> max speed = 0.1 m/s
H2 = rt2tr(initRot, [0.45 0.2 0.08]'); 
H3 = rt2tr(initRot, [0.65 0.05 0.08]'); 
[~, ~, ~, dxe_des, ~, Hf2, Rd2, Hd2]=...
    genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T,  dt);
aR_init{2} = dxe_des(1:H,2:3)*dt;

initPredVar = 1e-4;
startStateInterval = [0 0 0 0 0 0 0]';                                        % interval for which start states may vary [x y z]
mu0Sim = mu0(dyno);
S0Sim = S0(dyno,dyno);

if plotting.verbosity > 1
    figure(15);
    subplot(3,1,1)
    plot(t,xe_des(:,2:4),'Linewidth',1.5);
    [hleg1, hobj1] =legend('x','y','z','Location','NorthEast');
    textobj = findobj(hobj1, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', 16);
    title('\fontsize{16}Motion Plans');  ylabel('\fontsize{14}Positions [m]');
    grid on; axis tight
    subplot(3,1,2)
    plot(t,dxe_des(:,2:4),'Linewidth',1.5)
    ylabel('\fontsize{14}Velocity [m/s]')
    grid on; axis tight;
    subplot(3,1,3)
    plot(t,ddxe_des(:,2:4),'Linewidth',1.5)
    ylabel('\fontsize{14}Acceleration [m/s^2]'); xlabel(strcat('\fontsize{14}Time (d_t=',num2str(dt), ') [s]'),'Interpreter','Tex');
    grid on; axis tight;
end


%% 3. Set up the plant structure
plant.dt = dt_pilco;
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
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  
nc = 25;
maxVel = 0.25;         % maximum Cartesian velocity
policy.maxU  = [500/2 500/2, dt*maxVel  dt*maxVel];
policy.minU  = [25    25,   -dt*maxVel -dt*maxVel];
policy.impIdx = [1, 2]; 			% non-negative indices of policy outputs (saturate + translate)
policy.refIdx = [3, 4]; % reference indices (only saturate)
policy.SNR = 100;
Du = length(policy.maxU);
translVec = [ones(size(policy.impIdx)).*2, ones(size(policy.refIdx))];

a = repmat(mu0(dyno(poli)),nc,1); 
b = repmat([xhole(poli(1:2)), 0 0],nc,1);

policy.p.inputs = a + (b-a).*rand(nc,length(poli));
policy.p.targets = randn(nc, length(policy.maxU));                                         % init. policy targets 
policy.p.hyp = ...                                                                         % GP-log hyperparameters [(d+2) x  nU ]
    repmat(log([ones(1,length(poli))*1, 1, 1/policy.SNR]'), 1, length(policy.maxU));

% a_init = genInitActions(policy, J, 2, actionTitles, t_pilco, 6);     % 1=gaussian, 2=uniform, 3=orhnstein-uhlenbeck
aK_init = genInitActions(policy, J, 3, actionTitles, t_pilco, 5);     % 1=gaussian, 2=uniform, 3=orhnstein-uhlenbeck

for i=1:J
    a_init{i} = [aK_init{i}(1:H,:), aR_init{i}(1:H,:)];
end
    
%% 5. Set up the cost structure
cost.fcn   = @my_lossAdd2;                     % cost function
cost.gamma = 1;                                % discount factor  =1 for finite horizon
cost.expl  = 0;                            % exploration parameter (UCB) smoothes the value function out and simplifies the optimization problem.
cost.ep  = 0.001;                              % weight

idx = policy.impIdx;
nonIdx = find(~ismember(1:Du,idx));
normalizer = (policy.maxU*diag(translVec)).^2;
quadraticWidth  = diag(normalizer);          % normalization matrix for quadratic ep
iT = inv(quadraticWidth);
iT(nonIdx,nonIdx)= 0;    
cost.iT = iT;

cost.sub{1}.fcn     = @my_lossSat;
cost.sub{1}.losi    = [1];                            % indicies for saturating cost states
cost.sub{1}.target  = [0.8];                     % target state  xe=[0.75 0.125 0.095]
cost.sub{1}.width   = 0.1;
cost.sub{1}.angle   = plant.angi;

cost.sub{2}.fcn     = @my_lossSat;
cost.sub{2}.losi    = [5 6];                            % indicies for saturating cost states
cost.sub{2}.target  = [0 0];                            % target state
cost.sub{2}.width   = 25;
cost.sub{2}.angle   = plant.angi;

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
opt.length = 25;                        % (+): max. number of line searches
opt.MFEPLS = 25;                        % max. number of function evaluations per linesearch
% opt.MSR = 100;                        % max. slope ratio (default=100)

%% Final inits
prevLength = 0;
x = []; y = []; r = [];
fantasy.mean = cell(1,N+J); fantasy.std = cell(1,N+J);
realCost = cell(1,N+J);  latent = cell(1,N+J);  realAcumCost = zeros(1,N+J);
M = cell(N,1); Sigma = cell(N,1); Mcon = cell(N,1); Scon = cell(N,1);
Mfull = cell(N,1); Sfull = cell(N,1);
testLat = cell(N,1);    testCost = cell(N,1);
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