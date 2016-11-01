function [t_train, t_learn, t_pred] = timePILCOT_ST(nc, poli, nii, x, y, r, mu0Sim, S0Sim)
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

opt.verbosity = 0;                      % optimization verbosity      [0-3]
plotting.verbosity = 2;                 % plotting verbosity          [0-3]

%% 1. Define state indices
stateLength = 25;
indices = 1:1:stateLength;
n=3;
odei = 1:1:2*n+6;
augi    = [];                                           % augi  indicies for variables augmented to the ode variables
angi    = [];                                           % angi  indicies for variables treated as angles (using sin/cos representation) (subset of indices)
dyno    = [7 8 13 14 19 20];                          % dyno  indicies for the output from the dynamics model and indicies to loss    (subset of indices)
dyni    = [1 2 3 4];                              % dyni  indicies for inputs to the dynamics model                               (subset of dyno)
difi    = [1 2 3 4];                              % difi  indicies for training targets that are differences                      (subset of dyno)
refi    = [1 2 3 4];                                 % indices for which to encode a reference as  prior mean
REF_PRIOR  = 0;

ref_select = [1 2 7 8 13 14];                           % indices of reference corresponding to dyno    [xe dxe F]

dynoTitles = stateNames(indices(dyno));
actionTitles = {'Kp_x  [N/m]', 'Kp_{y}  [N/m]', 'x_{ref}','y_{ref}'};%, 'Kp_{rot} [Nm/rad]'};

%% 3. Set up the plant structure
dt_pilco = 0.1;          % [s] PILCO sampling rate 
dt = 0.005;
H = 100;

xhole = [0.7, 0.2, 0];   % center hole location [x, y, phi/z]

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
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  % linear saturating controller
maxVel = 0.25;         % maximum Cartesian velocity
policy.maxU  = [250/2 250/2, dt*maxVel  dt*maxVel];
policy.minU  = [10    10,   -dt*maxVel -dt*maxVel];
policy.impIdx = [1, 2]; 			% non-negative indices of policy outputs (saturate + translate)
policy.refIdx = [3, 4];             % reference indices (only saturate)
policy.SNR = 100;
Du = length(policy.maxU);
translVec = [ones(size(policy.impIdx)).*2, ones(size(policy.refIdx))];

% GP Controller:
policy.fcn = @(policy,m,s)my_mixedConCat(@congp,@my_mixedGSat,policy,m,s);  % linear saturating controller
policy.p.targets = 0.01*randn(nc, length(policy.maxU));    % init. policy targets
policy.p.inputs  = gaussian(zeros(1,length(poli)), diag(ones(1,length(poli))*0.01), nc)';                % policy pseudo inputs   [ N  x  d ]
policy.p.hyp = ...                                                                             % GP-log hyperparameters [(d+2) x  D ]
    repmat(log([ones(1,length(poli))*1, 1, 0.01]'), 1, length(policy.maxU));

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

cost.sub{1}.fcn     = @my_lossSat;
cost.sub{1}.losi    = [1 2];                        % indicies for saturating cost states
cost.sub{1}.target  = (xhole(1:2) + [0.05 0])';   % target state
cost.sub{1}.width   = 0.05;
cost.sub{1}.angle   = plant.angi;

cost.sub{2}.fcn     = @my_lossSat;
cost.sub{2}.losi    = 5;                        % indicies for saturating cost states
cost.sub{2}.target  = 0;   % target state
cost.sub{2}.width   = 20;
cost.sub{2}.angle   = plant.angi;

%% 6. Set up the GP dynamics model structure
dynmodel.fcn    = @my_gp1d;                    % function for GP predictions
dynmodel.train  = @my_train;                % function to train dynamics model
dynmodel.induce = zeros(nii,0,1);           % shared/individual inducing inputs per target dim (sparse GP)
dynmodel.parallel = false;                  % train individual target dimensions in parellel
trainOpt        = [100 100];                % max. number of line searches [full, sparse]

%% 7. Parameters for policy optimization
opt.fh = 1;
opt.method = 'BFGS';                    % 'BFGS' (default), 'LBFGS' (x>1000), 'CG'
opt.length = 50;                        % (+): max. number of line searches
opt.MFEPLS = 25;                        % max. number of function evaluations per linesearch
% opt.MSR = 100;                        % max. slope ratio (default=100)

%% Load dataset
% load('dataForTiming.mat');
x = x(1:500,:); y = y(1:500,:); r = r(1:500,:);
% dynmodel.ref = r(1:H,:);

%% 3. Controlled learning (N iterations)
tic
timedModelLearning;
t_train = toc;

tic
[policy.p, ~] = minimize(policy.p, 'my_value', opt, mu0Sim', S0Sim, ...
    dynmodel, policy, plant, cost, H);
t_learn = toc;

tic
my_pred(policy, plant, dynmodel, mu0Sim', S0Sim, H);
t_pred = toc;
   
end

