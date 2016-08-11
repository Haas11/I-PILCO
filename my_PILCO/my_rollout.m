%% rollout.m
% *Summary:* Generate a state trajectory using a simulink model (and any additional
% dynamics) from a particular initial state by applying either a particular
% policy or random actions.
%
%   function [x y L latent] = my_rollout(start, policy, H, plant, robot)
%
% *Input arguments:*
%
%   start       vector containing initial states (without controls)   [nX  x  1]
%   policy      policy structure
%     .fcn        policy function
%     .p          parameter structure (if empty: use random actions)
%     .maxU       vector of control input saturation values           [nU  x  1]
%   H           rollout horizon in steps
%   plant       the dynamical system structure
%     .subplant   (opt) additional discrete-time dynamics
%     .augment    (opt) augment state using a known mapping
%     .constraint (opt) stop rollout if violated
%     .poli       indices for states passed to the policy
%     .dyno       indices for states passed to cost
%     .odei       indices for states passed to the ode solver
%     .subi       (opt) indices for states passed to subplant function
%     .augi       (opt) indices for states passed to augment function
%   robot       robot model
%
% *Output arguments:*
%
%   x          matrix of observed states                           [H   x nX+nU]
%   y          matrix of corresponding observed successor states   [H   x   nX ]
%   L          cost incurred at each time step                     [H   x   1  ]
%   latent     matrix of latent states                             [H+1 x   nX ]
%   ref        matrix of references                                [H+1 x   nX ]
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-05-21
%
% Modified by Victor van Spaandonk, 2016.
%
%% High-Level Steps
%
% # Compute control signal $u$ from state $x$:
% either apply policy or random actions
% # Simulate the true dynamics for one time step using the current pair $(x,u)$
% # Check whether any constraints are violated (stop if true)
% # Apply random noise to the successor state
% # Compute cost (optional)
% # Repeat until end of horizon



function [x, y, L, latent, ref] = my_rollout(mu0, policy, H, plant, robot)
%% Code
global initTrial

% if isempty(initialTrial)
%     initialTrial = true;
% end

augi = plant.augi;
if isfield(plant,'subplant')
    subi = plant.subi;
else
    plant.subplant = inline('[]',1);  %#ok<DINLN>
    subi = [];
end

odei = plant.odei;
poli = plant.poli;  %#ok<*NASGU>
dyno = plant.dyno;
angi = plant.angi;
simi = sort([odei subi]);
dt = plant.dt;

% nX = length(simi)+length(augi);
nX = length(plant.indices);
nU = length(policy.maxU);
nA = length(angi);

% Simulate Dynamical Model using Simulink
goodSim = false;
while ~goodSim
    x0trial = [mu0(plant.dyno([1 2])), 0]';
    clear ref dref
    if initTrial
        % do nothing
    else
        x0trial = x0trial + plant.startStateInterval(1:3).*(2*rand(3,1)-1);  % perturb by uniform random amount [-S0Trial, +S0Trial]
    end
    Hd0 = transl(x0trial);
    
    % Inverse Kinematics for starting angles:
    if x0trial(2)<0
        q0est = [-pi/2 pi/2 0];
    else
        q0est = [pi/2 -pi/2 0];
    end
    q0 = robot.ikcon(Hd0,q0est);     % w/ joint limits
    %         q0 = robot.ikine(Hd0,q0est,[1 1 0 0 0 1]);     % without joint limits
    Hd0 = robot.fkine(q0);
    x0trial = transl(Hd0);
    
    assignin('base','x0trial',x0trial);
    assignin('base','q0',q0);
    assignin('base','Hd0',Hd0);
    
    simOut = sim(plant.rollout_model,'SaveOutput','on','SimulationMode','normal');
    x       = simOut.get('x');      % state vector  [H+1  x  nX]
    latent  = simOut.get('latent'); % latents       [H+1  x  nX+2*nU]
    a       = simOut.get('a');      % Actions       [H  x  nU]
    L       = simOut.get('L');      % cost          [H  x  1]
    ref     = simOut.get('ref');    % reference     [H+1  x  6]
    
    x = squeeze(x);    % reshape due to odd simulink quirk
    ref = squeeze(ref);
    a = squeeze(a);
    latent = squeeze(latent);
    
    if size(x,2) ~= nX      % ugly debug for simulink version issues
        x = x';
    end
    if size(a,2) ~= nU
        a = a';
    end
    if size(latent,2) ~= nX+nU
        latent = latent';
    end
    
    samPerSec = ceil(1/dt);
    Ha = size(x,1);
    if Ha ~= H+1            % if slight mismatch
        H1 = Ha-1;           % set straight
        if H1 < H-5          % if simulation aborted
            H = max(H1 - 0.5*samPerSec,1);     % delete last 0.5 second of data
        else
            H = H1;          % use full dataset
        end
    end
    
    if H > 10
        goodSim = true;
    end
end
% Prepare outputs:
dref = [zeros(1,6); ((ref(2:H+1,:) - ref(1:H,:))./dt)];
ref = [ref(1:H+1,:), dref, zeros(H+1,6), (0:dt:H*dt)'];

y = x(2:H+1,1:nX);          % successor states  [H x nX] (targets)
x = [x(1:H,:) a(1:H,:)];    % observed  states  [H x nX+2*nU]

latent(H+1,1:nX) = latent(end,1:nX);  % noiseless state
latent           = latent(1:H+1,:);

L = L(1:H,1);
end
