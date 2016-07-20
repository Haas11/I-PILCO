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
%   L          cost incurred at each time step                     [ 1  x    H ]
%   latent     matrix of latent states                             [H+1 x   nX ]
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



function [x, y, L, latent, ref] = my_rollout(start, policy, H, plant, robot)
%% Code
augi = plant.augi;
% if isfield(plant,'augment')
%     augi = plant.augi;             % sort out indices!
% else
%     plant.augment = inline('[]'); 
%     augi = []; 
% end

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
q0 = start(1:robot.n);
dq0 = start(robot.n+1:2*robot.n);
Te0 = robot.fkine(q0);      trans0 = transl(Te0);          Ang0 = tr2eul(Te0);
F0 = zeros(1,6);
% x0 = [q0', dq0', trans0', Ang0, F0];
% a0 = policy.fcn(policy, x0(dyno(poli)), zeros(length(poli)));

simOut = sim(plant.rollout_model,'SaveOutput','on','SimulationMode','normal');
x       = simOut.get('x');      % state vector  [H+1  x  nX]
latent  = simOut.get('latent'); % latents       [H+1  x  nX+2*nU]
a       = simOut.get('a');      % Actions       [H  x  nU]
L       = simOut.get('L');      % cost          [H  x  1]
ref     = simOut.get('ref');    % reference     [H  x  6]

x = squeeze(x)';    % reshape due to odd simulink quirk
ref = squeeze(ref);
a = squeeze(a)';
latent = squeeze(latent)';
% L = squeeze(L);

samPerSec = ceil(1/dt);
Ha = size(x,1);
if Ha ~= H+1            % if slight mismatch
   H1 = Ha-1;           % set straight 
   if H1 < H-5          % if simulation aborted
       H = max(H1 - 2*samPerSec,1);     % delete last 2 seconds of data
   else 
       H = H1;          % use full dataset
   end
end

dref = (ref(2:H+1,:) - ref(1:H,:))./dt;
ref = [ref(2:H+1,:), dref, zeros(H,6)];    % reference targets

y = x(2:H+1,1:nX);          % successor states  [H x nX] (targets)
x = [x(1:H,:) a(1:H,:)];    % observed  states  [H x nX+2*nU]

latent(H+1,1:nX) = latent(end,1:nX);  % noiseless state
latent           = latent(1:H+1,:); 

L = L(1:H,1);