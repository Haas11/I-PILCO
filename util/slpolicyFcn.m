

function [sys, x0, str,ts] = slpolicyFcn(t, x, u, flag, policy, plant)

switch flag   
    case 0  %setup
        [sys, x0, str, ts] = mdlInitializePolicySizes(policy, plant);	% Init
        
    case {3}       
        sys = policy.fcn(policy, u, zeros(length(plant.poli)));
        %sys = conCat_det(@congp,@gSat, policy, u, zeros(length(plant.poli)));
        
    case {1, 2, 4, 9}
        sys = [];
end

function [sys, x0, str, ts] = mdlInitializePolicySizes(policy, plant)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = length(policy.maxU);                                       
sizes.NumInputs      = length(plant.poli);                                  % poli states
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [-1 0];
%end mdlInitializeSizes