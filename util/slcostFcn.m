% S-function cost wrapper

%   cost            cost structure
%     .p            lengths of the 2 pendulums                      [2 x  1 ]
%     .width        array of widths of the cost (summed together)
%     .expl         (optional) exploration parameter
%     .angle        (optional) array of angle indices
%     .target       target state                                    [D x  1 ]
%   m               mean of state distribution                      [D x  1 ]
%   s               covariance matrix for the state distribution    [D x  D ]
%
% *Output arguments:*
%
%   L     expected cost                                             [1 x  1 ]
%   dLdm  derivative of expected cost wrt. state mean vector        [1 x  D ]
%   dLds  derivative of expected cost wrt. state covariance matrix  [1 x D^2]
%   S2    variance of cost                                          [1 x  1 ]

function [sys,x0,str,ts] = slcostFcn(t, x, u, flag, cost, plant)

switch flag   
    case 0  %setup
        [sys,x0,str,ts] = mdlInitializeCostSizes(plant);	% Init
        
    case {3}
        sys = cost.fcn(cost, u', zeros(length(plant.dyno)));
        
    case {1, 2, 4, 9}
        sys = [];
end

function [sys,x0,str,ts]=mdlInitializeCostSizes(plant)

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
sizes.NumOutputs     = 1;%+ length(plant.dyno) + length(plant.dyno)^2 + 1  % 1 + D + D^2 + 1
sizes.NumInputs      = length(plant.dyno);                                 % dyno states
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