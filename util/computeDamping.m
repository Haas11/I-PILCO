function Kd = computeDamping(Kp_env, Kp, Kpo)
% computes damping ratio based on known environment stiffness
% and desired cartesian manipulator stiffness and inherent manipulator inertia.

zeta = 0.7;                         % damping ratio
Kp = diag([Kp(1), Kp(2), Kp(2), Kpo]);
Kp_env = diag(Kp_env);

% M = robot.cinertia(q');              % robot cartesian inertia matrix

Kd = 2*zeta.*sqrt((Kp + Kp_env));  % damping matrix