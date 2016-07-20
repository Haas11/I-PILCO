function [mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf] = genInitTrajectory(robot, T0, T1, Tf, T, dt)
% Test trajectory:
t1 = T/3;
t2 = 2*T/3;

q0_est  = ones(1,robot.n)*-0.5;
dq0 = zeros(1,robot.n);    % initial config

fprintf('\nGenerating insertion trajectory...\n');
t1 = T/3;
t2 = 3*T/5;
t3 = 4*T/5;
t_total     = [0:dt:T]'; %#ok<*NBRAK>
t_approach  = [0:dt:t1]';
t_contact1  = [t1:dt:t2]';
t_final    = [t2:dt:T]';

q0 = robot.ikine(T0,q0_est,[1 1 0 0 0 1]);

Tcart1 = ctraj(T0, T1, length(t_approach));           % Cartesian trajectory generation
Tcart2 = ctraj(T1, Tf, length(t_contact1));
Tcart3 = ctraj(Tf, T1, length(t_final));
Ttot = cat(3,Tcart1,Tcart2,Tcart3);

xe_des = [transl(Ttot), tr2eul(Ttot)];          % extract des position
dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration

xe_des   = [t_total, xe_des(1:length(t_total),:)];
dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];

Hf = Ttot(:,:,end);
F0 = zeros(1,6);
xe0 = [transl(T0)', tr2eul(T0)];
dxe0 = zeros(1,6);

mu0 = [q0, dq0, xe0, dxe0, F0, 0]; %[q dq xe dxe Fext t]

S0  = diag([0*ones(1,robot.n) 0*ones(1,robot.n) 0*ones(1,length(F0))].^2); % initial state covariance
end

