%% 3-link planar manipulator
% this m-file describes the kinematic and dynamic definition of a 2-link
% planar maniplulator using the Robotics Toolbox.

clear; clc; close all;

%% Constants & Globals
dt_sim = 0.02;
dt = dt_sim;

% Environment
Kp_env = [5000, 5000, 0, 0, 0, 0];            %[N/m]  stiffness  (x, y, z, rotx, roty, rotz)
Kd_env = [200, 200, 0, 0, 0, 0];              %[Ns/m] damping

% Impedance parameters:
Kp = [25, 100, 0, 0, 0, 100];          %[N/m]  Stiffness
Kd = [50, 50,  0, 0, 0,  50];          %[Ns/m] Damping

Mdes = diag([1 1 0 0 0 1]);                %[kg]   Inertia matrix
Mdes_inv = pinv(Mdes);                     %[kg^-1]

% Robot model definition
fprintf('\nInitializing robot model...\n\n');
l1 = 0.25;      l2 = 0.25;      l3 = 0.25;          % link lengths
m1 = 1;         m2 = 1;         m3 = 1;             % link masses
Tc1 = [0 0];    Tc2 = [0 0];    Tc3 = [0 0];        % link coulomb friction
Jm1 = 1e-4;     Jm2 = 1e-4;     Jm3 = 1e-4;         % motor inertia
B1 = 1e-4;      B2 = 1e-4;      B3 = 1e-4;          % motor viscous friction

link1 = Revolute('d',0,'a',l1,'alpha',0,...
    'm',m1,'r',[-l1/2, 0, 0],'I',[0, 0, m1*l1/12],'B',B1,'Tc',Tc1,'G',1,'Jm',Jm1);
link2 = Revolute('d',0,'a',l2,'alpha',0,...
    'm',m2,'r',[-l2/2, 0, 0] ,'I',[0, 0, m2*l2/12],'B',B2,'Tc',Tc2,'G',1,'Jm',Jm2);
link3 = Revolute('d',0,'a',l3,'alpha',0,...
    'm',m2,'r',[-l3/2, 0, 0] ,'I',[0, 0, m3*l3/12],'B',B3,'Tc',Tc3,'G',1,'Jm',Jm3);
robot = SerialLink([link1, link2, link3], 'name', '3link-robot',...
    'base',[eye(3), zeros(3,1); zeros(1,3), 1], 'qlim',[-pi, pi; -pi, pi; -pi, pi],...
    'tool',[eye(3), zeros(3,1); zeros(1,3), 1]);
robot.gravity = [0, -9.81, 0];
robot.plotopt = {'jaxes','delay',dt,'trail','y--','zoom',0.75,'scale',0.75};
n = robot.n; m = 3;
robot %#ok<NOPTS>

perturbedRobot = robot.perturb(0.10); % perturbed robot model
perturbedRobot.gravity = robot.gravity;

% 2. Set up the scenario
T = 10.0;                                       % [s] prediction horizon
peg=1;
xhole = [0.5, 0.3, 0];   % center hole location [x, y, phi/z]
xc    = [0.45, 2, 2, 10, 10, 10]';  % [m] environment constraint location
x0    = [0.4 0 0];
T0    = transl(x0);      % start pose end-effector
T11   = transl([0.5 0 0]);
[mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf]...
    = genTrajectory(robot, peg, T0, T11, xhole, xc, T, dt);             
fprintf('\nFinal transformation: \nH^0_n(T) = \n\n');

% disp(Hf);
% robot.plot(mu0(1:robot.n)');


%% Simulate

% ideal simulation:
fprintf('\nSimulating perfect impedance control...\n')
simOut = sim('CartesianImpControl','SaveState','on','SaveOutput','on');
XYPhi_sim = get(simOut,'XYPhi');
q_sim = get(simOut,'q_sim');

% perturbed simulation:
fprintf('\nSimulating perturbed impedance control...\n')
simOutPert = sim('perturbedCartImpControl','SaveState','on','SaveOutput','on');
q_sim_pert = get(simOutPert,'q_sim');
XYPhi_sim_pert = get(simOutPert,'XYPhi');

% plot result
figure;
plot(XYPhi_sim.Time, XYPhi_sim.Data(:,1),'r');
hold on; grid on;
plot(XYPhi_sim_pert.Time, XYPhi_sim_pert.Data(:,1),'b');
legend('ideal','perturbed');
hold on;
plot(XYPhi_sim.Time, XYPhi_sim.Data(:,2),'r');
hold on;
plot(XYPhi_sim_pert.Time, XYPhi_sim_pert.Data(:,2),'b');
hold on;
plot(XYPhi_sim.Time, XYPhi_sim.Data(:,3),'r');
hold on;
plot(XYPhi_sim_pert.Time, XYPhi_sim_pert.Data(:,3),'b');
title('Ideal vs Perturbed Inverse Dynamics');
xlabel('Time [s]'); ylabel('Cartesian position');

% Decimate data for video output:
dt_anim = 0.1;
ratio = dt_anim/dt;
q_simd=downsample(q_sim,ratio);
robot.plotopt = {'jaxes','delay',dt_anim,'trail','y--','zoom',0.75,'scale',0.75};
% render & store animation:
animateVideo(robot, q_simd, 'robotVid')

