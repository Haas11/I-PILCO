l1 = 0.4;       l2 = 0.4;       l3 = 0.4;          % link lengths
m1 = 1;      m2 = 1;      m3 = 1;             % link masses
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
    'base',[eye(3), zeros(3,1); zeros(1,3), 1], 'qlim',[-1, 1; -1, 1; -1, 1]*99*pi/100,...
    'tool',[eye(3), zeros(3,1); zeros(1,3), 1]);
robot.gravity = [0, -9.81, 0];
robot.plotopt = {'jaxes','delay',dt_pilco,'trail','y--','zoom',0.75,'scale',0.75};
n = robot.n; m = 3;
robot.fast = 1;

fprintf('\nDynamics perturbed by %2.2f percent.\n',dynPert);    % perturb robot by dynPert [%]
perturbedRobot = my_perturb(robot,dynPert); % perturbed robot model
perturbedRobot.gravity = robot.gravity;

perturbedRobotNF = perturbedRobot;%.nofriction('all');
perturbedRobotNF.gravity = robot.gravity;