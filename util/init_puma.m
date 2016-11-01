run mdl_puma560.m
robot = p560;
robot.tool = transl([0 0 0]);
robot.plotopt = {'jaxes','delay',dt_pilco,'trail','y--','zoom',0.75,'scale',0.75};
robot.plotopt3d = {'jaxes','delay',dt_pilco,'scale',1,'noa','base'};
n = robot.n; m = n;
robot.fast = 1;
robot.gravity = [0, -9.81, 0];

dynPert = 0.15;          % [%] Perturbation of dynamics during simulation w.r.t ID
fprintf('\nDynamics perturbed by %2.0f percent.\n',dynPert*100);
perturbedRobot = robot.perturb(dynPert); % perturbed robot model
perturbedRobotNF = perturbedRobot.nofriction('all');
perturbedRobot.gravity = robot.gravity;
perturbedRobotNF.gravity = robot.gravity;