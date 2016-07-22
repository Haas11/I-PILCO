function [mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf, Rd, Td] = genTrajectory(robot, mode, T0, T1, xhole, xc, T,  dt)
global vert fac
vert = ones(8,3); vert(1,1) = xc(1); vert(4:5,1) = xc(1); vert(8,1) = xc(1);
vert(1:2,2:3) = -1; vert(3:4,3) = -1; vert(5:6,2) = -1;
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

% Test trajectory:
t1 = T/3;
t2 = 2*T/4;
t3 = 3*T/4;
t_total     = (0:dt:T)';
t_approach  = [0:dt:t1]';
t_contact1  = (t1:dt:t2)';
t_contact2  = (t2:dt:t3)';
t_final     = (t3:dt:T)';

q0_est  = ones(1,robot.n)*-0.5;
dq0 = zeros(1,robot.n);    % initial config

xhole2 = xhole + [0.05 0 0];
Thole1 = transl(xhole); Thole2 = transl(xhole2);

if mode==0
    fprintf('Generating compliance test trajectory...\n');
    T0 = transl(0.4, 0, 0); q0 = robot.ikine(T0,q0_est,[1 1 0 0 0 1]);
    robot.plot(q0); patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
    
    T1 = transl(0.5, 0, 0);                         % end config in task space
    T2 = transl(0.5, 0.15, 0);
    T3 = transl(0.5, 0.05, 0);
    Tf = transl(0.5, 0.1, 0);
    
    Tcart1 = ctraj(T0,T1,length(t_approach));           % Cartesian trajectory generation
    Tcart2 = ctraj(T1, T2, length(t_contact1));
    Tcart3 = ctraj(T2, T3, length(t_contact2));
    Tcart4 = ctraj(T3, Tf, length(t_final));
    Ttot = cat(3,Tcart1,Tcart2,Tcart3, Tcart4);
    xe_des = [transl(Ttot), tr2eul(Ttot)];        % extract position/pose
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position
    dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
    ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration
    
    xe_des   = [t_total, xe_des(1:length(t_total),:)];
    dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
    ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];
    
elseif mode==1
    fprintf('\nGenerating peg insertion trajectory...\n');
    t1 = T/3;
    t2 = 3*T/5;
    t3 = 4*T/5;
    t_total     = [0:dt:T]'; %#ok<*NBRAK>
    t_approach  = [0:dt:t1]';
    t_contact1  = [t1:dt:t2]';
    t_insert    = [t2:dt:t3]';
    t_rest      = [t3:dt:T]';
    
    if robot.isspherical
        q0 = robot.ikine6s(T0);
    elseif robot.n < 6
        q0 = robot.ikine(T0,q0_est,[1 1 0 0 0 1],'pinv');  % numerical inverse kinematics
    elseif robot.n > 3
        [q0,~] = robot.ikcon(T0);
    end
    
    %T1 = transl(xhole(1), 0, 0);                          % next config in task space
    
    Tcart1 = ctraj(T0, T1, length(t_approach));           % Cartesian trajectory generation
    Tcart2 = ctraj(T1, Thole1, length(t_contact1));
    Tcart3 = ctraj(Thole1, Thole2, length(t_insert));
    Tcart4 = ctraj(Thole2, Thole2, length(t_rest));
    
    Ttot = cat(3,Tcart1,Tcart2,Tcart3,Tcart4);
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position
    dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
    ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration
    
    xe_des   = [t_total, xe_des(1:length(t_total),:)];
    dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
    ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];
elseif mode==2
    t_01 = T/5/dt;
    t_12 = T/5/dt;
    t_23 = T*2/5/dt;
    t_33 = T/5/dt;
    
    initRot = t2r(quat2tform([0 1 0 0]));
    firstRot = t2r(quat2tform([0.1 0.9 0.1 0.1]));
    
    T0 = rt2tr(initRot,[0.5 0 0.4]');
    T1 = rt2tr(firstRot,[0.5 0 0.2]');
    T2 = rt2tr(initRot,[0.6 0 0.1]');
    T3 = rt2tr(initRot,[0.6 0.1 0.1]');
    
    Tcart1 = ctraj(T0, T1, t_01);           % Cartesian trajectory generation
    Tcart2 = ctraj(T1, T2, t_12);
    Tcart3 = ctraj(T2, T3, t_23);
    Tcart4 = ctraj(T3, T3, t_33);
    
    Ttot = cat(3,Tcart1,Tcart2,Tcart3,Tcart4,Tcart4);
           
    Href1 = cat(3,T0,T1,T2,T3);
end

Hf = Ttot(:,:,end);
F0 = zeros(1,6);
xe0 = [transl(T0)', tr2rpy(T0)];
dxe0 = zeros(1,6);

mu0 = [q0, dq0, xe0, dxe0, F0, 0]; %[q dq xe dxe Fext t]

S0  = diag([0*ones(1,robot.n) 0*ones(1,robot.n) 0*ones(1,length(F0))].^2); % initial state covariance

Td = Ttot;
Rd = tr2rt(Ttot);


% Desired angular acceleration
% dT_d = delta2tr(xe_des);
% dEUL_d =

% domega_d = T_d*ddEUL_d + dT_d*dEUL_d;   % desired angular acceleration


end

