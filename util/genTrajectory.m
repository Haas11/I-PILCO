function [mu0, S0, xe_des, dxe_des, ddxe_des, Hf, Rd, Hd] = genTrajectory(robot, mode, H0, H1, H2, H3, xhole, xc, T,  dt)
global vert fac
vert = ones(8,3); vert(1,1) = xc(1); vert(4:5,1) = xc(1); vert(8,1) = xc(1);
vert(1:2,2:3) = -1; vert(3:4,3) = -1; vert(5:6,2) = -1;
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

if isscalar(robot)
    n=robot;
else
    n=n;
end

% Test trajectory:
t1 = 2*T/5;
t2 = 2*T/4;
t3 = 3*T/4;
t_total     = (0:dt:T)';
t_approach  = [0:dt:t1]';
t_contact1  = (t1:dt:t2)';
t_contact2  = (t2:dt:t3)';
t_final     = (t3:dt:T)';

q0_est  = ones(1,n)*-0.5;
dq0 = zeros(1,n);    % initial config

xhole2 = xhole + [0.075 0 0];
Thole1 = transl(xhole); Thole2 = transl(xhole2);

if mode==0
    fprintf('Generating compliance test trajectory...\n');
    if n==6
        q0 = robot.ikine6s(H0);
    else
        q0 = robot.ikine(H0,q0_est);
    end
    H0 = robot.fkine(q0);
    robot.plot(q0);
    patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
    
    Hf = transl(xhole);
    
    [xe1, ~, ~] = mtraj(@tpoly, transl(H0)', transl(H1)', t_approach);
    [xe2, ~, ~] = mtraj(@tpoly, transl(H1)', transl(H2)', t_contact1);
    [xe3, ~, ~] = mtraj(@tpoly, transl(H2)', transl(H3)', t_contact2);
    [xe4, ~, ~] = mtraj(@tpoly, transl(H3)', transl(Hf)', t_final);
    
    Tcart1 = transl(xe1);    Tcart2 = transl(xe2);    Tcart3 = transl(xe3);    Tcart4 = transl(xe4);
    Ttot = cat(3,Tcart1,Tcart2,Tcart3, Tcart4);
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position  (xyz)
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
        q0 = robot.ikine6s(H0);
    elseif n < 6
        q0 = robot.ikine(H0,q0_est,[1 1 0 0 0 1],'pinv');  % numerical inverse kinematics
    elseif n > 3
        [q0,~] = robot.ikcon(H0);
    end
    H0 = robot.fkine(q0);
    
    %     Thole1 = T1; Thole2=T1;
    [xe1, ~, ~] = mtraj(@tpoly, transl(H0)', transl(H1)', t_approach);
    [xe2, ~, ~] = mtraj(@tpoly, transl(H1)', transl(Thole1)', t_contact1);
    [xe3, ~, ~] = mtraj(@tpoly, transl(Thole1)', transl(Thole2)', t_insert);
    [xe4, ~, ~] = mtraj(@tpoly, transl(Thole2)', transl(Thole2)', t_rest);
    
    Tcart1 = transl(xe1);
    Tcart2 = transl(xe2);
    Tcart3 = transl(xe3);
    Tcart4 = transl(xe4);
    
    Ttot = cat(3,Tcart1,Tcart2,Tcart3,Tcart4);
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position
    dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
    ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration
    
    xe_des   = [t_total, xe_des(1:length(t_total),:)];
    dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
    ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];
elseif mode==3
    
    t1 = T/6;
    t2 = 3*T/6;
    t3 = 5*T/6;
    t_total     = [0:dt:T]'; %#ok<*NBRAK>
    t_approach  = [0:dt:t1]';
    t_contact1  = [t1:dt:t2]';
    t_contact2  = (t2:dt:t3)';
    t_rest      = [t3:dt:T]';
    
%     q0 = [0; 0.25; 0; -1.43; 0; 1.46; 0];   % [0.45 0 0.4]         
    q0 = [0; 11.75; 0; -100; 0; 68; 0];     % [0.45 0 0.3]
    
    xe0 = [transl(H0)', tform2quat(H0)];
    dxe0 = zeros(1,7);
    F0 = zeros(1,6);
    
    initRot = t2r(quat2tform([0 1 0 0]));
    firstRot = t2r(quat2tform([0.1 0.9 0.1 0.1]));
    
    %     T1 = rt2tr(initRot,[0.5 0 0.05]');
    %     T2 = rt2tr(initRot,[0.5 0 0.05]');
    %     T3 = rt2tr(initRot,[0.6 0 0.05]');
    
    [xe1, ~, ~] = mtraj(@tpoly, transl(H0)', transl(H1)', t_approach);
    [xe2, ~, ~] = mtraj(@tpoly, transl(H1)', transl(H2)', t_contact1);
    [xe3, ~, ~] = mtraj(@tpoly, transl(H2)', transl(H3)', t_contact2);
    [xe4, ~, ~] = mtraj(@tpoly, transl(H3)', transl(H3)', t_rest);
    
    Tcart1 = transl(xe1);
    Tcart2 = transl(xe2);
    Tcart3 = transl(xe3);
    Tcart4 = transl(xe4);
    
    Ttot = cat(3,Tcart1,Tcart2,Tcart3,Tcart4);
    Ttot(1:3,1:3,:) = repmat(initRot,1,1,size(Ttot,3));
    
    mu0 = [xe0, dxe0, F0, 0]; %[q dq xe dxe Fext t]
    S0  = diag([1e-3*ones(1,7), 1e-5*ones(1,7) 0.1*ones(1,length(F0)), 1e-5].^2); % initial state covariance
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position
    dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
    ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration
    
    xe_des   = [t_total, xe_des(1:length(t_total),:)];
    dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
    ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];
elseif mode==2
    
    t1 = T/5;
    t2 = 4*T/5;
    t_total     = [0:dt:T]'; %#ok<*NBRAK>
    t_approach  = [0:dt:t1]';
    t_contact   = [t1:dt:t2]';
    t_rest      = [t2:dt:T]';
    
%     q0 = [0; 0.25; 0; -1.43; 0; 1.46; 0];   % [0.45 0 0.4]         
    q0 = [0; 11.75; 0; -100; 0; 68; 0];     % [0.45 0 0.3]
    
    xe0 = [transl(H0)', tform2quat(H0)];
    dxe0 = zeros(1,7);
    F0 = zeros(1,6);
    
    initRot = t2r(quat2tform([0 1 0 0]));
    firstRot = t2r(quat2tform([0.1 0.9 0.1 0.1]));
    
    %     T1 = rt2tr(initRot,[0.5 0 0.05]');
    %     T2 = rt2tr(initRot,[0.5 0 0.05]');
    %     T3 = rt2tr(initRot,[0.6 0 0.05]');
    
    [xe1, ~, ~] = mtraj(@tpoly, transl(H0)', transl(H1)', t_approach);
    [xe2, ~, ~] = mtraj(@tpoly, transl(H1)', transl(H2)', t_contact);
    [xe4, ~, ~] = mtraj(@tpoly, transl(H2)', transl(H2)', t_rest);
    
    Tcart1 = transl(xe1);
    Tcart2 = transl(xe2);
    Tcart4 = transl(xe4);
    
    Ttot = cat(3,Tcart1,Tcart2,Tcart4);
    Ttot(1:3,1:3,:) = repmat(initRot,1,1,size(Ttot,3));
    
    mu0 = [xe0, dxe0, F0]; %[q dq xe dxe Fext t]
    S0  = diag([0.001*ones(1,7) 0*ones(1,7) 0.1*ones(1,length(F0))].^2); % initial state covariance
    
    xe_des = [transl(Ttot), tr2rpy(Ttot)];          % extract des position
    dxe_des = [zeros(1,6); diff(xe_des)/dt];        % differentiate for des velocity
    ddxe_des = [zeros(1,6); diff(dxe_des)/dt];      % second difference for des acceleration
    
    xe_des   = [t_total, xe_des(1:length(t_total),:)];
    dxe_des  = [t_total, dxe_des(1:length(t_total),:)];
    ddxe_des = [t_total, ddxe_des(1:length(t_total),:)];
end

Hf = Ttot(:,:,end);
F0 = zeros(1,6);
xe0 = [transl(H0)', tr2rpy(H0)];
dxe0 = zeros(1,6);

if mode~=2 && mode~=3
    mu0 = [q0, dq0, xe0, dxe0, F0, 0]; %[q dq xe dxe Fext t]
    S0  = diag([0*ones(1,n) 0*ones(1,n) 0*ones(1,length(F0))].^2); % initial state covariance
end

Hd = Ttot;
Rd = tr2rt(Ttot);

end

