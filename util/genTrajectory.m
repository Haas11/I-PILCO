function [mu0, S0, xe_des, dxe_des, ddxe_des, T, Hf, Rd, Td] = genTrajectory(robot, peg, H0, H1, H2, H3, xhole, xc, T,  dt)
global vert fac
vert = ones(8,3); vert(1,1) = xc(1); vert(4:5,1) = xc(1); vert(8,1) = xc(1);
vert(1:2,2:3) = -1; vert(3:4,3) = -1; vert(5:6,2) = -1;
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
% if ~ishandle(5)         % robot animation
%     figure(5);
%     set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9])
% else
%     set(0,'CurrentFigure',5);
% end
% patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');

% Test trajectory:
t1 = 4*T/10;
t2 = 7*T/10;
t3 = 9*T/10;
t_total     = (0:dt:T)';
t_approach  = [0:dt:t1]';
t_contact1  = (t1:dt:t2)';
t_contact2  = (t2:dt:t3)';
t_final     = (t3:dt:T)';

% q0_est  = ones(1,robot.n)*-0.5;
x0 = transl(H0);
if x0(2)<0
    q0_est = [-pi/2 pi/2 0];
else
    q0_est = [pi/2 -pi/2 0];
end

dq0 = zeros(1,robot.n);    % initial config

xhole2 = xhole + [0.075 0 0];
Thole1 = transl(xhole); Thole2 = transl(xhole2);

if ~peg
    fprintf('Generating compliance test trajectory...\n');
    if robot.n==6
        q0 = robot.ikine6s(H0);
    else
%         q0 = robot.ikine(H0,q0_est,[1 1 0 0 0 1],'pinv');
        q0 = robot.ikcon(H0,q0_est);%,[1 1 0 0 0 1],'pinv');
    end
    
    H0 = robot.fkine(q0);   
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
    
elseif peg
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
    elseif robot.n < 6
        q0 = robot.ikine(H0,q0_est,[1 1 0 0 0 1],'pinv');  % numerical inverse kinematics
    elseif robot.n < 6 && jointLimit
        [q0,~] = robot.ikcon(H0);
    end
    H0 = robot.fkine(q0);
    
    %T1 = transl(xhole(1), 0, 0);                          % next config in task space
    
    % @lspb  = Trapizoidal function
    % @tpoly = 5th-order polynomial function
    %
    %     Tcart1 = ctraj(T0, T1, length(t_approach));           % Cartesian trajectory generation
    %     Tcart2 = ctraj(T1, Thole1, length(t_contact1));
    %     Tcart3 = ctraj(Thole1, Thole2, length(t_insert));
    %     Tcart4 = ctraj(Thole2, Thole2, length(t_rest));
    
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
end

Hf = Ttot(:,:,end);
F0 = zeros(1,6);
xe0 = [transl(H0)', tr2rpy(H0)];
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

