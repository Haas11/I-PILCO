function [x, y, L, latent, r] = my_iiwaRollout(policy, plant, cost, H, Href, init)
%
% Inputs:
%   - handles   handles for ROS communication
%   - policy    controller
%   - plant     plant structure
%   - H         Number of timesteps
%   - Href      Reference
%
% Outputs:
%   - x         observed states         [H x nX+nU]
%   - y         successor states        [H x nX]
%   - L         Cost                    [H x 1]
%   - latent    latent states           [H x nX]
%   - r         reference mean          [H+1 x nX]
% Code by Victor van Spaandonk
% July 2016


%% Inits
if nargin > 5
    initialRollout = true;  % precomputed stiffness trajectory
else
    initialRollout = false;
end

poli = plant.poli;
dyno = plant.dyno;

nX = 7+7+6;                 % number of states = xe + dxe + F
nU = length(policy.maxU);   % number of actions

x = zeros(H+1, nX);
a = zeros(H,nU);
L = zeros(H,1);
r = zeros(H,nX);

pose1Send = false;
pose2Send = false;
pose3Send = false;

% Initialize the ROS infrastructure:
bufferSize = 10;
relativeVelocity = 0.1;
handles = initROSiiwa(bufferSize, plant.dt, relativeVelocity, 'slip');

poseCommandMsg = handles.poseCommandMsg;
jointCommandMsg = handles.jointCommandMsg;

% state vector:
%x = zeros(handles.N,18);    % state vector [x dx W]

% Configure Cartesian Impedance Mode:
handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.CARTESIANIMPEDANCE;

handles.config.Mode.CartesianStiffness.Stiffness.X = 2500;
handles.config.Mode.CartesianStiffness.Stiffness.Y = 2500;
handles.config.Mode.CartesianStiffness.Stiffness.Z = 2500;
handles.config.Mode.CartesianStiffness.Stiffness.A = 300;
handles.config.Mode.CartesianStiffness.Stiffness.B = 300;
handles.config.Mode.CartesianStiffness.Stiffness.C = 300;

handles.config.Mode.CartesianDamping.Damping.X = 0.7;
handles.config.Mode.CartesianDamping.Damping.Y = 0.7;
handles.config.Mode.CartesianDamping.Damping.Z = 0.7;
handles.config.Mode.CartesianDamping.Damping.A = 0.7;
handles.config.Mode.CartesianDamping.Damping.B = 0.7;
handles.config.Mode.CartesianDamping.Damping.C = 0.7;

handles.config.Mode.NullspaceStiffness = 50;
handles.config.Mode.NullspaceDamping = 0.7;


try
    
    fprintf('\n --- Moving to start state --- \n');
    pause(2);
    
    % Initialize the pose of the robot:
    q0 = [-0.0012;
        0.2596;
        0.0014;
        -1.4267;
        -0.0003;
        1.4555;
        0.0002];
    
    jointCommandMsg.Position.A1 = q0(1);
    jointCommandMsg.Position.A2 = q0(2);
    jointCommandMsg.Position.A3 = q0(3);
    jointCommandMsg.Position.A4 = q0(4);
    jointCommandMsg.Position.A5 = q0(5);
    jointCommandMsg.Position.A6 = q0(6);
    jointCommandMsg.Position.A7 = q0(7);
    send(handles.pubJointPosition, jointCommandMsg);
    
    %     pos = transl(Href(:,:,1));              % initial position
    %     quat = Quaternion(Href(:,:,1)).double;  % initial orientation
    %     xe_0 = [pos; quat'];                    % initial state
    %     poseCommandMsg.Pose.Position.X = xe_0(1);
    %     poseCommandMsg.Pose.Position.Y = xe_0(2);
    %     poseCommandMsg.Pose.Position.Z = xe_0(3);
    %     poseCommandMsg.Pose.Orientation.X = xe_0(4);
    %     poseCommandMsg.Pose.Orientation.Y = xe_0(5);
    %     poseCommandMsg.Pose.Orientation.Z = xe_0(6);
    %     poseCommandMsg.Pose.Orientation.W = xe_0(7);
    %     send(handles.pubCartesianPose, poseCommandMsg);
    
    % Switch to Cartesian Impedance mode:
    response = call(handles.client, handles.config, 'Timeout', 5);
    pause(2);
    if response.Success
        fprintf('\n --- Cartesian Impedance Mode activated --- \n');
        handles.config.Mode.Mode = -1;      % indicate that control mode is the same (time saving)
    else
        error('Failed to start Cartesian impedance mode.: %s', response.Error);
    end
    
    %% Perform Rollout
%     validAnswer = false;
%     while ~validAnswer
%         reply = input('\nReady to start rollout? [y/n]...','s');
%         
%         if strcmpi(reply,'y')
%             validAnswer = true;
%         elseif strcmpi(reply,'n')
%             error('Rollout Aborted by User');
%         else
%             warning('please hit "y" or "n"');
%         end
%     end
    
    fprintf('Starting rollout in...\n');
    pause(1);
    fprintf('\n5 \n'); pause(1); fprintf('\n4 \n'); pause(1); fprintf('\n3 \n');
    pause(1); fprintf('\n2 \n'); pause(1); fprintf('\n1 \n'); pause(1);
    
    reset(handles.r);
    
    % loop
    t = rostime('now');
    startTime = double(t.Sec) + double(t.Nsec)/1e9;
    time = zeros(H,1);
    
    for i=1:H
        
        % read state information
        [time(i), pose, wrench] = readRobotState(handles, startTime, 'Quaternion');
        x(i,1:length(pose)) = pose; 
        x(i,length(pose)*2+1:length(pose)*2+6) = wrench;
        
        % compute stiffness values
        if initialRollout
            a(i,:) = init(i,:);
        else
            a(i,:) = policy.fcn(policy, x(i,dyno(poli))', zeros(length(poli)));
        end
        
        handles.config.Mode.CartesianStiffness.Stiffness.X = a(i,1);    % TODO: change for different action dimensionality!
        handles.config.Mode.CartesianStiffness.Stiffness.Y = a(i,1);
        handles.config.Mode.CartesianStiffness.Stiffness.Z = a(i,2);
        handles.config.Mode.CartesianStiffness.Stiffness.A = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.B = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.C = 100;
        
        % send service request
        response = call(handles.client, handles.config, 'Timeout', 0.2);
        if response.Success
            
                        % New reference per point of interest:
%                         if i < round(H/3) && ~pose1Send
%                             pos = transl(Href(:,:,1));              % position
%                             quat = Quaternion(Href(:,:,1)).double;  % orientation
%                             poseCommandMsg.Pose.Position.X = pos(1);
%                             poseCommandMsg.Pose.Position.Y = pos(2);
%                             poseCommandMsg.Pose.Position.Z = pos(3);
%                             poseCommandMsg.Pose.Orientation.X = quat(1);
%                             poseCommandMsg.Pose.Orientation.Y = quat(2);
%                             poseCommandMsg.Pose.Orientation.Z = quat(3);
%                             poseCommandMsg.Pose.Orientation.W = quat(4);
%                             send(handles.pubCartesianPose,poseCommandMsg);
%                             pose1Send = true;
%                         elseif i>round(H/3) && i<round(H*2/3) && ~pose2Send
%                             pos = transl(Href(:,:,2));              % position
%                             quat = Quaternion(Href(:,:,2)).double;  % orientation
%                             poseCommandMsg.Pose.Position.X = pos(1);
%                             poseCommandMsg.Pose.Position.Y = pos(2);
%                             poseCommandMsg.Pose.Position.Z = pos(3);
%                             poseCommandMsg.Pose.Orientation.X = quat(1);
%                             poseCommandMsg.Pose.Orientation.Y = quat(2);
%                             poseCommandMsg.Pose.Orientation.Z = quat(3);
%                             poseCommandMsg.Pose.Orientation.W = quat(4);
%                             send(handles.pubCartesianPose,poseCommandMsg);
%                             pose2Send = true;
%                         elseif i>round(H*2/3) && ~pose3Send
%                             pos = transl(Href(:,:,3));              % position
%                             quat = Quaternion(Href(:,:,3)).double;  % orientation
%                             poseCommandMsg.Pose.Position.X = pos(1);
%                             poseCommandMsg.Pose.Position.Y = pos(2);
%                             poseCommandMsg.Pose.Position.Z = pos(3);
%                             poseCommandMsg.Pose.Orientation.X = quat(1);
%                             poseCommandMsg.Pose.Orientation.Y = quat(2);
%                             poseCommandMsg.Pose.Orientation.Z = quat(3);
%                             poseCommandMsg.Pose.Orientation.W = quat(4);
%                             send(handles.pubCartesianPose,poseCommandMsg);
%                             pose3Send = true;
%                         end
            
            % New Reference time step:
            posCommand = transl(Href(:,:,i));              % position Vector
            quatCommand = Quaternion(Href(:,:,i)).double;  % orientation Quaternion
            r(i,1:length([posCommand', quatCommand])) = [posCommand', quatCommand];
            
            poseCommandMsg.Pose.Position.X = posCommand(1);
            poseCommandMsg.Pose.Position.Y = posCommand(2);
            poseCommandMsg.Pose.Position.Z = posCommand(3);
            poseCommandMsg.Pose.Orientation.X = quatCommand(1);
            poseCommandMsg.Pose.Orientation.Y = quatCommand(2);
            poseCommandMsg.Pose.Orientation.Z = quatCommand(3);
            poseCommandMsg.Pose.Orientation.W = quatCommand(4);
            
            send(handles.pubCartesianPose,poseCommandMsg);
        elseif ~isempty(response.Error)
            error('SmartServo Service not reached in time: %s', response.Error);
        else
            warning('SmartServo not responding...');
        end
        
        waitfor(handles.r);
    end
catch ME
    assignin('base', 'me', ME);
    disp(ME);
    rosshutdown;
end

[time(H+1), pose, wrench] = readState(handles, startTime);
x(H+1,1:length(pose)) = pose; x(H+1,length(pose)+1:length(pose)+6) = wrench;
r(H+1,1:length([posCommand', quatCommand])) = [posCommand', quatCommand];

% move end-effector up 10 cm to relax:
latestPose   = receive(handles.subCartesianPose);
latestPose.Pose.Position.Z = latestPose.Pose.Position.Z + 0.1;
send(handles.pubCartesianPose,latestPose);

rosshutdown;

fprintf('Timing data:\n')
if handles.r.statistics.NumOverruns > 0
    warning('The rollout encountered %i timing violations',handles.r.statistics.NumOverruns);
end
disp(handles.r.statistics)

%% Data processing

% % Filter position data:
% lpFilt = designfilt('lowpassfir', ...
%     'FilterOrder', 25, ...
%     'PassbandFrequency', 50, ...
%     'StopbandFrequency', 60,...
%     'DesignMethod', 'ls', ...
%     'SampleRate', 250);
% fvtool(lpFilt) % visualize filter response
% filtPose = filter(lpFilt,x(:,1:6)); % apply filter to your data
%
% [~, cartVelocity] = gradient(filtPose,plant.dt);

% non-filtered velocity:
[~,cartVelocity] = gradient(x(:,1:length(pose)),plant.dt);
x(:,length(pose)+1:length(pose)*2) = cartVelocity;

for i=1:H
    L(i,1) = cost.fcn(cost, x(i,dyno), zeros(length(dyno)), a(i,:)', policy);    % compute rollout cost w/ energy penalty
end
L = L(1:H,1);               % Cost

y = x(2:H+1,1:nX);          % successor states  [H x nX]
x = [x(1:H,:) a(1:H,:)];    % observed  states  [H x nX+nU]

latent = x;                 % actual latent states unknown, so measured state is passed for congruency.

r = r(1:H+1,:);

end




function [time, pose, wrench] = readState(handles, startTime)
%readState Extract pose and wrench from iiwa message
%
%%
% Read latest state messages:
latestPose   = receive(handles.subCartesianPose);
latestWrench = receive(handles.subCartesianWrench);

% compute relative time since start of rollout:
absTime = double(latestPose.Header.Stamp.Sec) + double(latestPose.Header.Stamp.Nsec)/1e9;
time = absTime-startTime;

% Extract numerical pose data:
poseMsg = latestPose.Pose;
xpos = poseMsg.Position.X; ypos = poseMsg.Position.Y; zpos = poseMsg.Position.Z;
quat = poseMsg.Orientation;
%angles = quat2eul([quat.W quat.X quat.Y quat.Z],'ZYZ');
angles = [quat.W quat.X quat.Y quat.Z];
pose = double([xpos, ypos, zpos, angles]);

% extract numerical wrench data:
wrenchMsg = latestWrench.Wrench;
xforce = wrenchMsg.Force.X; yforce = wrenchMsg.Force.Y; zforce = wrenchMsg.Force.Z;
xtorque = wrenchMsg.Torque.X; ytorque = wrenchMsg.Torque.Y; ztorque = wrenchMsg.Torque.Z;
wrench = double([xforce, yforce, zforce, xtorque, ytorque, ztorque]);
end
