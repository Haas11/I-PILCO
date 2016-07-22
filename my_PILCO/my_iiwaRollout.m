function [x, y, L, latent] = my_iiwaRollout(policy, plant, H, Href, init)
%
% Inputs:
%   - handles   handles for ROS communication
%   - policy    controller
%   - plant     plant structure
%   - H         Number of timesteps
%   - Href      Reference
%
% Outputs:
%   - x         observed states
%   - y         successor states
%   - L         Cost
%
%
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

handles = initRosiiwa(10, plant.dt, 0.25, 'slip');

% state vector:
x = zeros(handles.N,18);    % state vector [x dx W]

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
    
    fprintf('\n...\n======= MOVING TO START STATE ======== \n...\n');
    pause(2);
    
    % Initialize the pose of the robot:
    q0 = [-0.00122388906311;
        0.259632468224;
        0.00148364703637;
        -1.42674195766;
        -0.000345433305483;
        1.455529809;
        0.000291096803267];
    
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
        fprintf('\n\n======= CARTESIAN IMPEDANCE MODE ACTIVE ======== \n\n');
    else
        error('Failed to start CARTESIAN IMPEDANCE MODE: %s', response.Error);
    end
    
    %% Perform Rollout
    validAnswer = false;
    while ~validAnswer
        reply = input('Ready to start rollout? [Y/N] \n','s');
        
        if strcmpi(reply,'y')
            validAnswer = true;
        elseif strcmpi(reply,'n')
            error('Rollout Aborted by User');
        else
            warning('please hit "y" or "n"');
        end
    end
    
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
        [time(i), pose, wrench] = readState(handles, startTime);
        x(i,1:6) = pose; x(i,13:18) = wrench;
        
        % compute stiffness values
        if initialRollout
            a = init(i,:);
        else
            a = policy.fcn(policy, x(i,dyno(poli))', zeros(length(poli)));
        end
        
        handles.config.Mode.CartesianStiffness.Stiffness.X = a(1);    % TODO: change for different action dimensionality!
        handles.config.Mode.CartesianStiffness.Stiffness.Y = a(1);
        handles.config.Mode.CartesianStiffness.Stiffness.Z = a(2);
        handles.config.Mode.CartesianStiffness.Stiffness.A = 150;
        handles.config.Mode.CartesianStiffness.Stiffness.B = 150;
        handles.config.Mode.CartesianStiffness.Stiffness.C = 150;
        handles.config.Mode.Mode = -1;      % indicate that control mode is the same (time saving)
        
        % send service request
        response = call(handles.client, handles.config, 'Timeout', 0.1);
        if response.Success
            posCommand = transl(Href(:,:,i));              % position Vector
            quatCommand = Quaternion(Href(:,:,i)).double;  % orientation Quaternion
            
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
        
        %if handles.r.LastPeriod > 1/handles.frequency
        %    warning('Sampling time constraint was violated in iteration %1i by %4.5f [s].', i, handles.r.LastPeriod - 1/handles.frequency)
        %end
        waitfor(handles.r);
    end
catch ME
    ME %#ok<NOPRT>
    rosshutdown;
end


%% Data processing

% Filter position data:
lpFilt = designfilt('lowpassfir', ...
    'FilterOrder', 25, ...
    'PassbandFrequency', 50, ...
    'StopbandFrequency', 60,...
    'DesignMethod', 'ls', ...
    'SampleRate', 250);
fvtool(lpFilt) % visualize filter response
filtPose = filter(lpFilt,x(:,1:6)); % apply filter to your data

[~, cartVelocity] = gradient(filtPose,plant.dt);
x(:,7:12) = cartVelocity;

% dimensionality check on outputs:
y = x(2:H+1,1:nX);          % successor states  [H x nX]
x = [x(1:H,:) a(1:H,:)];    % observed  states  [H x nX+2*nU]

L = L(1:H,1);               % Cost

latent = x;                 % actual latent states unknown, so measured state is passed for congruency.

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