function [time, data, handles] = iiwaRolloutStart1(handles, Href)
%
%
%
%
% Code by Victor van Spaandonk
% July 2016

%% Initialize
handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.CARTESIANIMPEDANCE;

handles.config.Mode.CartesianStiffness.Stiffness.X = 1500;
handles.config.Mode.CartesianStiffness.Stiffness.Y = 1500;
handles.config.Mode.CartesianStiffness.Stiffness.Z = 1500;
handles.config.Mode.CartesianStiffness.Stiffness.A = 100;
handles.config.Mode.CartesianStiffness.Stiffness.B = 100;
handles.config.Mode.CartesianStiffness.Stiffness.C = 100;

handles.config.Mode.CartesianDamping.Damping.X = 0.7;
handles.config.Mode.CartesianDamping.Damping.Y = 0.7;
handles.config.Mode.CartesianDamping.Damping.Z = 0.7;
handles.config.Mode.CartesianDamping.Damping.A = 0.7;
handles.config.Mode.CartesianDamping.Damping.B = 0.7;
handles.config.Mode.CartesianDamping.Damping.C = 0.7;

pose = cell(handles.N,1);
wrench = cell(handles.N,1);

pose1Send = false;  pose2Send = false;  pose3Send = false;

try
    response = call(handles.client, handles.config, 'Timeout', 5);
    pause(5);
    if response.Success
        fprintf('\n\n======= CARTESIAN IMPEDANCE MODE ACTIVE ======== \n\n');
    end
    
    poseCommandMsg = rosmessage(handles.pubCartesianPose);
    poseCommandMsg.Header.FrameId = 'peg_link_ee_kuka';
    
    pos = transl(Href(:,:,1));              % initial position
    quat = Quaternion(Href(:,:,1)).double;  % initial orientation
    xe_0 = [pos; quat'];                    % initial state
    
    % TODO: consider doing in joint position mode
    
    poseCommandMsg.Pose.Position.X = xe_0(1);
    poseCommandMsg.Pose.Position.Y = xe_0(2);
    poseCommandMsg.Pose.Position.Z = xe_0(3);
    poseCommandMsg.Pose.Orientation.X = xe_0(4);
    poseCommandMsg.Pose.Orientation.Y = xe_0(5);
    poseCommandMsg.Pose.Orientation.Z = xe_0(6);
    poseCommandMsg.Pose.Orientation.W = xe_0(7);
    
    send(handles.pubCartesianPose, poseCommandMsg);
    
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
    t = rostime('now');
    startTime = double(t.Sec) + double(t.Nsec)/1e9;
    
    % loop
    for i=1:handles.N
        
        % read state information
        [time(i), pose, wrench] = readState(handles, startTime);
        
        % compute stiffness values
        handles.config.Mode.Mode = -1;      % indicate that control mode is the same    (saves a lot of time!)
        handles.config.Mode.CartesianStiffness.Stiffness.X = 1500;
        handles.config.Mode.CartesianStiffness.Stiffness.Y = 1500;
        handles.config.Mode.CartesianStiffness.Stiffness.Z = 1500;
        handles.config.Mode.CartesianStiffness.Stiffness.A = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.B = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.C = 100;
        
        % send service request
        response = call(handles.client, handles.config, 'Timeout', 0.2);
        if response.Success
            if i < round(handles.N/3) && ~pose1Send
                pos = transl(Href(:,:,2));              % position
                quat = Quaternion(Href(:,:,2)).double;  % orientation
                poseCommandMsg.Pose.Position.X = pos(1);
                poseCommandMsg.Pose.Position.Y = pos(2);
                poseCommandMsg.Pose.Position.Z = pos(3);
                poseCommandMsg.Pose.Orientation.X = quat(1);
                poseCommandMsg.Pose.Orientation.Y = quat(2);
                poseCommandMsg.Pose.Orientation.Z = quat(3);
                poseCommandMsg.Pose.Orientation.W = quat(4);                
                send(handles.pubCartesianPose,poseCommandMsg);
                pose1Send = true;
            elseif i>round(handles.N/3) && i<round(handles.N*2/3) && ~pose2Send
                pos = transl(Href(:,:,3));              % position
                quat = Quaternion(Href(:,:,3)).double;  % orientation
                poseCommandMsg.Pose.Position.X = pos(1);
                poseCommandMsg.Pose.Position.Y = pos(2);
                poseCommandMsg.Pose.Position.Z = pos(3);
                poseCommandMsg.Pose.Orientation.X = quat(1);
                poseCommandMsg.Pose.Orientation.Y = quat(2);
                poseCommandMsg.Pose.Orientation.Z = quat(3);
                poseCommandMsg.Pose.Orientation.W = quat(4);                
                send(handles.pubCartesianPose,poseCommandMsg);
                pose2Send = true;
            elseif i>round(handles.N*2/3) && ~pose3Send
                pos = transl(Href(:,:,4));              % position
                quat = Quaternion(Href(:,:,4)).double;  % orientation
                poseCommandMsg.Pose.Position.X = pos(1);
                poseCommandMsg.Pose.Position.Y = pos(2);
                poseCommandMsg.Pose.Position.Z = pos(3);
                poseCommandMsg.Pose.Orientation.X = quat(1);
                poseCommandMsg.Pose.Orientation.Y = quat(2);
                poseCommandMsg.Pose.Orientation.Z = quat(3);
                poseCommandMsg.Pose.Orientation.W = quat(4);                
                send(handles.pubCartesianPose,poseCommandMsg);
                pose3Send = true;
            end
 
        elseif ~isempty(response.Error)
            error('SmartServo Service not reached in time %s', response.Error);
        else
            warning('SmartServo not responding...');
        end
        
        waitfor(handles.r);
    end
catch ME
    ME
    rosshutdown;
end

data.pose = pose;
data.wrench = wrench;

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