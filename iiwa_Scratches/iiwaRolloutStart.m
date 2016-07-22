function [data, handles] = iiwaRolloutStart(handles, Href)
%
%
%
%
% Code by Victor van Spaandonk
% July 2016

%% Initialize
pose = zeros(handles.N,8);
wrench = cell(handles.N,1);

    poseCommandMsg = rosmessage(handles.pubCartesianPose);
    poseCommandMsg.Header.FrameId = 'peg_link_ee_kuka';
    
    jointCommandMsg = rosmessage(handles.pubJointPosition);
    jointCommandMsg.Header.FrameId = 'Robot';
        
    pos = transl(Href(:,:,1));              % initial position
    quat = Quaternion(Href(:,:,1)).double;  % initial orientation
    xe_0 = [pos; quat'];                    % initial state
    
    % Initialize the pose of the robot:
    q0 = [-0.00122388906311;
          0.259632468224;
          0.00148364703637;
          -1.42674195766;
          -0.000345433305483;
          1.455529809;
          0.000291096803267];
    
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
    % Switch to Cartesian Impedance mode:
    response = call(handles.client, handles.config, 'Timeout', 5);
    pause(2);
    if response.Success
        fprintf('\n\n======= CARTESIAN IMPEDANCE MODE ACTIVE ======== \n\n');
    else
        error('Failed to start CARTESIAN IMPEDANCE MODE: %s', response.Error);
    end
    
    fprintf('\n\n======= MOVING TO START STATE ======== \n\n');
    pause(2);
    if handles.initJoints
        jointCommandMsg.Position.A1 = q0(1);
        jointCommandMsg.Position.A2 = q0(2);
        jointCommandMsg.Position.A3 = q0(3);
        jointCommandMsg.Position.A4 = q0(4);
        jointCommandMsg.Position.A5 = q0(5);
        jointCommandMsg.Position.A6 = q0(6);
        jointCommandMsg.Position.A7 = q0(7);       
        send(handles.pubJointPosition, jointCommandMsg);
    else
        poseCommandMsg.Pose.Position.X = xe_0(1);
        poseCommandMsg.Pose.Position.Y = xe_0(2);
        poseCommandMsg.Pose.Position.Z = xe_0(3);
        poseCommandMsg.Pose.Orientation.X = xe_0(4);
        poseCommandMsg.Pose.Orientation.Y = xe_0(5);
        poseCommandMsg.Pose.Orientation.Z = xe_0(6);
        poseCommandMsg.Pose.Orientation.W = xe_0(7);
        send(handles.pubCartesianPose, poseCommandMsg);
    end
    
    %% Perform Rollout
    validAnswer = false;
    while ~validAnswer
        reply = input('Ready to start rollout? [Y/N]...','s');
        
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
    for i=1:handles.N        
        % read state information
        pose(i,:) = readPose(handles.subCartesianPose.LatestMessage, startTime);
        wrench{i} = readWrench(handles.subCartesianWrench.LatestMessage);
        
%         latestPose   = receive(handles.subCartesianPose);
%         latestWrench = receive(handles.subCartesianWrench);
%         
%         pose(i,:) = readPose(latestPose, startTime);
%         wrench{i} = readWrench(latestWrench);
        
        %disp(pose); disp(wrench);
        
        % compute stiffness values
        handles.config.Mode.Mode = -1;      % indicate that control mode is the same
        handles.config.Mode.CartesianStiffness.Stiffness.X = 1500;
        handles.config.Mode.CartesianStiffness.Stiffness.Y = 1000;
        handles.config.Mode.CartesianStiffness.Stiffness.Z = 50;
        handles.config.Mode.CartesianStiffness.Stiffness.A = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.B = 100;
        handles.config.Mode.CartesianStiffness.Stiffness.C = 100;
                               
        % send service request
        response = call(handles.client, handles.config, 'Timeout', 0.2);
        if response.Success
            
            pos = transl(Href(:,:,i));              % position
            quat = Quaternion(Href(:,:,i)).double;  % orientation
            
            poseCommandMsg.Pose.Position.X = pos(1);
            poseCommandMsg.Pose.Position.Y = pos(2);
            poseCommandMsg.Pose.Position.Z = pos(3);
            poseCommandMsg.Pose.Orientation.X = quat(1);
            poseCommandMsg.Pose.Orientation.Y = quat(2);
            poseCommandMsg.Pose.Orientation.Z = quat(3);
            poseCommandMsg.Pose.Orientation.W = quat(4);
            
            send(handles.pubCartesianPose,poseCommandMsg);
        elseif ~isempty(response.Error)
            error('SmartServo Service not reached in time: %s', response.Error);
        else
            warning('SmartServo not responding...');
        end
        
        % publish reference command
        %if handles.r.LastPeriod > 1/handles.frequency
        %    warning('Sampling time constraint was violated in iteration %1i by %4.5f [s].', i, handles.r.LastPeriod - 1/handles.frequency)
        %end
        waitfor(handles.r);
    end
catch ME
    ME %#ok<NOPRT>
    rosshutdown;
end

pose(:,1) = pose(:,1);
data.pose = pose;
data.wrench = wrench;

end




function pose = readPose(poseMsg_in, startTime)
%readPose Extract pose from iiwa message
pose = zeros(1,7);

time = double(poseMsg_in.Header.Stamp.Sec) + double(poseMsg_in.Header.Stamp.Nsec)/1e9;

poseMsg = poseMsg_in.Pose;

xpos = poseMsg.Position.X;
ypos = poseMsg.Position.Y;
zpos = poseMsg.Position.Z;

quat = poseMsg.Orientation;
%angles = quat2eul([quat.W quat.X quat.Y quat.Z],'ZYZ');
angles = [quat.W quat.X quat.Y quat.Z];

pose(1,2:8) = [xpos, ypos, zpos, angles];
pose(1,1) = time-startTime;
pose = double(pose);
end

function wrench = readWrench(wrenchMsg_in)
%readPose Extract wrench from iiwa message

wrenchMsg = wrenchMsg_in.Wrench;

xforce = wrenchMsg.Force.X;
yforce = wrenchMsg.Force.Y;
zforce = wrenchMsg.Force.Z;

xtorque = wrenchMsg.Torque.X;
ytorque = wrenchMsg.Torque.Y;
ztorque = wrenchMsg.Torque.Z;

wrench.data = [xforce, yforce, zforce, xtorque, ytorque, ztorque];
wrench.time = wrenchMsg_in.Header.Stamp.Sec;
end