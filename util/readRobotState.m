function [time, states_out] = readRobotState(handles, startTime, desiredStates, varargin)
%readState Extract desired state information from KUKA LBR7 iiwa
%communicating through matlabROS.
% 'desiredStates' is a 1x4 vector containing either ones or zeros for the
% states that should be queried. The output contains a relevant time
% instance based on the start of a rollout and a structure 'states_out'. 
% This structure only contains the fields of the states that were queried.
%
%
% Inputs are :
%   - handles passed by initRosiiwa.m
%
%   - startTime = time just before control loop was started to make the time
%
%   - states = [jointPosition, jointTorque, CartesianPose, CartesianWrench]
%
% associated to the states relevant to it (instead of january 1st 1970).
%   - (optional)
%     orientationType : 'QUAT'/'XYZ'/'ZYZ'
%
% outputs:
%       - time              [1 x 1]
%       - state  (structure)
%           .JointPosition
%           .JointTorque
%           .CartesianPose
%           .CartesianWrench
%       - pose
%           Euler:          [6 x 1]
%           Quaternion:     [7 x 1]
%       - wrench:           [6 x 1]
%
% code by Victor van Spaandonk.
% July 2016

%% Code
types = {'JointPosition','JointTorque','CartesianPose','CartesianWrench'};
time = [];

if desiredStates(1)
    latestJointAngles = receive(handles.subJointPosition);
    
    q = zeros(1,7);
    for i=1:7
        joint = strcat('A',num2str(i));
        q(i) = latestJointAngles.Position.(joint);
    end
    
    % compute relative time since start of rollout:
    absTime = double(latestJointAngles.Header.Stamp.Sec) + double(latestJointAngles.Header.Stamp.Nsec)/1e9;
    time = absTime-startTime;
    
    states_out.(types{1}) = q;
end

if desiredStates(2)    
    latestJointTorques = receive(handles.handles.subJointTorque);
    
    tau = zeros(7,1);
    for i=1:7
        joint = strcat('A',num2str(i));
        tau(i) = latestJointTorques.Torque.(joint);
    end
    states_out.(types{2}) = tau;
    
    if isempty(time)
        % compute relative time since start of rollout:
        absTime = double(latestJointTorques.Header.Stamp.Sec) + double(latestJointTorques.Header.Stamp.Nsec)/1e9;
        time = absTime-startTime;
    end
    
end

% Cartesian Pose
if desiredStates(3)
    latestPose = receive(handles.subCartesianPose);
    
    % Extract numerical pose data:
    poseMsg = latestPose.Pose;
    xpos = poseMsg.Position.X; ypos = poseMsg.Position.Y; zpos = poseMsg.Position.Z;
    quat = poseMsg.Orientation;
    angles = [quat.W quat.X quat.Y quat.Z];
    if nargin > 3
        if strcmpi(varargin{1},'ZYZ')
            angles = quat2eul([quat.W quat.X quat.Y quat.Z],'ZYZ');
        elseif strcmpi(varargin{1},'ZYX')
            angles = quat2eul([quat.W quat.X quat.Y quat.Z],'ZYX');
        elseif strcmpi(varargin{1},'QUAT')
            % Default
        else
            warning('Invalid angle type. Defaulting to quaternions.');
        end
    end
    pose = double([xpos, ypos, zpos, angles]);
    states_out.(types{3}) = pose;
    
    if isempty(time)
        % compute relative time since start of rollout:
        absTime = double(latestPose.Header.Stamp.Sec) + double(latestPose.Header.Stamp.Nsec)/1e9;
        time = absTime-startTime;
    end
    
end

% Cartesian Wrench
if desiredStates(4)
    latestWrench = receive(handles.subCartesianWrench);
    
    % compute relative time since start of rollout:
    absTime = double(latestWrench.Header.Stamp.Sec) + double(latestWrench.Header.Stamp.Nsec)/1e9;
    time = absTime-startTime;        
    
    % extract numerical wrench data:
    wrenchMsg = latestWrench.Wrench;
    wrench = double([wrenchMsg.Force.X, wrenchMsg.Force.Y, wrenchMsg.Force.Z, wrenchMsg.Torque.X, wrenchMsg.Torque.Y, wrenchMsg.Torque.Z]);
    
    if isempty(time)
        % compute relative time since start of rollout:
        absTime = double(latestWrench.Header.Stamp.Sec) + double(latestJointTorques.Header.Stamp.Nsec)/1e9;
        time = absTime-startTime;
    end
    
    states_out.(types{4}) = wrench;
end
end
