function handles = initROSiiwa(bufferSize, dt, relVelocity, overrunAction, varargin)
% Function to initialize the ROS network to interface with the KUKA LBR
% iiwa robot. Maximum attainable sampling frequency around 10 Hz.
% For faster sampling consider switching to C or C++ scripts.
% Relative velocity only works in joint position mode. Ignored for
% Cartesian pose commands.
%
% INPUTS:
%       bufferSize  - buffer size for the ROS messages  [1,inf)
%       dt          - sampling frequency for ROS rate object    (0,inf)
%       relVelocity - relative velocity of joint position motions (0,1]
%       overrunAction - What to do when overshooting sampling instance
%       ('slip' or 'drop')
%
% Intializes:
%       - ROS network connecting to iiwa LBR7
%
%       - Subscribers:
%           End-effector wrench
%           Cartesian pose
%           Joint angles
%
%       - Publishers:
%           Cartesian Pose commands
%           Joint angle commands
%
%       - Services:
%           ConfigureSmartServo
%               Cartesian Impedance = 0
%               Joint Impedance = 1
%               Joint Position  = 2
%
%       - Rate Object
%
%
% Code by Victor van Spaandonk
% July 2016

%% Code 

fprintf('\n\n --- Initializing ROS w/ LBR iiwa --- \n\n');

handles.ipaddress = 'http://192.170.10.5:11311';
rosshutdown;            % [bugfix] shutdown first for debug
rosinit(handles.ipaddress);     % start MATLAB ROS Node

handles.frequency = 1/dt;             % [Hz]

% subscribers:
handles.subCartesianPose = rossubscriber('/iiwa/state/CartesianPose','BufferSize', bufferSize);
handles.subCartesianWrench = rossubscriber('/iiwa/state/CartesianWrench','BufferSize', bufferSize);

% publishers
handles.pubCartesianPose  = rospublisher('/iiwa/command/CartesianPose', 'IsLatching', false);   % Cartesian Pose
poseCommandMsg = rosmessage(handles.pubCartesianPose);  % initialize pose command message
poseCommandMsg.Header.FrameId = 'peg_link_ee_kuka';     % tool frame

handles.pubJointPosition = rospublisher('/iiwa/command/JointPosition', 'IsLatching', false);    % Joint Position
jointCommandMsg = rosmessage(handles.pubJointPosition);     % joint position command message
jointCommandMsg.Header.FrameId = 'Robot';                   % (don't change)

% services
handles.client = rossvcclient('/iiwa/configuration/configureSmartServo');
handles.config = rosmessage(handles.client);
handles.config.Mode.RelativeVelocity = relVelocity;

% ROS Rate Timer:
handles.r = robotics.Rate(handles.frequency);
handles.r.OverrunAction = overrunAction;               % slip = don't skip step when violating sampling constraint, drop = skip one step

fprintf('\nKUKA iiwa succesfully initialized\n');

if nargin > 5
    mode = varargin{1};
    if strcmpi(mode,'CartesianImpedance')
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.CARTESIANIMPEDANCE;
        response = call(handles.client, handles.config, 'Timeout', 5);
    elseif strcmpi(mode,'JointImpedance')
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.JOINTIMPEDANCE;
        response = call(handles.client, handles.config, 'Timeout', 5);

    elseif strcmpi(mode,'JointPosition')
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.JOINTPOSITION;
        response = call(handles.client, handles.config, 'Timeout', 5);
    end
    
    if response.Success
        fprintf('\nSuccessfully switched control mode.\n');
    else
        error('Failed to start desired mode: %s', response.Error);
    end
end
