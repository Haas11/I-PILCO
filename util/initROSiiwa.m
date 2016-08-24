function handles = initROSiiwa(bufferSize, dt, relVelocity, overrunAction, tool_frame, varargin)
% Function to initialize the ROS network to interface with the KUKA LBR
% iiwa robot. Maximum attainable sampling frequency around 10 Hz.
% For faster sampling consider switching to C or C++ scripts.
% Relative velocity only works in joint position mode. Ignored for
% Cartesian pose commands.
%
% When providing an additional input for control mode, stiffnesses are set
% to midpoint values (2500 2500 2500 150 150 150) and a damping ratio of
% 0.7 for all Cartesian Directions.
%
% INPUTS:
%       bufferSize  - buffer size for the ROS messages  [1,inf)
%       dt          - sampling frequency for ROS rate object    (0,inf)
%       relVelocity - relative velocity of joint position motions (0,1]
%       overrunAction - What to do when overshooting sampling instance
%       ('slip' or 'drop')
% OPTIONAL INPUTS:
%       Mode        - Choose specific mode to initialize into.
%           0 = Joint Position
%           1 = Joint Impedance
%           2 = Cartesian Impedance
%
% Intializes:
%       - ROS network connecting to iiwa LBR7
%
%       - Subscribers:
%           End-effector wrench
%           Cartesian pose
%           Joint angles
%           Joint torques
%
%       - Publishers:
%           Cartesian Pose commands
%           Joint angle commands
%
%       - Services:
%           ConfigureSmartServo
%               'CartesianImpedance'
%               'JointImpedance'
%               'JointPosition'
%
%       - Rate Object
%
%
% Code by Victor van Spaandonk
% July 2016

%% Code
persistent CONNECTED
if isempty(CONNECTED)
    CONNECTED = false;
end
    
if ~CONNECTED
    validAnswer = false;
    while ~validAnswer
        reply = input('\nROSSmartServo Application Running on Cabinet and roscore active? [y/n]...','s');
        
        if strcmpi(reply,'y')
            validAnswer = true;
        elseif strcmpi(reply,'n')
            warning('Run the application ROSSmartServo on the KUKA Cabinet and optionally launch a roscore with additional nodes such as RVIZ.)');
        else
            warning('please hit "y" or "n"');
        end
    end
end

fprintf('\n\n --- Initializing ROS w/ LBR iiwa --- \n\n');

handles.ipaddress = 'http://192.170.10.5:11311';
rosshutdown;            % [bugfix] shutdown first for debug
rosinit(handles.ipaddress);     % start MATLAB ROS Node

handles.frequency = 1/dt;             % [Hz]

% subscribers:
handles.subCartesianPose   = rossubscriber('/iiwa/state/CartesianPose','BufferSize', bufferSize);
handles.subCartesianWrench = rossubscriber('/iiwa/state/CartesianWrench','BufferSize', bufferSize);
handles.subJointPosition   = rossubscriber('/iiwa/state/JointPosition','BufferSize', bufferSize);
handles.subJointTorque     = rossubscriber('/iiwa/state/JointTorque','BufferSize', bufferSize);

% publishers
handles.pubCartesianPose = rospublisher('/iiwa/command/CartesianPose', 'IsLatching', false);   % Cartesian Pose
handles.poseCommandMsg   = rosmessage(handles.pubCartesianPose);    % initialize pose command message
handles.poseCommandMsg.Header.FrameId = tool_frame;         % tool frame

handles.pubJointPosition = rospublisher('/iiwa/command/JointPosition', 'IsLatching', false);    % Joint Position
handles.jointCommandMsg = rosmessage(handles.pubJointPosition);     % joint position command message
handles.jointCommandMsg.Header.FrameId = 'Robot';                   % (don't change)

% services
handles.client = rossvcclient('/iiwa/configuration/configureSmartServo');
handles.config = rosmessage(handles.client);
handles.config.Mode.RelativeVelocity = relVelocity;

% ROS Rate Timer:
handles.r = robotics.Rate(handles.frequency);
handles.r.OverrunAction = overrunAction;               % slip = don't skip step when violating sampling constraint, drop = skip one step

fprintf('\n --- KUKA iiwa succesfully initialized --- \n');
CONNECTED = true;

handles.config.Mode.CartesianStiffness.Stiffness.X = 2500;      % [0 5000] N/m
handles.config.Mode.CartesianStiffness.Stiffness.Y = 2500;
handles.config.Mode.CartesianStiffness.Stiffness.Z = 2500;
handles.config.Mode.CartesianStiffness.Stiffness.A = 150;       % [0 300] Nm/rad
handles.config.Mode.CartesianStiffness.Stiffness.B = 150;
handles.config.Mode.CartesianStiffness.Stiffness.C = 150;

handles.config.Mode.CartesianDamping.Damping.X = 0.7;
handles.config.Mode.CartesianDamping.Damping.Y = 0.7;
handles.config.Mode.CartesianDamping.Damping.Z = 1;
handles.config.Mode.CartesianDamping.Damping.A = 0.7;
handles.config.Mode.CartesianDamping.Damping.B = 0.7;
handles.config.Mode.CartesianDamping.Damping.C = 0.7;

handles.config.Mode.NullspaceStiffness = 150;
handles.config.Mode.NullspaceDamping = 0.7;

if nargin > 5
    mode = varargin{1};
    if mode==0
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.JOINTIMPEDANCE;
        response = call(handles.client, handles.config, 'Timeout', 5);
   
    elseif mode==1
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.CARTESIANIMPEDANCE;                
        response = call(handles.client, handles.config, 'Timeout', 5);      % Set Control Mode
                
    elseif mode==2
        handles.config.Mode.Mode = robotics.ros.custom.msggen.iiwa_msgs.SmartServoMode.CONSTANTFORCE;
        response = call(handles.client, handles.config, 'Timeout', 5);
    else
        error('Invalid ControlMode. 0=JOINT IMPEDANCE, 1=CARTESIAN IMPEDANCE, 2=CONSTANT FORCE');
    end
    
    if response.Success
        fprintf('\n --- Successfully initialized desired control mode --- \n');
    else
        error('Failed to start desired mode: %s', response.Error);
    end
end

%% Display Result
fprintf('\n');
disp(handles);
