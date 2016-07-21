clc;
clear;

ipaddress = 'http://192.170.10.5:11311';
rosshutdown;
rosinit(ipaddress);

genKUKATrajectory;

handles.N = floor(T/dt);

handles.frequency = 1/dt;             % [Hz]

handles.initJoints = true;     % init in joint position mode.

% subscribers:
handles.subCartesianPose = rossubscriber('/iiwa/state/CartesianPose','BufferSize', 25);
handles.subCartesianWrench = rossubscriber('/iiwa/state/CartesianWrench','BufferSize', 25);

% publishers
handles.pubCartesianPose  = rospublisher('/iiwa/command/CartesianPose', 'IsLatching', false);
handles.pubJointPosition = rospublisher('/iiwa/command/JointPosition', 'IsLatching', false);

% service
handles.client = rossvcclient('/iiwa/configuration/configureSmartServo');
handles.config = rosmessage(handles.client);
handles.config.Mode.RelativeVelocity = 0.1;

handles.r = robotics.Rate(handles.frequency);
handles.r.OverrunAction = 'drop';
% handles.r = robotics.ros.Rate('/iiwa/robot_state_publisher', frequency);
% handles.r = rosrate(10, 'OverrunAction', 'drop');                     % loop rate
try
    [data, handles] = iiwaRolloutStart(handles, Href);          % perform rollout:
    %[data, handles] = iiwaRolloutStart1(handles, Href1);          % perform rollout:
    
    rosshutdown;
    
    fprintf('Time statistics:');
    handles.r
    handles.r.statistics
    
    figure
    %stairs(data.pose(:,1),data.pose(:,2:end));
    stairs(t(1:size(data.pose,1)),data.pose(:,2:end));
    hold on
    stairs(t(1:size(data.pose,1)),transl(Href(:,:,1:size(data.pose,1))),'k:');
    grid on;
    title('Pose')
    legend('p_x','p_y','p_z','\phi','\nu','\theta','w','ref');
    xlabel('Time [s]')
    
    figure
    plot(handles.r.statistics.Periods,'r*--');
    grid on;
    xlabel('Time step'); ylabel('Sample time [s]');
    title('Loop timing');
    
catch ME
    ME
    rosshutdown;                        % shutdown the ros node
end
%%