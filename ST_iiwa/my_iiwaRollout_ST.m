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
global TIME_INPUT
if nargin > 5
    initialRollout = true;  % precomputed stiffness trajectory
else
    initialRollout = false;
end
if TIME_INPUT
    dyno = [plant.dyno, plant.indices(end)];
else
    dyno = plant.dyno;
end
poli = plant.poli;
damping = ones(1,6)*0.7;

nX = 7+7+6+1;                 % number of states = xe + dxe + F +t
nU = length(policy.maxU);   % number of actions
x = zeros(H+1,nX);
a = zeros(H,  nU);
L = zeros(H,  1);
r = zeros(H,  nX);

prevPosRef = transl(Href(:,:,1)); posRef = prevPosRef;

%% Initialize the ROS infrastructure:
validRollout = false;
while ~validRollout
    bufferSize = 1;                     % number of messages to ho
    relativeVelocity = 0.3;             % [%] ratio of full capacity
    mode = 1;                           % cartesian impedance
    tool_frame = 'peg_link_ee_kuka';    % frame used for pose commands
    handles = initROSiiwa(bufferSize, plant.dt, relativeVelocity, 'slip', tool_frame, mode);
    
    poseCommandMsg = handles.poseCommandMsg;
    jointCommandMsg = handles.jointCommandMsg;
    
    %%
    try
        fprintf('\n --- Moving to start state --- \n');
        
        [~, state] = readRobotState(handles, 0, [1 0 0 0]);
        
        for i=1:7
            joint = strcat('A',num2str(i));
            jointCommandMsg.Position.(joint) = state.JointPosition(i);
        end
        send(handles.pubJointPosition, jointCommandMsg);
        
        %     q0_up = [0; 0.21; 0; -1.75; 0; 1.19; 0];              % ~[0.45 0 0.3]
        q0_down = [0; 28; 0; -108; 0; 40; 0]./360.*2.*pi;       % ~[0.45 0 0.095]
        
        jointCommandMsg.Position.A1 = q0_down(1);
        jointCommandMsg.Position.A2 = q0_down(2);
        jointCommandMsg.Position.A3 = q0_down(3);
        jointCommandMsg.Position.A4 = q0_down(4);
        jointCommandMsg.Position.A5 = q0_down(5);
        jointCommandMsg.Position.A6 = q0_down(6);
        jointCommandMsg.Position.A7 = q0_down(7);
        send(handles.pubJointPosition, jointCommandMsg);
        pause(4);
        
        pos = transl(Href(:,:,1));  quat = Quaternion(Href(:,:,1)).double;
        xe_0 = [pos; quat'];
        poseCommandMsg.Pose.Position.X = xe_0(1);
        poseCommandMsg.Pose.Position.Y = xe_0(2);
        poseCommandMsg.Pose.Position.Z = 0.095;
        poseCommandMsg.Pose.Orientation.X = xe_0(4);
        poseCommandMsg.Pose.Orientation.Y = xe_0(5);
        poseCommandMsg.Pose.Orientation.Z = xe_0(6);
        poseCommandMsg.Pose.Orientation.W = xe_0(7);
        send(handles.pubCartesianPose, poseCommandMsg);
        pause(1);
        poseCommandMsg.Pose.Position.Z = xe_0(3);
        send(handles.pubCartesianPose, poseCommandMsg);
        
        %% Perform Rollout        
        fprintf('Starting rollout in...\n');
        pause(1);
        fprintf('\n5 \n'); beep; pause(1); fprintf('\n4 \n'); beep; pause(1); fprintf('\n3 \n');
        beep; pause(1); fprintf('\n2 \n'); beep; pause(1); fprintf('\n1 \n'); beep; pause(1);
        
        reset(handles.r);
        
        % loop
        t = rostime('now');
        startTime = double(t.Sec) + double(t.Nsec)/1e9;
        time = zeros(H+1,1);
        
        for i=1:H            
            % read state information
            [time(i), state] = readRobotState(handles, startTime, [0 0 1 1], 'QUAT');
            x(i,1:length(state.CartesianPose)) = state.CartesianPose;
            x(i,length(state.CartesianPose)*2+1:length(state.CartesianPose)*2+6) = state.CartesianWrench;
            x(i,end) = i*plant.dt+randn*1e-5;   % time slightly disturbed for numerical reasons in GP
            
            % compute stiffness values
            if initialRollout
                a(i,:) = init(i,:);
            else
                a(i,:) = policy.fcn(policy, x(i,dyno(poli))', zeros(length(poli)));
            end            
            stiffness(1:policy.impIdx) = a(i,policy.impIdx);
            [handles, response] = setRobotStiffDamp(handles, 1, stiffness, damping, 0.2, 200, 1);
                        
            % send pose command
            if response.Success               
                if isfield(policy,'refIdx') && ~isempty(policy.refIdx)
                    deltaRef = a(i,policy.refIdx)';                  % change in position reference  (X & Y)
                    posRef(1:2) = prevPosRef(1:2) + deltaRef;       % new absolute position reference   (X & Y)
                    prevPosRef = posRef;                                        % bookkeepping
                    posCommand = posRef; 
                    quatCommand = [0 1 0 0];
                else
                    posCommand = transl(Href(:,:,i));              % position Vector
                    quatCommand = [0 1 0 0];
                end
                r(i,1:length([posCommand', quatCommand])) = [posCommand', quatCommand];                
                poseCommandMsg.Pose.Position.X = posCommand(1);
                poseCommandMsg.Pose.Position.Y = posCommand(2);
                poseCommandMsg.Pose.Position.Z = posCommand(3);
                poseCommandMsg.Pose.Orientation.X = quatCommand(1);
                poseCommandMsg.Pose.Orientation.Y = quatCommand(2);
                poseCommandMsg.Pose.Orientation.Z = quatCommand(3);
                poseCommandMsg.Pose.Orientation.W = quatCommand(4);
                send(handles.pubCartesianPose,poseCommandMsg);              % send service request                
                
            elseif ~isempty(response.Error)
                error('SmartServo Service not reached in time: %s', response.Error);
            else
                warning('SmartServo not responding...');
            end
            
            waitfor(handles.r);
        end
    catch me
        assignin('base', 'me', me);
        disp(me);
        rosshutdown;
    end
    
    [time(H+1), state] = readRobotState(handles, startTime, [0 0 1 1], 'QUAT');
    x(H+1,1:length(state.CartesianPose)) = state.CartesianPose;
    x(H+1,length(state.CartesianPose)*2+1:length(state.CartesianPose)*2+6) = state.CartesianWrench;
    r(H+1,1:length([posCommand', quatCommand])) = [posCommand', quatCommand];
        
    response.Success = false;
    while ~response.Success
        fprintf('\n --- Moving away from end-pose --- \n');
        [handles, response] = setRobotStiffDamp(handles, 1, [500 500 500 150 150 150], damping, 1, 200, 1);
        
        if response.Success
            latestPose   = receive(handles.subCartesianPose);
            latestPose.Pose.Position.X = latestPose.Pose.Position.X - 0.2;
            latestPose.Pose.Position.Z = latestPose.Pose.Position.Z + 0.1;
            send(handles.pubCartesianPose,latestPose);
            pause(1);
        end
    end
    
    rosshutdown;
    
    fprintf('Timing data:\n')    
    disp(handles.r.statistics)    
    if handles.r.statistics.NumOverruns > 0
        warning('The rollout encountered %i timing violations',handles.r.statistics.NumOverruns);        
        if handles.r.statistics.StandardDeviation > 0.1
            beep; beep; beep;
            validAnswer = false;
            while ~validAnswer
                reply = input('\nAccept roll-out despite timing violations? [y/n]...','s');
                if strcmpi(reply,'y')
                    validAnswer = true;
                    validRollout = true;
                elseif strcmpi(reply,'n')
                    %restart rollout
                    validAnswer = true;
                else
                    warning('please hit "y" or "n"');
                end
            end
        else
            validRollout = true;
        end
    else
        validRollout = true;
    end
end

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
[~,cartVelocity] = gradient(x(:,1:length(state.CartesianPose)), plant.dt);
x(:,length(state.CartesianPose)+1:length(state.CartesianPose)*2) = cartVelocity;

% Reference Velocities:
[~,refVelocity] = gradient(r(:,1:length([posCommand', quatCommand])), plant.dt);
r(:,length([posCommand', quatCommand])+1:length([posCommand', quatCommand])*2) = refVelocity;

% Concatenate states and reference with absolute time:
t = (plant.dt:plant.dt:H*plant.dt+plant.dt)';
r = [r, t];

for i=1:H
    L(i,1) = cost.fcn(cost, x(i,dyno), zeros(length(dyno)), a(i,:)', zeros(size(a,2)));    % compute rollout cost w/ energy penalty
end
L = L(1:H,1);               % Cost

y = x(2:H+1,1:nX);          % successor states  [H x nX]
x = [x(1:H,:) a(1:H,:)];    % observed  states  [H x nX+nU]

latent = x;                 % actual latent states unknown, so measured state is passed for congruency.

r = r(1:H+1,:);

end