function [handles, response] = setRobotStiffDamp(handles, mode, stiffness, damping, time_out, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% Code
handles.config.Mode.Mode = mode;

if mode == 1
    cartDim = {'X','Y','Z','A','B','C'};
    
    % set desired stiffness and damping
    for i=1:6
        handles.config.Mode.CartesianStiffness.Stiffness.(cartDim{i}) = stiffness(i);
        handles.config.Mode.CartesianDamping.Damping.(cartDim{i}) = damping(i);
    end
    
    % set nullspace dynamics
    if nargin == 6
        handles.config.Mode.NullspaceStiffness = varargin{1};
    elseif nargin == 7
        handles.config.Mode.NullspaceStiffness = varargin{1};
        handles.config.Mode.NullspaceDamping = varargin{2};
    end
            
elseif mode==0
    jointDim = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'};
    for i=1:7
        % stiffness
        handles.config.Mode.JointStiffness.Stiffness.(jointDim{i}) = stiffness(i);
        
        % damping
        handles.config.Mode.JointDamping.Damping.(jointDim{i}) = damping(i);
    end
end

% send service request to Cabinet
response = call(handles.client, handles.config, 'Timeout', time_out);
end


