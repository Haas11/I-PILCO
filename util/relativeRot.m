function [EUL_de, dEUL_de, T_e, dT_e] = relativeRot(T_ee, T_d, prevT_e, prevEUL_de, dt)
% Computes modified euler angles EUL_de
%
% Outputs
%       EUL_de        - Euler angles of relative error
%       T_e(EUL_de)   - Transformation matrix of ee-frame as a function of de
%     relative euler angles
%       dT_e          - Time derivative of transf. matrix.

% Inputs:
%       T_ee          - Transformation matrix of current ee-frame       (referred to base)
%       T_d           - Transformation matrix of desired frame          (referred to base)
%       dRel          - Error state vector of positions and Euler angles
%
% Code by Victor van Spaandonk
% May 2016

%% Prepare inputs
R_e = tr2rot(T_ee);         % Rotation of ee-frame          (referred to base)
R_d = tr2rot(T_d);          % Rotation of desired frame     (referred to base)

%% Code
Re_d = R_e'*R_d;                        % mutual orientation of ee-frame w.r.t. desired frame
EUL_de = rotm2eul(Re_d,'ZYX');          % extract (non-standard) EULER angles (yaw, pitch, roll)

T = cardatol(EUL_de,3,2,1);              % Transformation matrix relating Euler angle derivatives to angular velocity   (Z=3,Y=2,X=1)
T_e = R_e*T;                             % Relative transformation                      (referred to base-frame)

dT_e = (T_e - prevT_e)./dt;
dEUL_de = (EUL_de - prevEUL_de)./dt;


%deltaOmegae_de = T*dEUL_de;             % relative angular velocity                    (referred to ee-frame)
% T_d = R_d*T;
% a_o = domega_d + dT_e*dEUL_de* + T_e*(Kd*dEUL_de + Kp*EUL_de);
% domega_d = T_d*ddEUL_d + dT_d*dEUL_d;   % desired angular acceleration
% --> WILL USUALLY BE ZERO!

%%
% Angular velocities
%T = rpy2tr(EUL_de, 'zyx');               % relative transformation matrix               (referred to ee-frame)

% dT_e = delta2tr([dPos_de; dEUL_de]);     % Derivative of relative transformation matrix (referred to base frame)
% omegaTilde = dR*R';
% omega = vex(omegaTilde);
% dRe_d = skew(deltaOmegae_de)*Re_d;      % time derivative of rotation matrix


% dR = skew(omega)*R;
% omega  = T*dEUL;
% domega = T*ddEUL + dT*dEUL;