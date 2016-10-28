function iiwaJointCommand(handles,q)
% Send position comamnd to iiwa joints
%
%%
handles.jointCommandMsg.Position.A1 = q(1);
handles.jointCommandMsg.Position.A2 = q(2);
handles.jointCommandMsg.Position.A3 = q(3);
handles.jointCommandMsg.Position.A4 = q(4);
handles.jointCommandMsg.Position.A5 = q(5);
handles.jointCommandMsg.Position.A6 = q(6);
handles.jointCommandMsg.Position.A7 = q(7);

send(handles.pubJointPosition, handles.jointCommandMsg);
end