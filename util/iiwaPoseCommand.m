function iiwaPoseCommand(handles, pose)
%iiwaPoseCommand send pose command to iiwa robot

%%
handles.poseCommandMsg.Pose.Position.X = pose(1);
handles.poseCommandMsg.Pose.Position.Y = pose(2);
handles.poseCommandMsg.Pose.Position.Z = pose(3);
handles.poseCommandMsg.Pose.Orientation.X = pose(4);
handles.poseCommandMsg.Pose.Orientation.Y = pose(5);
handles.poseCommandMsg.Pose.Orientation.Z = pose(6);
handles.poseCommandMsg.Pose.Orientation.W = pose(7);

send(handles.pubCartesianPose, handles.poseCommandMsg);

end

