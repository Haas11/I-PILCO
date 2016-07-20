function animateVideo(robot, q_sim, filename)
% animateVideo animates robot trajectory and outputs .mp4 video.
global vert fac
H=length(q_sim);

anime = Animate('2Dpeg',250);

f=figure;
set(gcf, 'Position', get(0,'Screensize'));
clf(f);
robot.plot(q_sim(1,1:robot.n));
hold on
patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
for i=1:H
    robot.plot(q_sim(i,1:robot.n));
    anime.add();                        % export frame as .png
end

renderCommand = strcat('ffmpeg -framerate 10 -i 2Dpeg/%04d.png -c:v libx264 -r 25 -pix_fmt yuv420p ',filename,'.mp4');
system(renderCommand);   % convert /.png images to .mp4 movie

end