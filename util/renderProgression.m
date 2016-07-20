% Script to render movie of certain metrics througout iterations.
clear all; clc; close all;
global vert fac

%% Input-Output Definitions
dataSet = 'conGP-15_ep-0p01_expl-0p2_25-Apr-2016-9_H100.mat';

filenameModel = strcat(dataSet,'_Model');
filenameCost = strcat(dataSet,'_Cost');
filenameAnim = strcat(dataSet,'_Anim');
filenameAct = strcat(dataSet,'_Act');

modelDirectory = 'Model';
costDirectory = 'Cost';
animeDirectory = 'anim';
actDirectory = 'action';

modelFigureNo = 1;
costFigureNo = 2;
animFig = 3;
actFig = 4;

vidTime = 20; %[s]

animateModel = Animate(modelDirectory,250);
animateCost  = Animate(costDirectory, 250);
animateVid = Animate(animeDirectory,  250);
animateAct = Animate(actDirectory, 250);

%% Record Figures
load(dataSet);

for kk=1:length(M)
    close all;
    figure(modelFigureNo);
    set(gcf, 'Position', get(0,'Screensize'));
    ldyno = length(dyno);
    for i=1:ldyno       % plot the rollouts on top of predicted error bars
        %         subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
        subplot(2,3,i); hold on;
        errorbar( 0:length(M{kk}(i,:))-1, M{kk}(i,:), ...
            2*sqrt(squeeze(Sigma{kk}(i,i,:))) );
        
        for ii=1:Ntest
            stairs( 0:size(lat{ii}(:,dyno(i)),1)-1, lat{ii}(:,dyno(i)), 'g' );        % recorded latent states in multiple robustness test-rollouts
        end
        stairs( 0:size(latent{kk+J}(:,dyno(i)),1)-1, latent{kk+J}(:,dyno(i)),'r');      % recorded latent states in apply_controller roll-out
        title(dynoTitles{i});
        axis tight
    end
    drawnow;
    animateModel.add();
    
    % plot Cost:
    figure(costFigureNo);
    set(gcf, 'Position', get(0,'Screensize'));
    errorbar(0:H,fantasy.mean{kk},2*fantasy.std{kk});
    hold on
    stairs(1:length(realCost{J+kk}),realCost{J+kk},'r');
    title('Predicted (uncertain) & Rollout (deterministic) Immediate Cost');
    xlabel('Timestep');     ylabel('Immediate Cost');
    drawnow;
    animateCost.add();    
end

%% Best policy:
[~,best] = min(realAcumCost(1:11);
best = 11;
[~,~,~, latent{best}] = ...
    my_rollout(mu0, policy, HH, plant, robot); %#ok<*IJCL>

q_sim = latent{best}(:,1:robot.n);
a     = latent{best}(:,end-numel(policy.maxU)+1:end);
H     = length(q_sim);

figure(actFig);
set(gcf, 'Position', get(0,'Screensize'));
subplot(1,2,1);
robot.plot(q_sim(1,1:robot.n));     hold on;
patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');     hold on;    
for i=1:H
    subplot(1,2,1);
    title('Rollout trajectory')
    set(gca,'fontsize', 16);
    robot.plot(q_sim(i,1:robot.n));
    
    subplot(1,2,2);
    stairs(1:i,a(1:i,1),'b');       hold on;     grid on;
    stairs(1:i,a(1:i,2),'r');
    legend('KpX','KpY');    title('Final policy');    xlabel('Timestep');    ylabel('Kp_{XY} [N/m]');
    set(gca,'fontsize', 16);
    axis([0 H 0 max(policy.maxU)*2]);
    drawnow;
    animateAct.add();                        % export frame as .png   
end

%% Render Videos
inputFramerate = 1/2;
outputFramerate = 1;
renderCommand = ['ffmpeg -framerate ',num2str(inputFramerate),' -i ',modelDirectory,'/%04d.png -c:v libx264 -r ',num2str(outputFramerate),' -pix_fmt yuv420p ',filenameModel,'.mp4'];
system(renderCommand);   % convert /.png images to .mp4 movie

renderCommand = ['ffmpeg -framerate ',num2str(inputFramerate),' -i ',costDirectory,'/%04d.png -c:v libx264 -r ',num2str(outputFramerate),' -pix_fmt yuv420p ',filenameCost,'.mp4'];
system(renderCommand);   % convert /.png images to .mp4 movie
% 
% renderCommand = ['ffmpeg -framerate ',num2str(inputFramerate),' -i ',animeDirectory,'/%04d.png -c:v libx264 -r ',num2str(outputFramerate),' -pix_fmt yuv420p ',filenameAct,'.mp4'];
% system(renderCommand);   % convert /.png images to .mp4 movie

inputFramerate = round(H/vidTime);
outputFramerate = 25;
renderCommand = ['ffmpeg -framerate ',num2str(inputFramerate),' -i ',actDirectory,'/%04d.png -c:v libx264 -r ',num2str(outputFramerate),' -pix_fmt yuv420p ','test.mp4'];
system(renderCommand);   % convert /.png images to .mp4 movie