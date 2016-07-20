figure;
for kk=1:length(M)
    set(gcf, 'Position', get(0,'Screensize'));
    ldyno = length(dyno);
    clf;
    for i=1:ldyno       % plot the rollouts on top of predicted error bars
        %         subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
        subplot(3,3,i); hold on;
        errorbar( 0:length(M{kk}(i,:))-1, M{kk}(i,:), ...
            2*sqrt(squeeze(Sigma{kk}(i,i,:))) );
        
        for ii=1:Ntest
            stairs( 0:size(testLat{ii}(:,dyno(i)),1)-1, testLat{ii}(:,dyno(i)), 'g' );        % recorded latent states in multiple robustness test-rollouts
        end
        stairs( 0:size(latent{kk+J}(:,dyno(i)),1)-1, latent{kk+J}(:,dyno(i)),'r');      % recorded latent states in apply_controller roll-out
        title(strcat(dynoTitles{i},num2str(kk)));
        axis tight
    end
    drawnow;
    pause;
end

%%
figure;
for kk=1:length(M)
    q_sim = latent{kk}(:,1:robot.n);
    if ~ishandle(5)         % robot animation
        figure(5);
    else
        set(0,'CurrentFigure',5);
    end
    clf(5);
%     patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
    robot.plot(q_sim);
end