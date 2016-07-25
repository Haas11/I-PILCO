figure;
for kk=1:length(M)
    set(gcf, 'Position', get(0,'Screensize'));
    ldyno = length(dyno);
    clf;
    for i=1:ldyno       % plot the rollouts on top of predicted error bars
        subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
        
        % Full/Sparse model:
        errorbar(0:length(M{kk}(i,:))-1, M{kk}(i,:), ...
            2*sqrt(squeeze(Sigma{kk}(i,i,:))),'b');
        
        % Full model:
        if compareToFullModel && ~isempty(Mfull{kk})
            errorbar(0:length(Mfull{kk}(i,:))-1, Mfull{kk}(i,:), ...
                2*sqrt(squeeze(Sfull{kk}(i,i,:))), 'y');
        end
        
        % Test trials:
        for ii=1:Ntest
            stairs(0:size(testLat{ii}(:,dyno(i)),1)-1, testLat{ii}(:,indices(dyno(i))), 'g' );        % recorded latent states in multiple robustness test-rollouts
        end
        
        % Real trial:
        stairs(0:size(latent{kk+J}(:,dyno(i)),1)-1, latent{kk+J}(:,indices(dyno(i))),'r');      % recorded latent states in apply_controller roll-out
        
        % All input locations:
        if i <= length(dyni)
            plot(ones(size(dynmodel.inputs,1),1)*H,dynmodel.inputs(:,i),'cx');
        end
        
        % Inducing inputs Locations:
        if i <= length(dyni) && ~isempty(Mfull{kk})
            plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
        end
        
        title(dynoTitles{i});
        if i==1
            if isempty(Mfull{kk})
                legend('Full Model','Trial','Training Inputs','Sparse Inputs','Location','Best');
            else
                legend('Sparse Model','Full Model','Trial','Training Inputs','Sparse Inputs','Location','Best');
            end
        end
        axis tight
    end
    drawnow;
    pause;
end

%%
figure;
for kk=1:length(M)
    plot(dynmodel.inputs(1:kk*H,1),dynmodel.targets(1:kk*H,1),'b--+');
    hold on
    errorbar(M{kk}(1,1:end-1)',M{kk}(1,2:end)',2*sqrt(squeeze(Sigma{kk}(1,1,2:end))),'rx');
    hold on
%     plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
    title(strcat('Iteration: ',num2str(kk)));
    xlabel('Inputs')
    ylabel('Outputs');
    legend('Observed Data','GP Model','Location','Best');
    drawnow;
    pause;
end


%% Full Model:
figure;
set(gcf, 'Position', get(0,'Screensize'));
ldyno = length(dyno);
clf;
for i=1:ldyno       % plot the rollouts on top of predicted error bars
    %         subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
    subplot(3,3,i); hold on;
    errorbar(0:length(Mfull(i,:))-1, Mfull(i,:), ...
        2*sqrt(squeeze(Sfull(i,i,:))), 'r');
    for ii=1:Ntest
        stairs( 0:size(testLat{ii}(:,dyno(i)),1)-1, testLat{ii}(:,dyno(i)), 'g' );        % recorded latent states in multiple robustness test-rollouts
    end
    stairs( 0:size(latent{kk+J}(:,dyno(i)),1)-1, latent{kk+J}(:,dyno(i)),'r');      % recorded latent states in apply_controller roll-out
    title(strcat(dynoTitles{i},num2str(kk)));
    axis tight
end
drawnow;

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
