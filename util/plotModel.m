if ~ishandle(4)
    figure(4);
else
    set(0,'CurrentFigure',4)
end
clf(4);
if dyno(end)~=stateLength
    ldyno = length(dyno);
else
    ldyno = length(dyno)-1;     % don't plot time
end
for i=1:ldyno
    hold on;
    if ldyno==6
        subplot(3,2,i);
    else
        subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i);
    end
    
    % Reference:
    if dyno(i)~=stateLength
        plot((1:size(rr,1))*plant.dt, rr(:,i),'k:');
    end
    hold on;
    
    % Sparse Model:
    %             errorbar( (0:length(M{j}(i,:))-1)*plant.dt, M{j}(i,:), ...
    %                 2*sqrt(squeeze(Sigma{j}(i,i,:))),'b','AlignVertexCenters','on'); %#ok<*IJCL>
    shadedplot((0:length(M{j}(i,:))-1)*plant.dt, M{j}(i,:) - 2*sqrt(squeeze(Sigma{j}(i,i,:)))',...
        M{j}(i,:) + 2*sqrt(squeeze(Sigma{j}(i,i,:)))', [1 0.75 0.75], 'r'); %#ok<*IJCL>
    plot((0:length(M{j}(i,:))-1)*plant.dt, M{j}(i,:),'r--');
    
    % Full model:
    if compareToFullModel && ~isempty(Mfull{j})
        errorbar(0:length(Mfull{j}(i,:))-1, Mfull{j}(i,:), ...
            2*sqrt(squeeze(Sfull{j}(i,i,:))), 'y');
    end
    
    % Test:
    for ii=1:Ntest
        stairs((0:size(testLati{ii}(:,dyno(i)),1)-1)*plant.dt, testLati{ii}(:,indices(dyno(i))), 'g','LineWidth',1.1);        % recorded latent states in multiple robustness test-rollouts
    end

    % Model trial:
    stairs((0:size(latent{j+J}(:,dyno(i)),1)-1)*plant.dt, latent{j+J}(:,indices(dyno(i))),'Color', [0 0 255]./255,'LineWidth',1.2);      % recorded latent states in apply_controller roll-out
        
    % Inducing inputs Locations:
    if i <= length(dyni) && numel(dynmodel.induce) ~= 0
        plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
    end
    
    title(strcat('\fontsize{14}',dynoTitles{i}),'Interpreter','Tex');
    if i==5 || i==6; xlabel('\fontsize{14}Time [s]'); end
    if i==1
        ml(1) = plot(0,0,'k:','visible','off','MarkerSize',10,'LineWidth', 1);     % predicted
        ml(2) = plot(0,0,'r--','visible','off','MarkerSize',10,'LineWidth', 1.5);     
        ml(3) = plot(0,0,'b','visible','off','MarkerSize',10,'LineWidth', 1);     
        ml(4) = plot(0,0,'g--','visible','off','MarkerSize',10,'LineWidth', 1);     
        
        if compareToFullModel && ~isempty(Mfull{j})
            legend('Reference','Planning Model','Full Model','Data Trial','Test Trials','Location','Best');
        else
            legend(ml,'Reference','Planning Model','Data Trial','Test Trials','Location','Best');
        end
    end
    axis tight; grid on
end
ha = axes('Position',[0 0 1 0.98],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\fontsize{16}Long-term Bayesian Inference (95% conf.) + Roll-Out Data','HorizontalAlignment','center','VerticalAlignment', 'top')
drawnow;