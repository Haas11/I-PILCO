actionFig = figure(12);
modelFig = figure(23);
clf(modelFig);
clf(actionFig);
for kk=1:length(M)
    set(0,'CurrentFigure',modelFig);
    %     set(gcf, 'Position', get(0,'Screensize'));
    ldyno = length(dyno);
    clf;
    for i=1:ldyno       % plot the rollouts on top of predicted error bars
        %         subplot(ceil(ldyno/sqrt(ldyno)),ceil(sqrt(ldyno)),i); hold on;
        subplot(3,4,i); hold on;
        
        % Reference:
%         plot(rr(:,i),'k:');
        
        % Sparse Model:
        errorbar( 0:length(M{kk}(i,:))-1, M{kk}(i,:), ...
            2*sqrt(squeeze(Sigma{kk}(i,i,:))),'b' );
        
        % Full Model:
        if ~isempty(Mfull{kk})
            errorbar( 0:length(Mfull{kk}(i,:))-1, Mfull{kk}(i,:), ...
                2*sqrt(squeeze(Sfull{kk}(i,i,:))),'y' );
        end
        
        % Test:
        for n=1:Ntest
            stairs( 0:size(testLat{kk}{n}(:,dyno(i)),1)-1, testLat{kk}{n}(:,dyno(i)), 'g' );        % recorded latent states in multiple robustness test-rollouts
        end
        
        % Trial:
        stairs( 0:size(latent{kk+J}(:,dyno(i)),1)-1, latent{kk+J}(:,dyno(i)),'r');      % recorded latent states in apply_controller roll-out
        
        % Inducing inputs Locations:
        if i <= length(dyni) && ~isempty(Mfull{kk})
            plot(zeros(nii,1),dynmodel.induce(:,i),'kx');
        end
        
        title(strcat(dynoTitles{i},num2str(kk))); grid on;
        axis tight
        if i==1
            if ~isempty(Mfull{kk})
                legend('Reference','Planning Model','Full Model','Test','Trial','Sparse Inputs','Location','Best');
            else
                legend('Reference','Planning Model','Test','Trial','Sparse Inputs','Location','Best');
            end
        end
    end
    
        set(0,'CurrentFigure',actionFig);
        clf;
        for i=1:2
            subplot(1,2,i);
            
            % sim:
            errorbar( 0:length(Mcon{kk}(i,:))-1, Mcon{kk}(i,:), ...
                2*sqrt(squeeze(Scon{kk}(i,i,:))), 'b');
            hold on;
            
            % model:
            stairs(x(1+H*(kk):H*(kk+1),end-Du+i),'r');
            
            % test:
%             for n=1:Ntest
%                 stairs(testLat{kk}{n}(:,end-Du+i),'g');
%             end
            
            axis tight; grid on;
            xlabel('Timestep');     ylabel(actionTitles{i});           
            title(strcat('Iteration: ',num2str(kk)));
            if i==1
                legend('Simulated','Model Trial','Test Trial','Location','Best');
            end
        end
    
    drawnow;
    pause;
end

% %%
% figure;
% for kk=1:length(M)
%     q_sim = latent{kk}(:,1:robot.n);
%     if ~ishandle(5)         % robot animation
%         figure(5);
%     else
%         set(0,'CurrentFigure',5);
%     end
%     clf(5);
%     %     patch('Vertices',vert,'Faces',fac,'FaceVertexCData',hsv(6),'FaceColor','flat');
%     robot.plot(q_sim);
% end