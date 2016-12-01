%% Linespecs:
color = [255, 100, 0;
    255,255,0;
    0, 255, 0]./255;
markerStyle = {'s','d','o'};

%% Plot Fantasy Costs
if ~ishandle(10)         % cost iterations
    figure(10);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9])
else
    set(0,'CurrentFigure',10);
end
hold on; grid on;
clf;

% variances:
if isfield(fantasy,'acumMean')
    if j==1 %#ok<*IJCL>
        errorbar(1, fantasy.acumMean(1), 2*fantasy.acumStd(1), 'r','LineWidth',1.2);
        grid on;
    else
        shadedplot(1:j, fantasy.acumMean(1:j)-2*fantasy.acumStd(1:j), fantasy.acumMean(1:j)+2*fantasy.acumStd(1:j), [1 0.75 0.75], 'r'); % 95[%] confidence intervals
    end
    % means:
    hold on;
    plot(1:j, fantasy.acumMean(1:j),'r--','LineWidth',1.1);
    
    %% Plot Real World Costs
    failMedSuc{1} = find(scoreCard(1:j+1)==0);  % none success
    failMedSuc{2} = find(scoreCard(1:j+1)==1);  % some success
    failMedSuc{3} = find(scoreCard(1:j+1)==2);  % full success
    
    % means:
    for k=1:3
        for i=1:length(failMedSuc{k})
            plot(failMedSuc{k}-1, realWorld.mean(failMedSuc{k}),...
                markerStyle{k},...
                'MarkerEdgeColor',[0 0 1],...
                'MarkerFaceColor',color(k,:),...
                'MarkerSize',10,...
                'LineWidth',1.5);
            hold on
        end
        hb(k+2) = plot(0,0,markerStyle{k},...
            'MarkerEdgeColor',[0 0 1],...
            'MarkerFaceColor',color(k,:),...
            'LineWidth',1.5, 'visible', 'off','MarkerSize',10);   %#ok<*SAGROW> % dummy plot for legend
    end
    hold on;
    axis tight
    
    % variance:
    h1 = errorbar(0,realWorld.mean(1),2*realWorld.std(1),'b--','LineWidth',1.2);
    for i=1:j
        errorbar(i,realWorld.mean(i+1),2*realWorld.std(i+1),'b','LineWidth',1.2,'AlignVertexCenters','on');
    end
    
    ax = gca;
    if numel(dynmodel.induce) ~= 0
        plot([(nii/H)-J, (nii/H)-J], ax.YLim,'k:');
        [x, y] = dsxy2figxy(ax, [(nii/H)-J-1, (nii/H)-J],[(ax.YLim(2) - ax.YLim(1))/7, (ax.YLim(2) - ax.YLim(1))/5]);
        anno = annotation('textarrow',x,y,'String', 'Beginning of sparse model');
        anno.FontSize = 14;
    end
    
    %%
    hb(1) = plot(0,0,'r','visible','off','MarkerSize',10,'LineWidth', 1.5);     % predicted
    hb(2) = plot(0,0,'b','visible','off','MarkerSize',10,'LineWidth', 1.5);     % recorded
    
    ax.XTick = 0:1:j;
%     set(ax,'YLim',[min(fantasy.acumMean(1:j)-2.1*fantasy.acumStd(1:j)), max(fantasy.acumMean(1:j)+2.1*fantasy.acumStd(1:j))+1e-10]);
    
    legend(hb,'Simulated (95% conf.)','Recorded (95%)','None Success','Some Success','Full Success','Location','NorthEast');
    title(strcat('\fontsize{14}Predicted and Recorded Cost per Iteration.  (\color{red}K=',num2str(K),'\color{black} and \color{blue}Ntest=',num2str(Ntest),'\color{black})'));
    xlabel('\fontsize{14}No. of Trials [#]');   ylabel('\fontsize{14}Cost Distribution [\mu, \sigma]');
else
    if ~ishandle(10)         % cost iterations
        figure(10);
    else
        set(0,'CurrentFigure',10);
    end
    clf(10);
    hold on; grid on;
    failMedSuc{1} = find(insertSuccess(1:j+J)==0);
    failMedSuc{2} = find(insertSuccess(1:j+J)==1);
    failMedSuc{3} = find(insertSuccess(1:j+J)==2);
    color = {'rx','b+','go'};
    
    % plot recorded accumulated cost:
    for k=1:3
        for i=1:length(failMedSuc{k})
            plot(failMedSuc{k},realAcumCost(failMedSuc{k}),color{k},'MarkerSize',10,'LineWidth',1.5);
            hold on
        end
        hb(k) = plot(0,0,color{k}, 'visible', 'off','MarkerSize',10,'LineWidth',1.5);   % dummy plot for legend
    end
    plot(realAcumCost(1:j+J),'r--','LineWidth',0.1);
    
    % plot predicted accumulated costs:
    for i=1:j, errorbar(i+J,sum(fantasy.mean{i}), 2*sum(fantasy.std{i}), 'k*');    end
    hb(4) = plot(0,0,'k*','visible','off','MarkerSize',10,'LineWidth', 1.5);
    legend(hb,'Aborted','Failed','Successful','Predicted');
    title('Accumulated Rollout Cost');   xlabel('Learning iteration');   ylabel('Total Cost');
    ax = gca; ax.XTick = 1:1:J+j;
    
    hb(k+2) = plot(0,0,markerStyle{k},...
        'MarkerEdgeColor',[0 0 1],...
        'MarkerFaceColor',color(k,:),...
        'LineWidth',1.5, 'visible', 'off','MarkerSize',10);   %#ok<*SAGROW> % dummy plot for legend
    hold on;
    
    % variance:
    if Ntest > 1
        h1 = errorbar(0,realWorld.mean(1),2*realWorld.std(1),'b--','LineWidth',1.2);
        for i=1:j
            errorbar(i,realWorld.mean(i+1),2*realWorld.std(i+1),'b','LineWidth',1.2,'AlignVertexCenters','on');
        end
    end
% <<<<<<< HEAD
% end
% % connecting line
% plot(0:j,realWorld.mean(1:j+1),'b-.','LineWidth',0.5);
% 
% 
% %%
% hb(1) = plot(0,0,'r--','visible','off','MarkerSize',10,'LineWidth', 1.5);     % predicted
% hb(2) = plot(0,0,'b-.','visible','off','MarkerSize',10,'LineWidth', 1.5);     % recorded
% 
% ax = gca;
% 
% ax.XTick = 0:1:j; 
% ax.YLim = [fantasy.acumMean(1)-2.1*fantasy.acumStd(1), fantasy.acumMean(1)+2.1*fantasy.acumStd(1)];
% ax.FontSize = 16;
% % start of inducing inputs:
% % if numel(dynmodel.induce) ~= 0
% %     plot([(nii/H)-J, (nii/H)-J],ax.YLim,'k:');
% % end
% axis tight
% legend(hb,'Simulated (95% conf.)','Recorded (95%)','None Success','Some Success','Full Success','Location','NorthEast');
% title(strcat('Predicted and Recorded Roll-Out Distributions (\color{red}K=',num2str(K),'\color{black} and \color{blue}Ntest=',num2str(Ntest),'\color{black})'));
% xlabel('Iteration #');   ylabel('Expected Returns');
% 
% % set(ax,'YLim',[10 ax.YLim(2)]);
% 
% =======
    % connecting line
    plot(0:j,realWorld.mean(1:j+1),'b-.','LineWidth',0.5);
    
    ax = gca;
    % start of inducing inputs:
    if numel(dynmodel.induce) ~= 0
        plot([(nii/H)-J, (nii/H)-J],ax.YLim,'k:');
    end
    %%
    hb(1) = plot(0,0,'r--','visible','off','MarkerSize',10,'LineWidth', 1.5);     % predicted
    hb(2) = plot(0,0,'b-.','visible','off','MarkerSize',10,'LineWidth', 1.5);     % recorded
    
    ax.XTick = 0:1:j;
    set(ax,'YLim',[min(fantasy.acumMean(1:j)-2.1*fantasy.acumStd(1:j)), max(fantasy.acumMean(1:j)+2.1*fantasy.acumStd(1:j))]);
    
    legend(hb,'Simulated (95% conf.)','Recorded (95%)','None Success','Some Success','Full Success','Location','NorthEast');
    title(strcat('\fontsize{14}Predicted and Recorded Accumulated Cost Distributions (\color{red}K=',num2str(K),'\color{black} and \color{blue}Ntest=',num2str(Ntest),'\color{black})'));
    xlabel('\fontsize{14}Iteration #');   ylabel('\fontsize{14}Cost Distribution [\mu, \sigma]');
end
% >>>>>>> d5abe78d437e0f83a378b815e55750f3b2c4c9c6
