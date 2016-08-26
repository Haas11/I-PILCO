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

% variance:
if Ntest > 1
    h1 = errorbar(0,realWorld.mean(1),2*realWorld.std(1),'b--','LineWidth',1.2);
    for i=1:j
        errorbar(i,realWorld.mean(i+1),2*realWorld.std(i+1),'b','LineWidth',1.2,'AlignVertexCenters','on');
    end
end
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
% ax.YLim = [fantasy.acumMean(1)-2.1*fantasy.acumStd(1), fantasy.acumMean(1)+2.1*fantasy.acumStd(1)];

legend(hb,'Simulated (95% conf.)','Recorded (95%)','None Success','Some Success','Full Success','Location','NorthEast');
title(strcat('\fontsize{14}Predicted and Recorded Accumulated Cost Distributions (\color{red}K=',num2str(K),'\color{black} and \color{blue}Ntest=',num2str(Ntest),'\color{black})'));
xlabel('\fontsize{14}Iteration #');   ylabel('\fontsize{14}Cost Distribution [\mu, \sigma]');
set(ax,'YLim',[40 ax.YLim(2)]);

