%% Linespecs:
color = [255, 100, 0;
    255,255,0;
    0, 255, 0]./255;
markerStyle = {'s','d','o'};

%% Plot Fantasy Costs
if ~ishandle(10)         % cost iterations
    figure(10);
else
    set(0,'CurrentFigure',10);
end
hold on; grid on;
clf;

% variances:
if j==1
    errorbar(1:j, fantasy.acumMean(1:j), 2*fantasy.acumStd(1:j), 'r','LineWidth',1.1);
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
        'LineWidth',1.5, 'visible', 'off','MarkerSize',10);   % dummy plot for legend
end
hold on;

% real variances:
h1 = errorbar(0,realWorld.mean(1),2*realWorld.std(1),'b--','LineWidth',1.2);
for i=1:j
    errorbar(i,realWorld.mean(i+1),2*realWorld.std(i+1),'b','LineWidth',1.2,'AlignVertexCenters','on');
end

%%
hb(1) = plot(0,0,'r','visible','off','MarkerSize',10,'LineWidth', 1.5);     % predicted
hb(2) = plot(0,0,'b','visible','off','MarkerSize',10,'LineWidth', 1.5);     % recorded
axis auto
ax = gca; ax.XTick = 0:1:j;
legend(hb,'Simulated (95% conf.)','Recorded (95%)','None Success','Some Success','Full Success','Location','NorthEast');
title(strcat('\fontsize{14}Predicted and Recorded Cost per Iteration.  (\color{red}K=',num2str(K),'\color{black} and \color{blue}Ntest=',num2str(Ntest),'\color{black})'));
xlabel('\fontsize{14}No. of Trials [#]'); ylabel('\fontsize{14}Cost Distribution [\mu, \sigma]');