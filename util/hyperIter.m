
dim = 4;  % target dimension to be plotted
close all
for dataIter=15:23
    dataset = strcat('_01-Aug-2016_conLin-_ep-0p001_expl-0p2-',num2str(dataIter),'_H50.mat');
    load(dataset);
    l = exp(dynmodel.hyp);
    lhyp = size(dynmodel.hyp,1);
    if ~ishandle(20)         % robot animation
        figure(20);
    else
        set(0,'CurrentFigure',20);
    end
    hold on;
    for hyperIter=1:lhyp
        subplot(2,5,hyperIter); hold on;
        grid on;
        
        plot(dataIter,l(hyperIter,dim),'x','MarkerSize',8,'LineWidth',2);
        title(strcat(hyperTitles{hyperIter},' (',num2str(dataIter),')'));
    end
    drawnow;
    
    figs = get(0,'children');
    figs(figs == 20) = []; % delete your current figure from the list
    close(figs);

%     pause;
    
end
