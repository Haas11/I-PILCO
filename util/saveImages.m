
figHandles = findobj('Type','figure');
for i=1:length(figHandles);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) 
    saveas(gcf,strcat('eps',num2str(i)), 'eps');
    saveas(gcf,strcat('png',num2str(i)), 'png');
    clf(gcf);
end