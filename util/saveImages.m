
figHandles = findobj('Type','figure');
for i=1:length(figHandles);
    set(figHandles(i),'units','normalized','outerposition',[0 0 1 1])
%     print2eps(['eps',num2str(i),'.eps'], figHandles(i));
%     fixlines(['eps',num2str(i),'.eps']);
%     saveas(figHandles(i),['eps',num2str(i),'.eps'], 'epsc');
    saveas(figHandles(i),['png',num2str(i)], 'png');
end