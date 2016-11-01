
figHandles = findobj('Type','figure');
for i=1:length(figHandles);
<<<<<<< HEAD
    set(figHandles(i),'units','normalized','outerposition',[0 0 1 1])
%     print2eps(['eps',num2str(i),'.eps'], figHandles(i));
%     fixlines(['eps',num2str(i),'.eps']);
%     saveas(figHandles(i),['eps',num2str(i),'.eps'], 'epsc');
    saveas(figHandles(i),['png',num2str(i)], 'png');
=======
    set(figHandles(i),'units','normalized','outerposition',[0 0 1 1]) 
    saveas(figHandles(i),strcat('eps',num2str(i)), 'eps');
    saveas(figHandles(i),strcat('png',num2str(i)), 'png');
    clf(figHandles(i));
>>>>>>> d5abe78d437e0f83a378b815e55750f3b2c4c9c6
end