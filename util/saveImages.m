figHandles = findobj('Type','figure');
for i=1:length(figHandles);
    set(figHandles(i),'units','normalized','outerposition',[0 0 1 1]) 
    saveas(figHandles(i),strcat('eps',num2str(i)), 'epsc');
    saveas(figHandles(i),strcat('png',num2str(i)), 'png');
end