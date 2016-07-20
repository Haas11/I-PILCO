function genLegend(Marker,text,Fignum,FONTSIZE,MARKERSIZE,LWD)
% Generate custom legend
%
% Example:
% clc;clear;close all;fclose all
% Marker = {'r -*' , 'b -o' , 'g --d'};
% text = {'First Fun ' , 'Second Fun ' , 'Third Fun ' };
% FONTSIZE = 6;
% Fignum = 26;

if length(Marker)~=length(text)
    error('You should have the same number of markers as texts');
end

%% FIGURE
set(0,'defaultaxesfontsize',FONTSIZE);
figure(Fignum);
% clf;
hold on;
% set(gcf,'OuterPosition',[60  90 200  250],'PaperPositionMode','auto');

Leg_TEXTS = [];
for i = 1:length(Marker)
    x = 0:0.1:1;
    y = 3*randn(1)*sin(x);
    p = plot(x,y,Marker{i},'MarkerSize',MARKERSIZE,'LineWidth',LWD);
    Leg_TEXTS = [ Leg_TEXTS ; {text{i}} ];
end

legend(Leg_TEXTS);
figure(Fignum);
set([gca;get(gca, 'children')], 'visible', 'off');

end