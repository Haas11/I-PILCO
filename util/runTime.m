% analysis of PILCOT runtime
clear all; clc; close all;
% policy:
%   nc = 5/10/25
%   input = [x y] / [x y dx dy]
%   output = [Kx Ky] / [Kx Ky Rx Ry]

% GP model:
%   nii = 100/200/300
%   dyni = [x y dx dy]
%   dyno = [x y dx dy fx fy]

% # of arrangements: 3*2*2*3*1*1 = 36


% dataset:
load('D:\victor\MATLAB\Results\ST\GP50\Data1\ST_28-Aug-2016_conGP-_ep-0p01_expl-0p25-K-125_H100.mat')
close all;
clear -x -y -r -dynmodel -mu0Sim -S0Sim

colorVec = {'r', 'b', 'm', 'c','y'};
lineSpecVec = {'o--', '^:', 's-.', 'x-', '<--'};

Vnc = [10 25 50 100 200];
Cpoli = {[1 2], [1 2 3 4], [1 2 3 4 5 6]};
lpoli = [2 4 6];

Vnii = [200 300 400 500];
Vdyno = [7 8 13 14 19 20];
Vdyni = [1 2 3 4];

numRuns = 1;

t_train = zeros(numRuns,length(Vnc),length(Cpoli),length(Vnii)); 
t_learn = zeros(numRuns,length(Vnc),length(Cpoli),length(Vnii)); 
t_pred  = zeros(numRuns,length(Vnc),length(Cpoli),length(Vnii)); 
% [run,p_nc,p_i,gp_nii]


%% Run
for runCount=1:numRuns
    for p_nc=1:length(Vnc)    % basis functions
        for p_i=1:length(Cpoli)     % inputs
            % GP variants
            for gp_nii=1:length(Vnii)  % inducing inputs
                % define variables
                nc = Vnc(p_nc);
                poli = Cpoli{p_i};
                nii = Vnii(gp_nii);
                
                % record execution times
                [t_train(runCount,p_nc,p_i,gp_nii), ...
                    t_learn(runCount,p_nc,p_i,gp_nii),...
                    t_pred(runCount,p_nc,p_i,gp_nii)] = timePILCOT_ST(nc, poli, nii, x, y, r, mu0Sim, S0Sim);
                p_nc 
                p_i
                runCount
                gp_nii
            end
        end
    end
end

% averages
t_trainAvg = mean(t_train(:,:,:,:),1);
t_learnAvg = mean(t_learn(:,:,:,:),1);
t_pedAvg = mean(t_pred(:,:,:,:),1);

save('cpuTimes','t_train','t_learn','t_pred');
save('avgCPUTimes','t_trainAvg','t_learnAvg','t_pedAvg');

load('cpuTimes');

%% plot training times
Vnii_inter = 200:20:500;
t_trainInter = zeros(length(Vnc),length(Vnii_inter),length(Cpoli));
n=2;
trainTimeFig = figure;
hold on;
grid minor;
% t_trainAvgNc = mean(t_train,2);
for i=1:length(Vnc)     % for #inducing inputs
    for k=1:length(Cpoli)   % for #policy inputs
        % inter/extrapolate values
        p_learn = polyfit(Vnii(1:end-1),squeeze(t_train(1,i,k,1:end-1))',n);   
        t_trainInter(i,:,k) = polyval(p_learn,Vnii_inter);
        
        % plot:
        plot(Vnii(1:end-1),squeeze(t_train(1,i,k,1:end-1))'./60,strcat(colorVec{i},lineSpecVec{k}(1)));       % data points
        hold on;
        plot(Vnii_inter,t_trainInter(i,:,k)'./60,strcat(colorVec{i},lineSpecVec{k}(2:end)))         % fitted curves
    end    
end
hb(1) = plot(0,0,'ko--','visible','off','MarkerSize',10,'LineWidth', 1);        % lpoli=2
hb(2) = plot(0,0,'k^:','visible','off','MarkerSize',10,'LineWidth', 1);         % lpoli=4
hb(3) = plot(0,0,'ks-.','visible','off','MarkerSize',10,'LineWidth', 1);        % lpoli=6

legend(hb,'s_{\pi}=[x y]','s_{\pi}=[x y dx dy]','s_{\pi}=[x y dx dy f_x f_y]','Interpreter','Tex');
title('\fontsize{13}Comparison in cpu time for VILTA GP training (\color{red}n_c=10\color{black}, \color{blue}n_c=25\color{black}, \color{magenta}n_c=50\color{black}, \color{cyan}n_c=100\color{black}, \color{yellow}n_c=200)','Interpreter','Tex')
xlabel('\fontsize{16}No. of Inducing Inputs n_{ii}','Interpreter','Tex'); ylabel('\fontsize{16}Computation time [minutes]');%% predicting

%% plot learning times
Vnc_inter = 0:10:300;
t_learnInter = zeros(length(Vnii),length(Vnc_inter),length(Cpoli));
n=2;
learnTimeFig = figure;
hold on;
grid minor;
for i=2:2%length(Vnii)     % for #inducing inputs
    for k=1:length(Cpoli)
        % inter/extrapolate values
        p_learn = polyfit(Vnc,t_learn(1,:,k,i),n);   
        t_learnInter(i,:,k) = polyval(p_learn,Vnc_inter);
        
        % plot:
        plot(Vnc,t_learn(1,:,k,i)./60,strcat(colorVec{k},lineSpecVec{k}(1)),'Linewidth',1.3);       % data points
        hold on;
        plot(Vnc_inter,t_learnInter(i,:,k)./60,strcat(colorVec{k},lineSpecVec{k}(2:end)),'Linewidth',1.3)         % fitted curves
    end    
end
hb(1) = plot(0,0,'ro--','visible','off','MarkerSize',10,'LineWidth', 1);        % nii=200
hb(2) = plot(0,0,'b^:','visible','off','MarkerSize',10,'LineWidth', 1);         % nii=300
hb(3) = plot(0,0,'ms-.','visible','off','MarkerSize',10,'LineWidth', 1);        % nii=400

legend(hb,'poli=[x y]','poli=[x y dx dy]','poli=[x y dx dy f_x f_y]');
title('\fontsize{13}Comparison in cpu time for VILTA policy update (n_{ii}=300\color{black})');%, \color{blue}n_{ii}=300\color{black}, \color{magenta}n_{ii}=400\color{black})','Interpreter','Tex')
xlabel('\fontsize{16}No. of Policy Kernels'); ylabel('\fontsize{16}Computation time [minutes]');%% predicting

%% plot learning times vs No of policy inputs
% [run,p_nc,p_i,gp_nii]
lpoli_inter = 1:1:10;
t_learnInter2 = zeros(length(Vnii),length(lpoli_inter),length(Vnii));
n=3;
learnTimeFig2 = figure;
hold on;
grid minor;
for i=1:length(Vnc)     % for #inducing inputs
    for k=1:length(Vnii)
        % inter/extrapolate values
        p_learn = polyfit(lpoli,squeeze(t_learn(1,i,:,k))',n);   
        t_learnInter2(i,:,k) = polyval(p_learn,lpoli_inter);
        
        % plot:
        plot(lpoli,squeeze(t_learn(1,i,:,k))'./60,strcat(colorVec{i},lineSpecVec{k}(1)));       % data points
        hold on;
        plot(lpoli_inter,t_learnInter2(i,:,k)./60,strcat(colorVec{i},lineSpecVec{k}(2:end)))         % fitted curves
    end    
end
hb(1) = plot(0,0,'ko--','visible','off','MarkerSize',10,'LineWidth', 1);        % nii=200
hb(2) = plot(0,0,'k^:','visible','off','MarkerSize',10,'LineWidth', 1);         % nii=300
hb(3) = plot(0,0,'ks-.','visible','off','MarkerSize',10,'LineWidth', 1);        % nii=400

legend(hb,'nii=200','nii=300','nii=400');
title('\fontsize{13}Comparison in cpu time for policy update (\color{red}n_c=10\color{black}, \color{blue}n_c=25\color{black}, \color{magenta}n_c=50\color{black}, \color{cyan}n_c=100\color{black}, \color{yellow}n_c=200)','Interpreter','Tex')
xlabel('\fontsize{16}No. of Policy Inputs (n_{\pi})','Interpreter','Tex'); ylabel('\fontsize{16}Computation time [minutes]');%% predicting
    
%% plot prediction times
Vnii_inter = 100:20:500;
t_predInter = zeros(length(Vnc),length(Vnii_inter),length(Cpoli));
n=4;
% [run,p_nc,p_i,gp_nii]

predTimeFig = figure;
hold on;
grid minor;
for i=1:length(Vnc)     % for #inducing inputs
    for k=1:length(Cpoli)
        % inter/extrapolate values
        p_learn = polyfit(Vnii,squeeze(t_pred(1,i,k,:))',n);
        t_predInter(i,:,k) = polyval(p_learn,Vnii_inter);
        
        % plot:
        plot(Vnii,squeeze(t_pred(1,i,k,:))',strcat(colorVec{i},lineSpecVec{k}(1)));       % data points
        hold on;
        plot(Vnii_inter,t_predInter(i,:,k)',strcat(colorVec{i},lineSpecVec{k}(2:end)))         % fitted curves
    end    
end
hb(1) = plot(0,0,'ko--','visible','off','MarkerSize',10,'LineWidth', 1);        % lpoli=2
hb(2) = plot(0,0,'k^:','visible','off','MarkerSize',10,'LineWidth', 1);         % lpoli=4
hb(3) = plot(0,0,'ks-.','visible','off','MarkerSize',10,'LineWidth', 1);        % lpoli=6

legend(hb,'s_{\pi}=[x y]','s_{\pi}=[x y dx dy]','s_{\pi}=[x y dx dy f_x f_y]','Interpreter','Tex');
title('\fontsize{13}Comparison in cpu time for VILTA GP prediction (\color{red}n_c=10\color{black}, \color{blue}n_c=25\color{black}, \color{magenta}n_c=50\color{black}, \color{cyan}n_c=100\color{black}, \color{yellow}n_c=200)','Interpreter','Tex')
xlabel('\fontsize{16}No. of Inducing Inputs n_{ii}','Interpreter','Tex'); ylabel('\fontsize{16}Computation time [s]');%% predicting

%%
saveImages
