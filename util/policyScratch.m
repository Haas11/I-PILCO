% Policy scratch
close all; clear all; clc
%% shared
amin =  -100;    % limits of the unsquashed control function
amax =   200;

x = linspace(-15,15,50);

maxLinOut = pi/2;

maxU = 100;

target = 50;
target_scaled = (target/maxU)*(maxLinOut);



%% linear controller
policy_lin.fcn      = @(policy_lin,m,s)my_conlin(policy_lin,m,s);

policy_lin.p.w = 0;
policy_lin.p.b = target - maxU;

lin_out = policy_lin.fcn(policy_lin,x,0);

figure; plot(x,lin_out); grid on;

%% Saturating linear controller
policy_linSat  = policy_lin;
policy_linSat.fcn   = @(policy_linSat,m,s)my_conCat(@my_conlin,@my_gSat,policy_linSat,m,s);  % linear saturating controller
policy_linSat.maxU = maxU;
policy_linSat.lim  = [amin, amax];


for i=1:length(x)
    linSat_out(i) = policy_linSat.fcn(policy_linSat,x(i),0);
end

hold on;    plot(x,linSat_out);