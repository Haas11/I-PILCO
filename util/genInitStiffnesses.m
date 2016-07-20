% Ohrnstein-Uhlenbeek Stochastic Process for initial rollouts:
ou_opts = sdeset('RandSeed',2);
K0 = cell(1,J);
th = [0.3; 0.2; 0.2];   % drift rate parameters
sig = [20; 12; 10];
startPoint = [20 20 20;
              40 40 40;
              60 125 125];
        mu = [150 150 150;
              20 75 75;
              40 50 50];
for i=1:J
    K0{i}(:,1) = t_pilco;
    K0{i}(:,2:1+length(policy.impIdx)) = min(policy.maxU(1)*2,max(policy.minU(1),sde_ou(th(i),mu(i,:),sig(i),t_pilco,startPoint(i,:),ou_opts)));
end