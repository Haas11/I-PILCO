function a_init = genInitActions(policy, J, type, actionTitles, varargin)
%genInitActions generate distributions for actions in initial rollouts.
%
%   1=gaussian, 2=Uniform, 3=Random Walk (Brownian), 4=sinusoid


%% Code
a_init = cell(1,J);

% Normally distributed initial rollouts:
if type==1
    t = varargin{1};
    H = varargin{2};
    initMean = policy.maxU;     initVar = policy.maxU.*2;
    for i=1:J
        a_init{i} = [t(1:H),gaussian(initMean,diag(initVar),H)'];
    end
    
elseif type==2
    % Uniformly distributed initial rollouts:
    t = varargin{1};
    H = varargin{2};
    for i=1:J
        a_init{i} = [t(1:H), repmat(policy.minU(policy.impIdx),H,1) + repmat(policy.maxU(policy.impIdx).*2-policy.minU(policy.impIdx),H,1).*rand(H,2)];
    end
    
elseif type==3
    % Ohrnstein-Uhlenbeek Stochastic Process initial rollouts:
    ou_opts = sdeset('RandSeed',2);
    
    t = varargin{1};
    if nargin > 5
        sig_ratio = varargin{2};
        if nargin > 6
            th = ones(1,J)*varargin{3};         % convergence speed
        else
            th = ones(1,J)*0.3;
        end
    else
        sig_ratio = 5;                          % diffusion ratio
        th = ones(1,J)*0.3;
    end
    sig = repmat(policy.maxU(1)*2/sig_ratio(1),J,1);    % spread
    
    startPoint  = [policy.minU(policy.impIdx);
        policy.maxU(policy.impIdx).*2;
        policy.maxU(policy.impIdx)*1.25];
    mu  = [policy.maxU(policy.impIdx).*2;
        policy.minU(policy.impIdx);
        policy.maxU(policy.impIdx)*1.25];
    
    for i=1:J
        a_init{i} = [t, min(policy.maxU(1)*2,max(policy.minU(1),sde_ou(th(i),mu(i,:),sig(i),t,startPoint(i,:),ou_opts)))];
    end
    
elseif type==4
    t = varargin{1};
    f = varargin{2};
    jitter = varargin{3};
    for i=1:J
        a_init{i} = [t, repmat(policy.maxU(policy.impIdx),length(t),1).*[sin(f*t + rand), cos(f*t + rand)]];
        a_init{i}(:,2:end) = a_init{i}(:,2:end) + repmat(policy.maxU(policy.impIdx),length(t),1) + rand(length(t),size(a_init{i},2)-1).*repmat(policy.maxU(policy.impIdx)./jitter,length(t),1);
    end
end

%% Plot Results:
if ~ishandle(22)
    figure(22);
else
    set(0,'CurrentFigure',22);
end
clf(22);
for i=1:J
    stairs(a_init{i}(:,2:end),'Linewidt',1.5); hold on;
end
axis tight
grid minor;
title('\fontsize{16}Initial Stiffness Trajectories','Interpreter','Tex');
xlabel('\fontsize{16}Time steps [d_t = 0.1]','Interpreter','Tex'); 
ylabel('\fontsize{16} K_P [N/m]','Interpreter','Tex');
xt = get(gca, 'XTick');
set(xt, 'FontSize', 16)
legendVec = cell(1,J*length(policy.maxU(policy.impIdx)));
for j=1:J
    for u=1:length(policy.maxU(policy.impIdx))
        temp{u+2*(j-1)} = strcat(actionTitles(u),' (J=',num2str(j),')');
        legendVec{u+2*(j-1)} = temp{u+2*(j-1)}{1};
    end
end
legend(legendVec);
end