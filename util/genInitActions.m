function a_init = genInitActions(policy, J, type, actionTitles, t, varargin)
%genInitActions generate distributions for actions in initial rollouts.
%
%   three types available: 1=gaussian, 2=Uniform, 3=Random Walk (Brownian)


%% Code
a_init = cell(1,J);
H = length(t);

% Normally distributed initial rollouts:
if type==1
    initMean = policy.maxU;     initVar = policy.maxU.*2;
    for i=1:J
        a_init{i} = gaussian(initMean,diag(initVar),H)';
    end
    
elseif type==2
    % Uniformly distributed initial rollouts:
    for i=1:J
        a_init{i} = repmat(policy.minU,H,1) + repmat(policy.maxU.*2-policy.minU,H,1).*rand(H,length(policy.impIdx));
    end
    
elseif type==3
    sig_ratio = varargin{1};
    
    % Ohrnstein-Uhlenbeek Stochastic Process initial rollouts:
    ou_opts = sdeset('RandSeed',2);
    th          = [0.25; 0.3];   % drift rate parameters
%     sig         = [150; 150];     % diffusion parameters
    
    sig = repmat(policy.maxU(1)*2/sig_ratio,J,1);
    
    startPoint  = [policy.minU(policy.impIdx);
        policy.maxU(policy.impIdx).*2];
    mu  = [policy.maxU(policy.impIdx).*2;
        policy.minU(policy.impIdx)];
    
    for i=1:J
        a_init{i} = min(policy.maxU(1)*2,max(policy.minU(1),sde_ou(th(i),mu(i,:),sig(i),t,startPoint(i,:),ou_opts)));
    end

end

if ~ishandle(12)
    figure(12); 
else
    set(0,'CurrentFigure',12); 
end
clf(12);
for i=1:J
    stairs((1:length(a_init{i})).*(t(2)-t(1)),a_init{i}); hold on;
end
grid on; axis tight;
title('Initial Rollout Actions')
xlabel('Time steps'); ylabel('Actions');
legendVec = cell(1,J*length(policy.maxU(policy.impIdx)));
for j=1:J
    for u=1:length(policy.maxU(policy.impIdx))
        temp{u+2*(j-1)} = strcat(actionTitles(u),' (J=',num2str(j),')');
        legendVec{u+2*(j-1)} = temp{u+2*(j-1)}{1};
    end
end
legend(legendVec);


end

