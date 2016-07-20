% Fit hyperparameters:
fprintf('\nFitting policy hyperparameters to initials targets...\n');
clear inputs
noOfTr = 10;
% inputs = randn(length(poli),noOfTr);
% inputs(end,:) = rand(1,noOfTr)*10;

step = floor(H/noOfTr);
for i=1:noOfTr
    inputs(:,i) = [x(i+step*(i-1),dyno(poli))]';%, x(i+step*(i-1),end-1:end)]';
end

[newP, ~, ~] = fminsearch(@(p) initPolicy(p, policy, inputs, targets), unwrap(policy.p));
policy.p = rewrap(policy.p, newP);

% [newP, ~, ~] = minimize('initPolicy',(p, policy, inputs, targets), unwrap(policy.p));


% D = length(poli);
% E = length(policy.maxU);
% policy.p.hyp = zeros(D+2,E);
% policy.p.hyp(D+1,:) = log(1);                        % signal std dev
% policy.p.hyp(D+2,:) = log(0.01);                     % noise std dev (SNR = 100)

noOfQ = 20;
pout = zeros(length(policy.maxU),noOfQ);
if plotting.verbosity > 1
    for i=1:noOfQ
        input = 4*rand(length(poli),1)-2; input(end) = rand*10;
        pout(:,i) = policy.fcn(policy,input,zeros(7));
    end
    
    if ~ishandle(8)
        figure(8)
    else
        set(0,'CurrentFigure',8);
    end
    clf(8);
    grid on; plot(pout','x');
    hold on; line(get(gca,'XLim'),[targets(1) targets(1)],'Color',[0.5 0.5 0.5],'LineStyle',':');
    hold on; line(get(gca,'XLim'),[targets(2) targets(2)],'Color',[0.5 0.5 0.5],'LineStyle',':');
    title('Policy initialization')
    xlabel('Query no.');    ylabel('Computed action');
end
drawnow;
