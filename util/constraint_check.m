function flag = constraint_check(dq, F0, F1, t)
% summary: check if any constraints are violated by the simualtion.
% global ROBOT
flag = 0;

if max(F1) > 500
    flag = 1;
    fprintf('\nsimulation stopped due to high interaction forces.\n');
elseif max(dq) > 50
    fprintf('\nsimulation stopped. maximum joint speed reached.\n');
    flag = 1;
elseif t > 5 && (F1(1) - F0(1) > 15)
    flag = 1;
    fprintf('\nsimulation stopped due to force oscillations.\n');
end

% Mt = ROBOT.maniplty(q','dof',[1 1 0 0 0 0])
% Mr = ROBOT.maniplty(q','dof',[0 0 0 0 0 1])
%
% if M < 0.1
%     flag = 1;
%     fprintf('Simulation ended due to low manipulability');
% end