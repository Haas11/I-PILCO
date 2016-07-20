function f = initPolicy(p, policy, inputs, targets)
%initPolicy Function that computes the the difference between the output of
% a given policy and a set of targets. Used for initializing the
% hyperparemeters of a policy.
%
% Input arguments:
% p: parameters of the policy
% policy: structure
% inputs: query points to use for fit
% targets: vector of initialization targets
%
% Outputs:
% f : function value, sum of differences between targets.

noOfQuery = size(inputs,2);
policy.p = rewrap(policy.p, p);
difference = zeros(1,noOfQuery);

for i=1:noOfQuery
    m = inputs(:,i);    
    [actions, ~, ~, ~] = policy.fcn(policy, m, zeros(length(m)));    
    difference(1,i) = norm(actions - targets', 2);
end

f = sum(difference);

%df = dMdm;      % gradient of the difference ( = grad of policy_out)

end

