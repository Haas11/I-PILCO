%% my_predcost.m
% *Summary:* Compute trajectory of expected costs for a given set of 
% state distributions
%
% inputs:
% m0          mean of states, D-by-1 or D-by-K for multiple means
% S           covariance matrix of state distributions
% dynmodel    (struct) for dynamics model (GP)
% plant	      (struct) of system parameters
% policy      (struct) for policy to be implemented
% cost        (struct) of cost function parameters
% H           length of optimization horizon
%
% outputs:
% L            expected cumulative (discounted) cost
% s            standard deviation of cost
% Mx           (cell) means of predicted state trajectories
% Sx           (cell) variances of predicted state trajectories
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2012-01-12
%
%% High-Level Steps
% # Predict successor state distribution
% # Predict corresponding cost distribution

function [L, s, Mx, Sx, Mcon, Scon] = my_predcost(m0, S, dynmodel, plant, policy, cost, H)
%% Code 
K = size(m0,2);
poli = plant.poli;

L = zeros(size(m0,2),H+1); 
s = zeros(size(m0,2),H+1);

Mx = cell(K,1); Mcon = Mx;
Sx = cell(K,1); Scon = Sx;

for k = 1:K
  m = m0(:,k);
  Mx{k}(:,1) = m;
  Sx{k}(:,:,1) = S;
  [ma, sa, ~] = policy.fcn(policy, m(poli), S(poli,poli));
  [L(k,1), ~, ~, v]  = cost.fcn(cost, m, S, ma, policy);
  s(k,1) = sqrt(v);
  Mcon{k}(:,1) = ma;
  Scon{k}(:,:,1) = sa;
  for t = 1:H
    [m, S] = plant.prop(m, S, plant, dynmodel, policy);	     % get next state
    
    [ma, sa, ~] = policy.fcn(policy, m(poli), S(poli,poli));
    [L(k,t+1), ~, ~, v]  = cost.fcn(cost, m, S, ma, policy);  
    s(k,t+1) = sqrt(v);
    
    Mx{k}(:,t+1) = m;
    Sx{k}(:,:,t+1) = S;
    
    Mcon{k}(:,t+1) = ma;
    Scon{k}(:,:,t+1) = sa;
  end
end
L = mean(L,1); s = mean(s,1); 