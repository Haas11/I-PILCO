% Print relevant tuning parameters to console.

fprintf('\n\n')

fprintf('GP Model parameters:\n')
fprintf('=====================================================\n');
fprintf('Outputs: \t\t\t')
disp(stateNames(indices(dyno)))
fprintf('Inputs: \t\t\t')
disp(stateNames(indices(dyno(dyni))))
fprintf('Outputs that are differences:')
disp(stateNames(indices(dyno(difi))))
fprintf('\nOutputs w/ prior reference:')
disp(stateNames(indices(dyno(refi))))
fprintf('No. of inducing inputs: \t%i \n', nii)
fprintf('Sampling interval: \t\t%4.2f [s] \n', dt_pilco)
fprintf('Noisy inputs: \t\t\t%i \n', noisyInputs)

fprintf('\n\nCost parameters:\n')
fprintf('=====================================================\n');
fprintf('Energy penalty: \t%4.5f \n', cost.ep)
fprintf('Exploration param: \t%4.5f \n\n', cost.expl)
if isfield(cost,'sub')
    for i=1:length(cost.sub)
        fprintf('Inputs %i: \t',i); disp(stateNames(indices(dyno(cost.sub{i}.losi))))
        fprintf('Cost width %i: \t\t%4.5f \n\n', i, cost.sub{i}.width)
    end
else
    fprintf('Inputs %i: \t',i); disp(stateNames(indices(dyno(cost.losi))))
    fprintf('Cost width \t\t%4.2f \n\n', cost.width)
end

fprintf('\n\nPolicy parameters:\n')
fprintf('=====================================================\n');
if isfield(policy.p,'w')
    fprintf('linear controller');
else
    fprintf('RBF controller\n');
    fprintf('No. of policy kernels: \t\t\t%i \n\n', nc)
end
fprintf('Learned parameters: \t\t');
disp(actionTitles)
fprintf('Inputs: \t\t\t')
disp(stateNames(indices(dyno(poli))))

fprintf('No. of learning iterations: \t %i \n', N)
fprintf('No. of linesearches: \t\t\t%i \n', opt.length)
fprintf('Function evaluations per linesearch: \t%i \n', opt.MFEPLS)
