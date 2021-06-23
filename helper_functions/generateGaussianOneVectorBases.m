function [ basis ] = generateGaussianOneVectorBases(neurons_num, mu, stds, trials_per_gaussian)
    % store the new basis
    basis = cell(1, 1); %just initialize a cell
    k = 1;
    for s = 1 : length(stds)
        cur_std = stds(s);
        for t = 1 : trials_per_gaussian
            evector = normrnd(mu, cur_std, neurons_num, 1);
            evector = evector/norm(evector);
            basis{k, 1} = evector;
            k = k + 1;
        end
    end
end

