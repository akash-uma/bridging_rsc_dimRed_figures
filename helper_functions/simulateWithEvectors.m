function [all_stats, all_correlations, all_sharedvar] = simulateWithEvectors(...
        evectors_basis, eigenvalues_coeffs, psv_coeffs, p_variances)
    
    [neurons_num, total_dims] = size(evectors_basis);
    correlations_num = (neurons_num*neurons_num - neurons_num)/2;
    
    % private variances remain constant throughout
    p_variances = diag(p_variances);
    dimensions = 1:total_dims;
    
    all_stats = cell(length(eigenvalues_coeffs), length(psv_coeffs));
    all_correlations = cell(length(eigenvalues_coeffs), length(psv_coeffs));
    all_sharedvar = cell(length(eigenvalues_coeffs), length(psv_coeffs));
    
    for k = 1 : length(eigenvalues_coeffs)
        if (eigenvalues_coeffs(k) == 0)
            all_evalues = ones(size(dimensions));
        else
            all_evalues = exp(-dimensions/eigenvalues_coeffs(k));
        end
        for p = 1 : length(psv_coeffs)
            psv_value = psv_coeffs(p);

            % storing the mean, std, psv, dshared and individual psvs
            dim_v_stats = nan(total_dims, 4);

            % save the correlations so we can plot the histograms to check the
            % distributions

            dim_saved_correlations = nan(total_dims, correlations_num);
            dim_saved_sharedvar = nan(total_dims,neurons_num);
            bs_index = 0;
            for dim = 1 : total_dims
                evectors_dim = evectors_basis(:, 1:dim);
                evalues_dim = all_evalues(1:dim);
                evalues_sum = sum(evalues_dim);
                evalues_dim = 20*(evalues_dim/evalues_sum);
                % adjust the eigenvalues to achieve the desired psv
                % binary search index - the index for the initial query
                [svariance_factor, bs_index, error] = adjustSharedVariances(...
                    evectors_dim, evalues_dim, psv_value, p_variances, bs_index);

                text = sprintf('dim = %d svariance_factor=%f \n', dim, ...
                    svariance_factor);
                if (error == 1)
                    disp('ERROR');
                    fprintf(text);
                end
                evalues_dim = diag(svariance_factor*evalues_dim);
                shared_covariance = evectors_dim*evalues_dim*evectors_dim';
                dim_saved_sharedvar(dim,:) = diag(shared_covariance);
                C = shared_covariance + p_variances;
                correlations = corrcov(C);
                correlations = correlations(triu(true(size(C)),1))';
                dim_saved_correlations(dim, :) = correlations;
                dim_v_stats(dim, 1) = mean(correlations);
                dim_v_stats(dim, 2) = std(correlations,1);
                dim_v_stats(dim, 3) = computePercentSharedVariance(...
                            shared_covariance, p_variances);
                dim_v_stats(dim, 4) = calculateDshared(diag(evalues_dim));
            end

            all_stats{k, p} = dim_v_stats;
            all_correlations{k, p} = dim_saved_correlations;
            all_sharedvar{k,p} = dim_saved_sharedvar;
        end
    end
end

