function [ all_basis_stats, all_correlations, all_sharedvar] = simulateOneVectorBases(one_vector_bases, ... 
        eigenvalues_coeff, psv_coeffs, p_variances)
    % Sweep through eigenvectors and psv values to generate different 
    % covariance matrices and observe their mean and std. 
    % We will use the same eigenspectrums and same psv values for all of 
    % the different bases. The idea is to check how the different 
    % eigenvectors change the mean and std even when the eigenspectrum and
    % psv is constant.
    
    neurons_num = size(one_vector_bases{1, 1}, 1);
    correlations_num = (neurons_num*neurons_num - neurons_num)/2;
    
    % we only have one eigenvector, so the eigenvalue_coeff must be a
    % scalar.
    all_basis_stats = cell(1, length(psv_coeffs));
    all_correlations = cell(1, length(psv_coeffs));
    all_sharedvar = cell(1, length(psv_coeffs));
    for p = 1 : length(psv_coeffs)
        % mean, std, psv, dshared are reported, 4 columns.
        stats_for_psv = nan(size(one_vector_bases, 1), 4);
        correlations_for_psv = nan(size(one_vector_bases, 1), correlations_num);
        sharedvar_for_psv = nan(size(one_vector_bases, 1), neurons_num);
        % the psv for each individual neuron
        all_basis_stats{1, p} = stats_for_psv;
        all_correlations{1, p} = correlations_for_psv;
        all_sharedvar{1,p} = sharedvar_for_psv;
    end

    for b = 1 : size(one_vector_bases, 1)
        cur_basis = one_vector_bases{b, 1};
        
        [cur_stats, cur_correlations, cur_sharedvar] = simulateWithEvectors(...
        cur_basis, eigenvalues_coeff, psv_coeffs, p_variances);
        
        for p = 1 : length(psv_coeffs)
            stats_for_psv = all_basis_stats{1, p};
            cur_stats_for_psv = cur_stats{1, p};
            stats_for_psv(b, :) = cur_stats_for_psv(:);
            all_basis_stats{1, p} = stats_for_psv;
            
            correlations_for_psv = all_correlations{1, p};
            cur_correlations_for_psv = cur_correlations{1, p};
            cur_sharedvar_for_psv = cur_sharedvar{1,p};
            correlations_for_psv(b, :) = cur_correlations_for_psv(:);
            sharedvar_for_psv(b,:) = cur_sharedvar_for_psv;
            all_correlations{1, p} = correlations_for_psv;
            all_sharedvar{1,p} = sharedvar_for_psv;
        end
    end
end

