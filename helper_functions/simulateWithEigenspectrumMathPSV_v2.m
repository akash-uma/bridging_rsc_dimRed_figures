function [ all_stats, all_correlations ] = simulateWithEigenspectrumMathPSV_v2(...
    basis, all_eigenspectrums, psv_coeffs, varargin)
    

    all_stats = cell(1, length(psv_coeffs));
    all_correlations = cell(size(all_eigenspectrums, 1), length(psv_coeffs));
    nVarargs = length(varargin);
    if (nVarargs == 1)
        equal_variances = diag(varargin{1});
    else
        equal_variances = [];
    end
    
    neurons_num = size(basis{1, 1}, 1);
    correlations_num = (neurons_num*neurons_num - neurons_num)/2;
    for p = 1 : length(psv_coeffs)
        for e = 1 : size(all_eigenspectrums, 1)
            cur_correlations = nan(size(basis, 1), correlations_num);
            all_correlations{e, p} = cur_correlations;
        end
    end
    
    % Will need to compute the correlations to plot the histograms
    for p = 1 : length(psv_coeffs)
        psv_value = psv_coeffs(p);
        % Given a particular eigenspectrum, use all the different basis
        % with that given eigenspectrum to generate a shared covariance
        % matrix. 
        espectrum_v_basis_stats_means = nan(size(all_eigenspectrums,1), ...
                                            size(basis, 1));
        espectrum_v_basis_stats_stds = nan(size(all_eigenspectrums,1), ...
                                           size(basis, 1));
                                       
        for b = 1 : size(basis, 1)
            evectors_basis = basis{b, 1};
            if isempty(equal_variances)
                p_variances = diag(basis{b, 2});
            else
                p_variances = equal_variances;
            end
            for e = 1 : size(all_eigenspectrums, 1)
                
                eigenspectrum = all_eigenspectrums(e, :);
                
                [svariance_factor, ~, ~] = adjustSharedVariances(...
                    evectors_basis, eigenspectrum, psv_value, p_variances, 0);
                
                evalues_dim = diag(svariance_factor*eigenspectrum);
                shared_covariance = evectors_basis*evalues_dim*evectors_basis';
                C = shared_covariance + p_variances;
                correlations = corrcov(C);
                correlations = correlations(triu(true(size(C)),1))';
                
                espectrum_psv_correlations = all_correlations{e, p};
                espectrum_psv_correlations(b, :) = correlations;
                all_correlations{e, p} = espectrum_psv_correlations;
                
                espectrum_v_basis_stats_means(e, b) = mean(correlations);
                espectrum_v_basis_stats_stds(e, b) = std(correlations);
            end
        end
        psv_stats = cell(2, 1);
        psv_stats{1, 1} = espectrum_v_basis_stats_means;
        psv_stats{2, 1} = espectrum_v_basis_stats_stds;
        all_stats{1, p} = psv_stats;
    end
end

