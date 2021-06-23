evector_num = 5; % number of eigenvectors in the basis
neurons_num = 30;

%% Generate one-vector basis
rng(0);
basis = cell(1, 1); %just initialize a cell
mu = 2.5; 
stds = 0.1:0.1:5.5;
trials_per_gaussian = 50;
one_vector_basis = generateGaussianOneVectorBases(neurons_num, mu, stds, trials_per_gaussian);

%% N-vector basis analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate basis with N eigenvectors each
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%e%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_bases = 1000;
n_vector_basis = generateNvectorBasis(one_vector_basis, evector_num, total_bases);

%% Stimulate with the n vector basis

% Now we have a basis with the desired number of eigenvectors. 
% Sweep through eigenvectors and psv values to generate different 
% covariance matrices and observe their mean and std. 
% We will use the same eigenspectrums and same psv values for all of the
% different basis. The idea is to check how the different eigenvectors
% change the mean and std even when the eigenspectrum and psv is constant.
eigenvalues_coeffs = [0];
psv_coeffs = [0.5 0.3];

all_basis_stats = cell(length(eigenvalues_coeffs), length(psv_coeffs));
all_correlations = cell(length(eigenvalues_coeffs), length(psv_coeffs));
neurons_num = 30;
p_variances = ones(1, neurons_num).*1; % all the same
for e = 1 : length(eigenvalues_coeffs)
    evalue_coeff = eigenvalues_coeffs(e);
    for p = 1 : length(psv_coeffs)
        psv_value = psv_coeffs(p);
        espectrum_means = nan(evector_num, total_bases);
        espectrum_stds = nan(evector_num, total_bases) ;
        espectrum_psvs = nan(evector_num, total_bases) ;
        espectrum_dshare = nan(evector_num, total_bases);
        espectrum_adjustments = nan(evector_num, total_bases);
        eval_psv_correlations = cell(total_bases, 1);
        for b = 1 : total_bases
%             fprintf('b = %d, p = %d, e= %d\n', b, p, e);
            evectors_base = n_vector_basis{b, 1};
            % using the given eigenvectors, psv value and eigenspectrum
            % compute statistics related to the covariance matrix simulated
            % from the FA model using the given parameters. These include
            % the rsc mean, std, %sv, dshared, and an adjusting factors for
            % the eigenspectrum used to keep the %sv at the desired level.
            [basis_stats, basis_correlations] = simulateWithEvectors(...
                evectors_base, evalue_coeff, psv_value, p_variances);
            
            basis_stats = basis_stats{1, 1};
            basis_correlations = basis_correlations{1, 1};
            espectrum_means(:, b) =  basis_stats(:, 1);
            espectrum_stds(:, b) =  basis_stats(:, 2);
            espectrum_psvs(:, b) =  basis_stats(:, 3);
            espectrum_dshare(:, b) =  basis_stats(:, 4);
            eval_psv_correlations{b, 1} = basis_correlations;
        end
        espectrum_stats = cell(5, 1);
        espectrum_stats{1, 1} = espectrum_means;
        espectrum_stats{2, 1} = espectrum_stds;
        espectrum_stats{3, 1} = espectrum_psvs;
        espectrum_stats{4, 1} = espectrum_dshare;
        all_basis_stats{e, p} = espectrum_stats;
        all_correlations{e, p} = eval_psv_correlations;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmap = colormap;
colors_raw_pink_3 = [128,0,128; 197,27,170; 247,104,161; 250,159,181; 255,205,118];
figure(2);
colors = colors_raw_pink_3 ./ 255;

ax = [];
for dshare_val = 1 : evector_num
    color = colors(dshare_val, :);
    for e = 1 : 1
        evalue = eigenvalues_coeffs(e);
        for dim = 1 : evector_num
            for p = 1 : 1
                psv_value = psv_coeffs(p);
                cur_marker = 'o';
                espectrum_stats = all_basis_stats{e, p};
                espectrum_means = espectrum_stats{1, 1};
                espectrum_stds = espectrum_stats{2, 1};
                espectrum_dshare = espectrum_stats{4, 1};
                dim_means = espectrum_means(dim, :);
                dim_stds = espectrum_stds(dim, :);
                dim_dshared = espectrum_dshare(dim, :);
                if dim_dshared(1) == dshare_val
                    ax = [ax, scatter(dim_means, dim_stds, 60, color, cur_marker, 'filled')];
                    hold on;
                end
            end
            hold on;
        end
    end
end

hold all;

radii = 0.1:0.1:0.3;
hold on;
for r = 1 :length(radii)
    radius = radii(r);
    fnCircle(0,0, radius);
    hold on;
end
legend(ax,{'1-d','2-d','3-d','4-d','5-d'},'Location','Best'); legend boxoff;
xlabel('r_{sc} mean'); ylabel('r_{sc} s.d.');
title('Figure 3g')
box off; axis tight; axis equal;
axis([0 0.51 0 0.51]);
set(gca,'fontsize', 18, 'linewidth', 1.5);
set(gca,'XTick',0:.1:.5,'YTick',0:.1:.5);



