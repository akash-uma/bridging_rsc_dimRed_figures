neurons_num = 30;

%% Generate many co-fluctuation patterns with different loading similarities

% Store the new basis
rng(0);
mu = 2.5; 
stds = 0.1:0.1:5.5;
trials_per_gaussian = 50; 
one_vector_bases = generateGaussianOneVectorBases(neurons_num, mu, stds, trials_per_gaussian);

variance_all = nan(length(one_vector_bases),1);
loading_similarity_all = nan(length(one_vector_bases),1);

% Compute loading similarity of all vectors
for b = 1 : length(one_vector_bases)
    vector = one_vector_bases{b};

    variance_all(b) = var(vector);
    loading_similarity_all(b) = (1/(neurons_num-1) - var(vector))*(neurons_num-1);
end

[sorted_variance, variance_idxs] = sort(variance_all);
one_vector_bases_variance = one_vector_bases(variance_idxs);
loading_similarity_all = loading_similarity_all(variance_idxs);

%% Simulate using the created basis
eigenvalues_coeffs = 1;
psv_coeffs = [0.3, 0.5];
p_variances = ones(1, neurons_num); % all the same

[all_basis_stats_variance, all_correlations_variance] = ...
    simulateOneVectorBases(one_vector_bases_variance, ... 
        eigenvalues_coeffs, psv_coeffs, p_variances);

%% Color based on the variance of the eigenvector components similarity
percent_sv_idx = 2;

target_stats_variance = all_basis_stats_variance{1, percent_sv_idx};
target_means_variance = target_stats_variance(:, 1);
target_stds_variance = target_stats_variance(:, 2);
target_correlations_variance_5 = all_correlations_variance{1, percent_sv_idx};

% set color to be based on variance
rounded_sorted_variance = round(sorted_variance, 6);
[values, IA, IC] = unique(rounded_sorted_variance);
total_points = length(values);
cmap = colormap(winter(total_points));
colors = cmap(IC, :);
colors = flipud(colors);

% do a random shuffle, so that there is no bias to have the points with
% high variance on top of those with low variance. 
random_indices = randperm(length(one_vector_bases));
target_means_variance = target_means_variance(random_indices);
target_stds_variance = target_stds_variance(random_indices);
colors = colors(random_indices, :);
target_correlations_variance_5 = target_correlations_variance_5(random_indices, :);

figure(1);
scatter(target_means_variance, target_stds_variance, 60, colors, 'filled');
hold all;

%% inner arc (%sv=30%)
percent_sv_idx = 1;

target_stats_variance = all_basis_stats_variance{1, percent_sv_idx};
target_means_variance = target_stats_variance(:, 1);
target_stds_variance = target_stats_variance(:, 2);
target_correlations_variance_3 = all_correlations_variance{1, percent_sv_idx};

target_correlations_variance_3 = target_correlations_variance_3(random_indices, :);
target_means_variance = target_means_variance(random_indices);
target_stds_variance = target_stds_variance(random_indices);

hold all;

scatter(target_means_variance, target_stds_variance, 60, colors, 'filled');
colorbar
radii = 0.1:0.1:0.3;
hold on;
for r = 1 :length(radii)
    radius = radii(r);
    fnCircle(0,0, radius);
end

xlabel('r_{sc} mean'); ylabel('r_{sc} s.d.');
title('Figure 3f')
box off; axis tight; axis equal;
axis([0 0.51 0 0.51]); 
set(gca,'fontsize', 18, 'linewidth', 1.5);
set(gca,'XTick',0:.1:.5,'YTick',0:.1:.5);
