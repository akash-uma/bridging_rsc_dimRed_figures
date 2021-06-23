neurons_num = 30;
dim_num = 2;
figure(3); pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 1]);

%% Generate one-vector basis
rng(10);
mu = 2.5; 
stds = 0.1:0.1:5.5;
trials_per_gaussian = 50;
one_vector_basis = generateGaussianOneVectorBases(neurons_num, mu, stds, trials_per_gaussian);
 
%% Analysis with dshared = 2
total_bases = 3000;
multi_vector_basis = generateNvectorBasis(one_vector_basis, dim_num, total_bases);
 
%% Compute loading similarity of the first and second eigenvector
[first_load_sims,second_load_sims] = deal(nan(1,length(multi_vector_basis)));
for i = 1 : length(multi_vector_basis)
    basis = multi_vector_basis{i};
    first_load_sims(i) = computeLoadSim(basis(:,1));
    second_load_sims(i) = computeLoadSim(basis(:,2));
end
 
%% Simulate with the 2-vector basis
psv_coeffs = 0.5;
p_variances = ones(1, neurons_num);

training_eigenspectrums = [95 5; 100 25; 100 100; 25 100; 5 95];
espectrum_sum = sum(training_eigenspectrums, 2);
espectrum_sum = repmat(espectrum_sum, 1, 2);
training_eigenspectrums = training_eigenspectrums ./ espectrum_sum;
% espectrums_dshared = calculateEspectrumsDshared(training_eigenspectrums);

[all_basis_stats, all_correlations] = simulateWithEigenspectrumMathPSV_v2(multi_vector_basis, ...
    training_eigenspectrums, psv_coeffs, p_variances);

%% Plot the results
all_basis_stats = all_basis_stats{1, 1};
all_means = all_basis_stats{1, 1};
all_stds = all_basis_stats{2, 1};
radii = [0.1 0.2 0.3];

first_high_ls = first_load_sims>.8;
first_low_ls = first_load_sims<.2;
second_high_ls = second_load_sims>.8;
second_low_ls = second_load_sims<.2;

%% plot for high 1st ls & low 2nd ls
subplot(1,2,1); hold on;
cols = [255 0 0; 161 33 33; 128 128 128; 76 76 76; 25 25 25]./255;
for e = 1:size(training_eigenspectrums,1)
    cur_means = all_means(e,:);
    cur_stds = all_stds(e,:);
    
    cur_idx = first_high_ls & second_low_ls;
    
    cur_first_ls = first_load_sims(cur_idx);
    cur_second_ls = second_load_sims(cur_idx);
    
    n = sum(cur_idx);
    
    cur_means = cur_means(cur_idx);
    cur_stds = cur_stds(cur_idx);
    
    rng(0);
    rng_idx = randperm(length(cur_means));
    cur_means = cur_means(rng_idx);
    cur_stds = cur_stds(rng_idx);
    
    cur_means = cur_means(1:min([200, n]));
    cur_stds = cur_stds(1:min([200, n]));
    
    scatter(cur_means,cur_stds,40,cols(e,:),'filled');
end

for r = 1 :length(radii)
    radius = radii(r);
    fnCircle(0,0, radius);
    hold on;
end
xlabel('r_{sc} mean'); ylabel('r_{sc} s.d.');
title('Figure 4a');
box off; axis tight; axis equal;
axis([-0.02 0.51 0 0.51]);
set(gca,'fontsize', 18, 'linewidth', 1.5);
legend('95/5','80/20','50/50','20/80','5/95','Location','Best'); legend boxoff;

%% plot for low 1st ls & low 2nd ls
subplot(1,2,2); hold on;
cols = [128 128 128; 76 76 76; 25 25 25]./255;
for e = 1:3
    cur_means = all_means(e,:);
    cur_stds = all_stds(e,:);
    
    cur_idx = first_low_ls & second_low_ls;
    
    cur_first_ls = first_load_sims(cur_idx);
    cur_second_ls = second_load_sims(cur_idx);
    
    n = sum(cur_idx);
    
    cur_means = cur_means(cur_idx);
    cur_stds = cur_stds(cur_idx);
    
    cur_means = cur_means(1:min([200, n]));
    cur_stds = cur_stds(1:min([200, n]));
    
    scatter(cur_means,cur_stds,40,cols(e,:),'filled');
end

for r = 1 :length(radii)
    radius = radii(r);
    fnCircle(0,0, radius);
    hold on;
end
xlabel('r_{sc} mean'); ylabel('r_{sc} s.d.');
title('Figure 4b');
box off; axis tight; axis equal;
axis([-0.02 0.51 0 0.51]);
set(gca,'fontsize', 18, 'linewidth', 1.5);

legend('95/5','80/20','50/50','Location','Best'); legend boxoff;
