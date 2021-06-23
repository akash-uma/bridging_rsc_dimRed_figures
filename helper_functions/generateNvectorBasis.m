function [ all_new_bases ] = generateNvectorBasis(one_vector_basis, n, total_bases)
    bases_choices = size(one_vector_basis, 1);
    neurons_num = size(one_vector_basis{1, 1}, 1);
    all_new_bases = cell(1, 1);
    k = 1;
    while k <= total_bases
        ind = randi([1, bases_choices], 1, n);
        % combine the n 1-vector basis into 1 basis with n orthonormal
        % vectors. This creates an orthonormal basis with n vectors.
        new_basis = nan(neurons_num, n);
        for i = 1 : length(ind)
            index = ind(i);
            cur_basis = one_vector_basis{index, 1};
            new_basis(:, i) = cur_basis;
        end
        if rank(new_basis) == n
            orthonormal_basis = gramSchmidt( new_basis );
            all_new_bases{k, 1} = orthonormal_basis;
            k = k + 1;
        end
    end

end

