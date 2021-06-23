function [ num_needed ] = calculateDshared(eigenvalues)
    total_variance = sum(eigenvalues);
    num_needed = length(eigenvalues);
    for i = 1:length(eigenvalues)
        num_needed = i;
        cur_variance = sum(eigenvalues(1:i));
        if (cur_variance/total_variance >= 0.95)
            return
        end
    end
end

