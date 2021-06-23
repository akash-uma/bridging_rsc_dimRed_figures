function [ svariance_factor, mid_point, error ] = adjustSharedVariances(...
        evectors, evalues, psv_value, private_variances, start_index)

    range = 0.01:0.001:5000;
    lower_bound = 1;
    upper_bound = length(range);
    if (start_index == 0)
        mid_point = floor((lower_bound + upper_bound) / 2);
    else
        mid_point = start_index;
    end
    
    original_evalues = evalues;
    
    svariance_factor = range(mid_point);
    evalues = diag(svariance_factor*original_evalues);
    shared_covariance = evectors*evalues*evectors';
    shared_variance = computePercentSharedVariance(shared_covariance, private_variances);
    
    epsilon = 0.001;
    % binary search
    while (upper_bound > lower_bound)
        if abs(shared_variance - psv_value) <= epsilon
            error = 0;
            return;
        end
        if (shared_variance > psv_value)
            upper_bound = mid_point;
        else
            lower_bound = mid_point + 1;
        end
        mid_point = floor((upper_bound + lower_bound) / 2);
        svariance_factor = range(mid_point);
        evalues = diag(svariance_factor*original_evalues);
        shared_covariance = evectors*evalues*evectors';
        shared_variance = computePercentSharedVariance(shared_covariance, ...
            private_variances);
    end
    
    error = 1;
    disp('Did NOT Converge!!');
end

