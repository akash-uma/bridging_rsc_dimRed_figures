function load_sim = computeLoadSim( x )
    
    x = x./norm(x);
    d = length(x);

    load_sim = d*(mean(x)^2);

end

