function t = minimum_euclidean_distance(frequency, measured_reflectivity, ks, thickness_step, E_oil, E_air, temp, salinity, theta)
    thickness = 0:thickness_step:10;
    
    theoretical_reflectivity = 10*log10(abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta))); 
    
    distance = sqrt(sum((theoretical_reflectivity-transpose(measured_reflectivity)).^2, 1));
    [~, indexOfMinimum] = min(distance);
    t = thickness(indexOfMinimum);
end