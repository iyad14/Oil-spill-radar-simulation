function scatterplot(frequency, measured_reflectivity, ks, thickness_step, t, variance, trials, E_oil, E_air, temp, salinity, theta)
    clf;
    thickness = 0:thickness_step:10;
    
        %% Creating the theoretical curve by which the measured reflectivities are compared
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    theoretical_reflectivity = abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        %% Generating the random reflectivities with given frequencies and at a given thickness to test the frequencies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    noise = sqrt(variance)*randn(length(frequency), trials);
    [~, indexOfThickness] = min(thickness-t);
    reflectivities = 10*log10(abs(theoretical_reflectivity(:, indexOfThickness) + noise));
    hold on;
    scatter(reflectivities(1, :), reflectivities(2, :), 'x');
    scatter(10*log10(abs(theoretical_reflectivity(1, :))), 10*log10(abs(theoretical_reflectivity(2, :))), '.');
    scatter(measured_reflectivity(1), measured_reflectivity(2), 'x', 'black');
    
   
    
    position = minimum_euclidean_distance(frequency, measured_reflectivity, ks, thickness_step, E_oil, E_air, temp, salinity, theta);
    [~, indexOfThickness] = min(abs(thickness-position));
    scatter(10*log10(theoretical_reflectivity(1, indexOfThickness)), 10*log10(theoretical_reflectivity(2, indexOfThickness)));
     hold off;
    xlim([-8 0]);
    ylim([-8 0]);
end