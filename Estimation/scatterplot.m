%%  frequency               --> frequencies used given (in GHz)
%%  measured_reflectivity   --> Reflectivitis measured by drone (in dB scale)
%%  ks                      -->  Surface roughness
%%  thickness_step          --> Thickness resolution (in mm)
%%  t                       --> Thickness at which the intended noisy reflectivities need to be generated
%%  variance                --> Noise variance
%%  trials                  --> Quantity of noisy reflectivity values needed to be generated
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)




function scatterplot(frequency, measured_reflectivity, ks, thickness_step, t, variance, trials, E_oil, E_air, temp, salinity, theta)
        
    clf;                                    % clear figures
    thickness = 0:thickness_step:10;        % thickness over which the reflectivities will be calculated
    
    
        %% Creating the theoretical curve by which the measured reflectivities are compared
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    theoretical_reflectivity = abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        %% Generating the random reflectivities with given frequencies and at a given thickness to test the frequencies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    noise = sqrt(variance)*randn(length(frequency), trials);
    [~, indexOfThickness] = min(thickness-t);
    noisy_reflectivities = 10*log10(abs(theoretical_reflectivity(:, indexOfThickness) + noise));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        %% Finding the point closest to the measured reflectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    position = minimum_euclidean_distance(frequency, measured_reflectivity, ks, thickness_step, E_oil, E_air, temp, salinity, theta);
    [~, indexOfThickness] = min(abs(thickness-position));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    
     %% Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if length(frequency) == 2
        
        % 2D Scatter plots 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
        hold on;
        scatter(noisy_reflectivities(1, :), noisy_reflectivities(2, :), 'x');
        scatter(10*log10(abs(theoretical_reflectivity(1, :))), 10*log10(abs(theoretical_reflectivity(2, :))), '.');
        scatter(measured_reflectivity(1), measured_reflectivity(2), 'x', 'black');
        scatter(10*log10(theoretical_reflectivity(1, indexOfThickness)), 10*log10(theoretical_reflectivity(2, indexOfThickness)));
        grid on;
        hold off;
        xlim([-8 0]);
        xlabel(strcat("Reflectivity(dB)@", num2str(frequency(1)), "GHz"));
        ylim([-8 0]);
        ylabel(strcat("Reflectivity(dB)@", num2str(frequency(2)), "GHz"));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           
    else      
        % 3D Scatter plots 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        scatter3(noisy_reflectivities(1, :), noisy_reflectivities(2, :), noisy_reflectivities(3, :), 'x');
        hold on;
        scatter3(10*log10(abs(theoretical_reflectivity(1, :))), 10*log10(abs(theoretical_reflectivity(2, :))), 10*log10(abs(theoretical_reflectivity(3, :))), '.');
        scatter3(measured_reflectivity(1), measured_reflectivity(2), measured_reflectivity(3), 'x', 'black');
        scatter3(10*log10(theoretical_reflectivity(1, indexOfThickness)), 10*log10(theoretical_reflectivity(2, indexOfThickness)), 10*log10(theoretical_reflectivity(3, indexOfThickness)));
        grid on;
        
        xlim([-8 0]);
        xlabel(strcat("Reflectivity(dB)@", num2str(frequency(1)), "GHz"));
        ylim([-8 0]);
        ylabel(strcat("Reflectivity(dB)@", num2str(frequency(2)), "GHz"));
        zlim([-8 0]);
        zlabel(strcat("Reflectivity(dB)@", num2str(frequency(3)), "GHz"));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
        
        %% Histograms 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    h1 = minimum_euclidean_distance(frequency, noisy_reflectivities, ks, thickness_step, E_oil, E_air, temp, salinity, theta);
    figure;
    histogram(h1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end