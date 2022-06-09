%%  M                       --> number of scans
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




function [error, probability_of_error] = scatterplot(M, frequency, ks, thickness_step, t, variance, trials, E_oil, E_air, temp, salinity, theta)
        
    clf;                                    % clear figures
    thickness = 0:thickness_step:10;        % thickness over which the reflectivities will be calculated
    
    
        %% Creating the theoretical curve by which the measured reflectivities are compared
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    theoretical_reflectivity = abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        %% Generating the random reflectivities with given frequencies and at a given thickness to test the frequencies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    noise = sqrt(variance)*randn(length(frequency), trials, M);
    [~, indexOfThickness] = min(abs(thickness-t));
    noisy_reflectivities = 10*log10(abs(theoretical_reflectivity(:, indexOfThickness) + noise));
    noisy_reflectivities = sum(noisy_reflectivities, 3)/M;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
     %% Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if length(frequency) == 2
        
        % 2D Scatter plots 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
        scatter(noisy_reflectivities(1, :), noisy_reflectivities(2, :), 'x');
        hold on;
        scatter(10*log10(abs(theoretical_reflectivity(1, :))), 10*log10(abs(theoretical_reflectivity(2, :))), '.');
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
        grid on;
        hold off;   
        xlim([-8 0]);
        xlabel(strcat("Reflectivity(dB)@", num2str(frequency(1)), "GHz"));
        ylim([-8 0]);
        ylabel(strcat("Reflectivity(dB)@", num2str(frequency(2)), "GHz"));
        zlim([-8 0]);
        zlabel(strcat("Reflectivity(dB)@", num2str(frequency(3)), "GHz"));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
    end    
    
    
        %% Histograms 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%     h = double.empty;
%     for i = 1:1:length(noisy_reflectivities)
%         h(i) = minimum_euclidean_distance(frequency, transpose(noisy_reflectivities(:, i)), ks, thickness_step, E_oil, E_air, temp, salinity, theta);
%     end
%     histogram(h);
%     xlim([0 10]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
          %% Error & Percentage of correctness caclulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%Calculate the error!!!!!!!!!!!!!!!!!!!!!!!!
                         
    
    % Percentage of correctness
    %Number_of_correctly_estimated_thickness = y(indexOfThickness);                                                                  % find the number of correctly estimated thicknesses        
    %probability_of_error = (length(noisy_reflectivities) - Number_of_correctly_estimated_thickness)/length(noisy_reflectivities)    %divide wrong estimations by total estimations
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end