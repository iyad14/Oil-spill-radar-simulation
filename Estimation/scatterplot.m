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




function [Estimated_thickness, MaximumError, probability_of_error] = scatterplot(M, frequency, ks, t, variance, trials, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)
        
    clf;                                         % clear figures
    thickness = tmin:thickness_step:tmax;        % thickness over which the reflectivities will be calculated
    
    
        %% Creating the theoretical curve by which the measured reflectivities are compared
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    R_oil = abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta)); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
        %% Generating the random reflectivities with given frequencies and at a given thickness to test the frequencies
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    noise = sqrt(variance)*randn(length(frequency), trials, M);
    [~, indexOfThickness] = min(abs(thickness-t));
    noisy_reflectivities = 10*log10(abs(R_oil(:, indexOfThickness) + noise));
    noisy_reflectivities = sum(noisy_reflectivities, 3)/M;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
     
     %% Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if length(frequency) == 2
        
        % 2D Scatter plots 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
       
            %% Plotting regions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        SymbolSet = zeros(size(R_oil,2),1);
        SymbolSet(:,1)= 10*log10(R_oil(1,:).');
        SymbolSet(:,2)= 10*log10(R_oil(2,:).');
        [NumOfSignal, ~] = size(SymbolSet);
        ProbSet = ones(1,NumOfSignal)/NumOfSignal; % probability of each symbol
        No = 2;
        %Generate the Decision Region Diagram of the input Signal Symbol Set
        hold on;
        SignalSymbolDecisionRegionGenerator(SymbolSet, ProbSet, No);
        hold on;
        plot(SymbolSet(:,1),SymbolSet(:,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        % Noise Scattering
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        scatter(noisy_reflectivities(1, :), noisy_reflectivities(2, :), 'x');
        hold on;
        scatter(10*log10(abs(R_oil(1, :))), 10*log10(abs(R_oil(2, :))), '.');
        grid on;
        hold off;
        xlim([-8 0]);
        xlabel(strcat("Reflectivity(dB)@", num2str(frequency(1)), "GHz"));
        ylim([-8 0]);
        ylabel(strcat("Reflectivity(dB)@", num2str(frequency(2)), "GHz"));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           
    elseif length(frequency) == 3  
        
        % 3D Scatter plots 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        scatter3(noisy_reflectivities(1, :), noisy_reflectivities(2, :), noisy_reflectivities(3, :), 'x');
        hold on;
        scatter3(10*log10(abs(R_oil(1, :))), 10*log10(abs(R_oil(2, :))), 10*log10(abs(R_oil(3, :))), '.');   
        grid on;
        hold off;   
        xlim([-8 0]);
        xlabel(strcat("Reflectivity(dB)@", num2str(frequency(1)), "GHz"));
        ylim([-8 0]);
        ylabel(strcat("Reflectivity(dB)@", num2str(frequency(2)), "GHz"));
        zlim([-8 0]);
        zlabel(strcat("Reflectivity(dB)@", num2str(frequency(3)), "GHz"));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
    else
        
    end    
    
    
        %% Histograms 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    figure;
    h = double.empty;
    for i = 1:1:length(noisy_reflectivities)
        h(i) = minimum_euclidean_distance(transpose(noisy_reflectivities(:, i)), frequency, ks, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax);
    end
    s = histogram(h);
    xlim([tmin tmax+1]);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
          %% Maximum Error & Probability of error caclulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Error Calculation
    BinCenters(1) = s.BinLimits(1) + s.BinWidth/2;          % Find the thickness at the which the farthest bin to the left from the actual thickness is plotted
    BinCenters(2) = s.BinLimits(2) - s.BinWidth/2;          % Find the thickness at the which the farthest bin to the right from the actual thickness is plotted
    distance = abs(t - BinCenters);                         % Find the deveation from the actual value
    [MaximumError, ~] = max(distance);                      % The error represents the value with maximum deveation from the actual value
    
    
    % Probability of error
    hist_thickness = s.BinEdges - s.BinWidth/2;                                                                                         % Find the centers of the bins created which represent each thickness
    hist_thickness(1) = [];                                                                                                             % Delete the first element since it represents nothing ( first binEdge-binWidth/2)
    [~, indexOfActualThickness] = min(abs(t- hist_thickness));                                                                          % Find the index of the actual thickness
    Number_of_correctly_estimated_thickness = s.Values(indexOfActualThickness);                                                         % find the number of correctly estimated thicknesses        
    probability_of_error = (length(noisy_reflectivities) - Number_of_correctly_estimated_thickness)/length(noisy_reflectivities);       %divide wrong estimations by total estimations
    
    % Estimated thickness
    [~, indexOfEstimatedThickness] = max(s.Values); 
    Estimated_thickness = hist_thickness(indexOfEstimatedThickness);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end