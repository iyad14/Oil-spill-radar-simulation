%%  measured_reflecitivity  --> reflectivities collected (linear scale)
%%  M                       --> Number of scans
%%  frequency               --> frequencies used given (in GHz)
%%  ks                      --> Surface roughness
%%  variance                --> Noise variance
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)]
%%  tmin & tmax             --> minimum and maximum value for thikness range
%%

function Detection = Detection_probability(R_oil, M, frequency, ks, variance, E_air, temp, salinity, theta, tmin, tmax)
   

    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water = ones(size(R_oil)) * ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water = abs(coherent_reflectivity(R_water, ks, theta));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Generating noise 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    trials = 1;          % Quantity of noisy reflectivity values needed to be generated
    
    % First & Second dimension  : image resolution
    % Third dimension  : number of frequencies used
    % Fourth dimension : Frequencies at which the scans are done
    
    % generates random values based on the given variance
    noise = sqrt(variance) * randn(length(R_oil), size(R_oil, 1), length(frequency), trials, M); 
  
    % Add the reflectivity values to the noise matrix, thus getting the noisy reflectivity values to be studied 
    noisy_reflectivity = noise + R_oil;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %%  Probability of detection calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the probabilty density functions
    h1 = pdf("Normal", noisy_reflectivity, R_oil, sqrt(variance));                                    % Oil
    h2 = pdf("Normal", noisy_reflectivity, R_water, sqrt(variance));                                                  % Water
   
    % Multiply the pdf values at all frequencies(f) and scans(M)
    x = prod(h1, 5);
    y = prod(h2, 5);
    
    x = prod(x, 3);
    y = prod(y, 3);

    % Check how many oil pdf values are greater than water pdf values
    result = gt(x, y);
    
    % Divide by the total number of trials to find the probability
    Detection(:, :) = sum(result, 4)/trials;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end