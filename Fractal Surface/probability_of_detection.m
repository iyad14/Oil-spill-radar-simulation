%%  M                       --> Number of scans
%%  frequency               --> frequencies used given (in GHz)
%%  ks                      --> Surface roughness
%%  thickness_step          --> Thickness resolution (in mm)
%%  variance                --> Noise variance
%%  trials                  --> Quantity of noisy reflectivity values needed to be generated
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)]
%%

function probability_of_detection = probability(M, frequency, ks, thickness_step, variance, trials, E_oil, E_air, temp, salinity, theta, tmin, tmax)

    thickness = tmin:thickness_step:tmax;      %thickness over which the reflectivities will be calculated
    
    
        %%  Calculate reflectivities for given frequencies, thicknesses, and surface roughness
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    theoretical_reflectivity = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water_fprob_planar = ones(size(theoretical_reflectivity)) * ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water_fprob_planar = coherent_reflectivity(R_water_fprob_planar, ks, theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Generating noise 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % First dimension  : thickness values
    % Second dimension : Noise values
    % Third dimension  : Number of scans (M)
    % Fourth dimension : Frequencies at which the scans are done
    
    % generates random values based on the given variance
    noise = sqrt(variance) * randn(length(theoretical_reflectivity), trials, M, size(theoretical_reflectivity, 1)); 
    
    % transpose the theoretical reflectivity to add the reflectivity values to the noise generated (each reflecctivity for its correspondant thickness)
    A = transpose(abs(theoretical_reflectivity));  
    
    % reshape the reflectivity values to suit the dimensions of the noise
    % matrix
    A = reshape(A , [length(theoretical_reflectivity), 1, 1, size(noise, 4)]);
    
    % Add the reflectivity values to the noise matrix, thus getting the noisy reflectivity values to be studied 
    noisy_reflectivity = noise + A;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %%  Probability of detection calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the probabilty density functions
    h1 = pdf("Normal", noisy_reflectivity, A, sqrt(variance));                                    % Oil
    h2 = pdf("Normal", noisy_reflectivity, abs(R_water_fprob_planar(1,1)), sqrt(variance));       % Water
    
    % Multiply the pdf values at all frequencies(f) and scans(M)
    x = prod(h1, 4);
    y = prod(h2, 4);

    x = prod(x, 3);
    y = prod(y, 3);

    % Check how many oil pdf values are greater than water pdf values
    result = gt(x, y);
    
    % Divide by the total number of trials to find the probability
    probability_of_detection(:, :) = sum(transpose(result(:, :)) == 1)/trials;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end