%%  M                       --> Number of scans
%%  frequency               --> frequencies used given (in GHz)
%%  ks                      --> Surface roughness
%%  thickness_step          --> Thickness resolution (in mm)
%%  variance                --> Noise variance
%%  samples                 --> Quantity of noisy reflectivity values needed to be generated
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)]
%%  tmin & tmax             --> minimum and maximum value for thikness range
%%

function probability = probability_of_detection(M, frequency, ks, thickness_step, variance, samples, E_oil, E_air, temp, salinity, theta, tmin, tmax)

    % Changes in cases of witnesses
    pO = 0.5;   % Assuming p(O) = 0.5
    pW = 0.5;   % Assuming p(W) = 0.5 
    
    
    thickness = tmin:thickness_step:tmax;      %thickness over which the reflectivities will be calculated
    
    
        %%  Calculate reflectivities for given frequencies, thicknesses, and surface roughness
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    R_oil = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water_fprob_planar = ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water_fprob_planar = coherent_reflectivity(R_water_fprob_planar, ks, theta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Generating noise 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % First dimension  : thickness values
    % Second dimension : Noise values
    % Third dimension  : Number of scans (M)
    % Fourth dimension : Frequencies at which the scans are done
    
    % generates random values based on the given variance
    noise = sqrt(variance) * randn(size(R_oil, 2), samples, M, size(R_oil, 1)); 
    
    % transpose the theoretical reflectivity to add the reflectivity values to the noise generated (each reflecctivity for its correspondant thickness)
    A = transpose(abs(R_oil));  
    
    % reshape the reflectivity values to suit the dimensions of the noise
    % matrix
    A = reshape(A , [size(R_oil, 2), 1, 1, size(noise, 4)]);
    
    % Add the reflectivity values to the noise matrix, thus getting the noisy reflectivity values to be studied 
    noisy_reflectivity = noise + A;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %%  Probability of detection calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the probabilty density functions
    h1 = pdf("Normal", noisy_reflectivity, A, sqrt(variance));                                    % Oil
    h2 = pdf("Normal", noisy_reflectivity, abs(R_water_fprob_planar), sqrt(variance));            % Water
    
    % Multiply the pdf values at all frequencies(f) and scans(M)
    x = prod(h1, 4) * pO;              % p(O) = 0.5
    y = prod(h2, 4) * pW;              % p(W) = 0.5

    x = prod(x, 3);
    y = prod(y, 3);

    % Check how many oil pdf values are greater than water pdf values
    result = gt(x, y);
    
    % Divide by the total number of samples to find the probability
    probability(:, :) = sum(transpose(result(:, :)) == 1)/samples;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end