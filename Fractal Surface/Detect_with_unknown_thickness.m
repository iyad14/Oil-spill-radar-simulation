%%  measured_reflecitivity  --> reflectivities collected (linear scale)
%%  M                       --> Number of scans
%%  frequency               --> frequencies used given (in GHz)
%%  ks                      --> Surface roughness
%%  variance                --> Noise variance
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)]
%%

function oil_found = Detect_with_unknown_thickness(measured_reflectivity, M, frequency, ks, variance, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)
   
    thickness = tmin:thickness_step:tmax;
    
    R_oil = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta);
    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water = ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water = abs(coherent_reflectivity(R_water, ks, theta));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
        %%  Generating noise 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % First & Second dimension  : image resolution
    % Third dimension  : number of frequencies used
    % Fourth dimension : Frequencies at which the scans are done
    
    % generates random values based on the given variance
    noise = sqrt(variance) * randn(size(measured_reflectivity, 1), size(measured_reflectivity, 2), M); 
  
    % Add the reflectivity values to the noise matrix, thus getting the noisy reflectivity values to be studied 
    noisy_reflectivity = noise + measured_reflectivity;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %%  Probability of detection calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculate the probabilty density functions
    h1 = pdf("Normal", noisy_reflectivity, R_oil, sqrt(variance));                                    % Oil
    h2 = pdf("Normal", noisy_reflectivity, R_water, sqrt(variance));                                                  % Water
   
    x = sum(h1, 2)/length(h1);
 
    % Multiply the pdf values at all frequencies(f) and scans(M)
    x = prod(x, 3);
    y = prod(h2, 3)*0.5;
    
    x = prod(x, 1);
    y = prod(y, 1);
    
    if x/y >= 1
        oil_found = 1;
    else
        oil_found = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     
 
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end