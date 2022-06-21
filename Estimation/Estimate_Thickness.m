%%  measured_reflecitivity  --> reflectivities collected (linear scale)
%%  M                       --> number of scans
%%  frequency               --> frequencies used given (in GHz)
%%  ks                      -->  Surface roughness
%%  thickness_step          --> Thickness resolution (in mm)
%%  variance                --> Noise variance
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)
%%  tmin & tmax             --> minimum and maximum value for thikness range
%%


function thicknesses = Estimate_Thickness(measured_reflectivity, M, frequency, ks, thickness_step, variance,E_oil, E_air, temp, salinity, theta, tmin, tmax)

    
    
        %% Simulate the M observations
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
     % generate noise M times to simulate the collection of different observations
     noise = sqrt(variance) * randn(length(measured_reflectivity), size(measured_reflectivity, 1), length(frequency), M); 
     measured_reflectivity = abs(measured_reflectivity + noise);
     measured_reflectivity = sum(measured_reflectivity, 4)/M;
 
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        %% Estimation & Histogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    thicknesses = double.empty;
    
    %   use minimum euclidean distance to find estimate the thickness on
    %   each reflectvity value collected
    for j = 1:1:size(measured_reflectivity,2)
        for i = 1:1:size(measured_reflectivity, 1)
            thicknesses(i, j) = minimum_euclidean_distance(frequency, measured_reflectivity(i, j), ks, thickness_step, E_oil, E_air, temp, salinity, theta, tmin, tmax);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end