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


function Estimated_thicknesses = Estimate_Thickness(measured_reflectivity, M, frequency, ks,  variance, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)



    %% Simulate the M observations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % generate noise M times to simulate the collection of different observations
    noise = sqrt(variance) * randn(size(measured_reflectivity, 1), size(measured_reflectivity, 2), length(frequency), M);
    measured_reflectivity = measured_reflectivity + noise;
    measured_reflectivity = sum(measured_reflectivity, 4)/M;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Estimation & Histogram
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Estimated_thicknesses = double.empty;
 
    %   use minimum euclidean distance to find estimate the thickness on
    %   each reflectvity value collected
    for i = 1:1:size(measured_reflectivity, 1)
       for j = 1:1:size(measured_reflectivity, 2)
        for k = 1:1:size(measured_reflectivity, 3)
            ref = 10*log10(measured_reflectivity(i, j, :));
            Estimated_thicknesses(i, j) = minimum_euclidean_distance(reshape(ref, [1, size(ref, 3)]), frequency, ks, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax);
        end
       end
    end
       
end