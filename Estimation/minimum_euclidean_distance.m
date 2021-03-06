%%  frequency               --> frequencies used given (in GHz)
%%  measured_reflectivity   --> Reflectivitis measured by drone (in dB scale)
%%  ks                      -->  Surface roughness
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (in ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)
%%  tmin & tmax             --> minimum and maximum value for thikness range
%%  thickness_step          --> Thickness resolution (in mm)
%%



function t = minimum_euclidean_distance(measured_reflectivity, frequency, ks, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)
    
    thickness = tmin:thickness_step:tmax;    %thickness over which the reflectivities will be calculated

        %% Creating the theoretical curve by which the measured reflectivities are compared
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    theoretical_reflectivity = 10*log10(abs(reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Finding the point closest to the measured reflectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    distance = sqrt(sum((theoretical_reflectivity-transpose(measured_reflectivity)).^2, 1));
    [~ , indexOfThickness] = min(distance);
    t = thickness(indexOfThickness);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end