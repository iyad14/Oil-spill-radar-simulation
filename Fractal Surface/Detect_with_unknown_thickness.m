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

function h1 = Detect_with_unknown_thickness(measured_reflectivities, M, frequency, ks, variance, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)
   
    thickness = tmin:thickness_step:tmax;      %thickness over which the reflectivities will be calculated
    R_oil = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta) .* ones(size(measured_reflectivities));
    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water = ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water = abs(coherent_reflectivity(R_water, ks, theta));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    

    % columns --> thickness, rows -- >  frequency, 3d --> M
    for i = 1:1:length(thickness)
        h1(i) = pdf("Normal", measured_reflectivities(i), R_oil(i), sqrt(variance));
        h2(i) = pdf("Normal", measured_reflectivities(i), R_water, sqrt(variance));
    end

    probability_of_detection = prod(prod(sum(h1, 2), 3), 1)/ prod(prod(sum(h2, 2), 3), 1);
     
 
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end