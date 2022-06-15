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

function probability_of_detection = Detect_with_unknown_thickness(measured_reflectivities, M, frequency, ks, variance, E_oil, E_air, temp, salinity, theta, tmin, thickness_step, tmax)
   
    pW = 0.5;
    thickness = tmin:thickness_step:tmax;      %thickness over which the reflectivities will be calculated
    R_oil = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta);
    
        %%  Dielectric constant of water at given frequencies & water reflectivity values at a given surface roughness ks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    E_water_prob = E_water(temp, salinity, frequency);
    R_water = ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water = abs(coherent_reflectivity(R_water, ks, theta));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trials = 5000;
    pO = probability(M, frequency, ks, thickness_step, variance, trials, E_oil, E_air, temp, salinity, theta, tmin, tmax) ;
    
    for m = 1:1:M      
            h1(:, :, m) = pdf("Normal", transpose(measured_reflectivities(m, :)), R_oil, sqrt(variance));
            h2(:, :, m) = pdf("Normal", transpose(measured_reflectivities(m, :)), R_water, sqrt(variance));
    end
  
    plot(h1(1, :, 1))
    hold on
    plot(h2(1, :, 1))
    
    h1 = h1.*pO;
    h2 = h2*pW;
    
    x = sum(h1, 2)/10;
    
    x = prod(x, 3);
    y = prod(h2, 3);
    
    x = prod(x, 1);
    y = prod(y, 1);

    probability_of_detection = x/y;
     
 
end



function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end