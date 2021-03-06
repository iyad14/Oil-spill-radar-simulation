%%  frequency               --> frequencies used given in GHz
%%  thickness               --> thickness range at which the system will operate (given in mm)
%%  ks                      -->  Surface roughness
%%  E_oil                   --> Dielectric constant of oil
%%  E_air                   --> Dielectric constant of air
%%  temp                    --> Temperature of water (Degrees Celsius)
%%  salinity                --> Salinity of water (ppt)
%%  theta                   --> Incident angle of the electromagnetic wave to interface (given in degrees)
%%



function R = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta) 
    

        %%  Water Dielectric constant and field reflection coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    E_Water = E_water(temp, salinity, frequency);
    raw12 = field_reflection_coefficient(sqrt(E_air), sqrt(E_oil));                 
    raw23  = field_reflection_coefficient(sqrt(E_oil), sqrt(E_Water));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
        %%  Phase shift & Reflectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta = (2*pi*10^9*sqrt(E_oil)*thickness*10^(-3))/(3*10^8);
    delta = transpose(frequency).*delta;
    
    num = raw12*exp(1i*delta) + transpose(raw23).*exp(-1i*delta);
    denum = exp(1i*delta) + (raw12 .* transpose(raw23) .* exp(-1i*delta));
    R = abs(num./denum).^2;
    
     %% Other form of reflectivity equation 
    %     num = (raw12.^2 + transpose(raw23.^2)  + 2*raw12.*transpose(raw23).*cos(2.*delta));
    %     denum = 1 + transpose((raw12.*raw23).^2)+2*raw12.*transpose(raw23).*cos(2.*delta);
    %     R = num./denum;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %% Coherent Reflectivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = R.*exp(-4*((ks*cosd(theta))^2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end

