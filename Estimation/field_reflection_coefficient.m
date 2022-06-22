%%  dielectric_constant1               --> dielectric constant of material on top layer
%%  dielectric_constant2               --> dielectric constant of material on bottom layer
%%  r                                  --> field reflection coefficient    

function r = field_reflection_coefficient(dielectric_constant1, dielectric_constant2)
    r = (dielectric_constant1 - dielectric_constant2)./(dielectric_constant1 + dielectric_constant2);
end