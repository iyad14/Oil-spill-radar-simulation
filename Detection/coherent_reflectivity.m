function R_coh = coherent_reflectivity(reflectivity, ks, theta)
    %ks: the surface roughness 
    %theta: incident angle of EM wave to interface
    x = ks*cosd(theta);
    R_coh = reflectivity.*exp(-4*(x^2));
end