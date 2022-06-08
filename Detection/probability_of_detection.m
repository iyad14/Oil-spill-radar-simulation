function probability_of_detection = probability_of_detection(M, frequency, ks, variance, E_oil, E_air, temp, salinity, theta)

    step = 0.01;
    trials = 2500;


    thickness = 1:step:10;

    

    desired_reflectivity = reflectivity(frequency, thickness, ks, E_oil, E_air, temp, salinity, theta) ;
    E_water_prob = E_water(20, 35, frequency);
    R_water_fprob_planar = ones(size(desired_reflectivity)) * ((sqrt(E_air) - sqrt(E_water_prob))/(sqrt(E_air) + sqrt(E_water_prob)))^2;
    R_water_fprob_planar = coherent_reflectivity(R_water_fprob_planar, ks, theta);

    noise = sqrt(variance) * randn(length(desired_reflectivity), trials, M, size(desired_reflectivity, 1));
    A = transpose(abs(desired_reflectivity));
    A = reshape(A , [length(desired_reflectivity), 1, 1, size(noise, 4)]);
    noise = noise + A;

    h1 = pdf("Normal", noise, A, variance);
    h2 = pdf("Normal", noise, abs(R_water_fprob_planar(1,1)), variance);

    x = prod(h1, 4);
    y = prod(h2, 4);

    x = prod(x, 3);
    y = prod(y, 3);

    result = gt(x, y);

    probability_of_detection(:, :) = sum(transpose(result(:, :)) == 1)/trials;
    plot(thickness, probability_of_detection);
end