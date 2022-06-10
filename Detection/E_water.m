function E = E_water(temp, salinity, f)
    %frequency is given in GHz
    t = temp;
    S = salinity;
    
    %Conductivity
    
    A = [2.903602, 8.607e-2, 4.738817e-4, -2.991e-6, 4.3041e-9];
    sig35 = A(1) + A(2)*t + A(3)*t*t + A(4)*t*t*t + A(5)*t*t*t*t;

    A = [37.5109, 5.45216, 0.014409, 1004.75, 182.283];
    P = S * ((A(1) + A(2)*S + A(3)*S*S) / (A(4) + A(5)*S + S*S));

    A = [6.9431, 3.2841, -0.099486, 84.85, 69.024];
    alpha0 = (A(1) + A(2)*S + A(3)*S*S) / (A(4) + A(5)*S + S*S);

    A = [49.843, -0.2276, 0.00198];
    alpha1 = A(1) + A(2)*S + A(3)*S*S;
    
    Q = 1 + ((alpha0*(t-15))/(t+alpha1));

    sigma = sig35*P*Q;

    %some used parameters
    a=[0.46606917e-2 -0.26087876e-4 -0.63926782e-5 0.63000075e1 0.26242021e-2 -0.42984155e-2 ...
   0.34414691e-4 0.17667420e-3 -0.20491560e-6 0.58366888e3 0.12634992e3 0.69227972e-4 ...
   0.38957681e-6 0.30742330e3 0.12634992e3 0.37245044e1 0.92609781e-2 -0.26093754e-1];


    epsS = 87.85306*exp(-0.00456992*t - a(1)*S - a(2)*S*S - a(3)*S*t);
    epsOne = a(4)*exp(-a(5)*t-a(6)*S-a(7)*S*t);
    tau1 = (a(8)+a(9)*S)*exp(a(10)/(t+a(11)));
    tau2 = (a(12)+a(13)*S)*exp(a(14)/(t+a(15)));
    epsInf = a(16) + a(17)*t + a(18)*S;

    %Complex Permitivity Calculation
    eps = ((epsS-epsOne)./(1-1i*2*pi.*f.*tau1)) + ((epsOne-epsInf)./(1-1i*2*pi.*f.*tau2)) + (epsInf) + 1i*((17.9751*sigma)./f);
    
    %Seperates real and imaginary components
    epsr = real(eps);
    epsi = imag(eps);
 
    E = epsr - epsi*1i;              % According to Microwave Radar and Radiometric Remote Sensing by David Gardner Long , Fawwaz T. Ulaby p.124