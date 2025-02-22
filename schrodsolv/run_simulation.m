function run_simulation()
    % Define constants and initial conditions
    a = 260;      %% Length
    M = 1057.27504; %% Mass
    N = 1535;
    x = linspace(-a/2,a/2,N); x = x';
    k = (N*linspace(0,1/2,N)); k = k';
    dt = 10;     %% Time step
    dx = x(2) - x(1);
    V = 0;
    x0 = 75;
    delta = 1;
    E0 = 10;
    
    % Initialize Phi0
    %Phi0=init(N, x, M, x0, E0, delta)';
    Phi0 = sqrt(dx / (delta * sqrt(pi))) * exp(-0.5 * ((x - x0) / delta).^2) .* exp(-1i * sqrt(2 * M * E0) * x);
    %Phi0c = conj(Phi0); %% Real(Phi0) - i * imag(Phi0)

    % Define propagators
    GK = exp((1i * dt / (4 * M)) * ((2 * pi / a)^2) * (k.^2)); %% dt/2 kinetic energy propagator
    GK2 = exp((1i * dt / (2 * M)) * ((2 * pi / a)^2) * (k.^2)); %% dt kinetic energy propagator
    GV = exp(1i * dt * V); %% Potential spatial interaction

    % Initialize Phi
    iPhi = fft(Phi0);
    Phi = ifft(iPhi .* GK);
    Phi = GV .* Phi;

    % Parameters for the potential function
    NPt = 20000;
    ch = -0.15;
    lh = 100;
    hw = 1 ./ (1 + exp(ch * ((-lh:lh))))';
    reflectionRegion = x < x(lh);  % Define regions
    transmissionRegion = x > x(lh);
    % Initialize variables for results
    R = 0;
    T = 0;
    % Open file for writing
    fileID = fopen(sprintf('RT%0.3f-%d.txt', ch, lh), 'w');
    % Time-stepping loop
    for nrn = 1:NPt
        iPhi = fft(Phi);
        Phi = ifft(iPhi .* GK2);
        Phi = GV .* Phi;
        if mod(nrn, 1000) == 0
        amp = sum(abs(Phi).^2);
        fprintf(fileID, '%f\t%f\n', T, R);
        plot(x,abs(Phi))        
        pause(0.05)
        hold on
        end      
        Phi(1:2*lh) = Phi(1:2*lh) .* hw(1:2*lh);
R = trapz(x(reflectionRegion), abs(Phi(reflectionRegion)).^2)/ trapz(x, abs(Phi0).^2);
T = trapz(x(transmissionRegion), abs(Phi(transmissionRegion)).^2)/ trapz(x, abs(Phi0).^2);
    end
        hold off

    % Close file
    fclose(fileID);
end


function W2=init(NQ1, RQ1, RMS, R0, ESK, DELTA)
    % Initialize variables and constants
    %WI = [0.0, 1.0] / pi;
    WI =1i/pi;
    
    % Convert constants
    PI = pi;
    
    % Calculate AK
    AK = sqrt(2.0 * RMS * abs(ESK));
    if ESK <= 0.0
        AK = -AK;
    end
       
    % Initialize AE0
    
    % Calculate constants for the loop
    DR = RQ1(2) - RQ1(1);
    FAC = 1.0 / (PI * DELTA^2)^(0.25) * sqrt(DR);
    %FAC0 = sqrt(2.0 * RMS / PI);
    
    % Loop over NQ1
    for I = 1:NQ1
        % Calculate FACE and W2
        FACE = exp(-((RQ1(I) - R0) / DELTA)^2 / 2.0);
        W2(I) = FAC * exp(-1i * AK * RQ1(I)) * FACE;    
    end
end
