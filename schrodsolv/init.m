function init(NQ1, RQ1, RMS, R0, ESK, DELTA, W2, NE0, ESK0, AE0, W1)
    % Initialize variables and constants
    WI = [0.0, 1.0] / pi;
    
    % Convert constants
    PI = pi;
    
    % Calculate AK
    AK = sqrt(2.0 * RMS * abs(ESK));
    if ESK <= 0.0
        AK = -AK;
    end
    fprintf('AK = %f\n', AK);
    
    % Initialize AE0
    AE0(:) = 0.0;
    
    % Calculate constants for the loop
    DR = RQ1(2) - RQ1(1);
    FAC = 1.0 / (PI * DELTA^2)^(0.25) * sqrt(DR);
    FAC0 = sqrt(2.0 * RMS / PI);
    
    % Loop over NQ1
    for I = 1:NQ1
        % Calculate FACE and W2
        FACE = exp(-((RQ1(I) - R0) / DELTA)^2 / 2.0);
        W2(I) = FAC * exp(-1i * AK * RQ1(I)) * FACE;
        
        % Loop over NE0
        for NE = 1:NE0
            AKI = sqrt(2.0 * RMS * ESK0(NE));
            AE0(NE) = AE0(NE) + FAC0 / sqrt(AKI) * W2(I) * exp(1i * AKI * RQ1(I)) / 2.0 * sqrt(DR);
        end
        
        % Optionally print small W2 values
        if abs(W2(I)) > 1e-20
            fprintf('%f %e\n', RQ1(I), abs(W2(I)));
        end
    end
    
    % Output AE0 values
    for I = 1:NE0
        fprintf(24, '%f %f\n', ESK0(I) * 27.21139, abs(AE0(I))^2);
    end
end
