function tdse_solver
    % Constants
    hbar = 1; % Planck's constant over 2*pi
    m = 1; % Mass of the particle
    L = 10; % Length of the spatial domain
    N = 100; % Number of spatial points
    dt = 0.01; % Time step
    T = 1; % Total simulation time
    k = 2 * pi / L; % Wave number

    % Spatial discretization
    x = linspace(0, L, N);
    dx = x(2)-x(1); % Grid spacing

    % Initial condition (e.g., Gaussian wave packet centered at x=L/2)
    Psi = exp(-(x - L/2).^2) .* cos(k * (x - L/2));

    % Potential function (e.g., zero potential)
    V = @(x) zeros(size(x));

    % Pre-allocate matrices for efficiency
    H = zeros(N, N);
    V_matrix = zeros(N, 1);
    KineticOp = 1i * hbar^2 / (2 * m * dx^2) * [1 -2 1]; % 3-point stencil difference operator

    % Time evolution
    t = 0;
    while t < T
        % Calculate potential energy vector
        V_matrix = V(x);

        % Transparent boundary conditions (simple exponential decay)
        V_matrix(1) = V_matrix(1) * exp(-100 * abs((x(1) - L/2) / L));
        V_matrix(end) = V_matrix(end) * exp(-100 * abs((x(end) - L/2) / L));

        % Construct the kinetic energy operator
        H = KineticOp;

        % Crank-Nicolson method
        for iter = 1:100 % Number of Crank-Nicolson iterations (can be adjusted)
            Psi_half = Psi + 0.5 * dt * (H * Psi + V_matrix);
            H = circshift(KineticOp, [0, 1]); % Shift the kinetic operator for the next step
            Psi = Psi_half - dt * (H * Psi_half + V_matrix);
        end

        % Update time
        t = t + dt;

        % Optional: Visualization
        if mod(t/dt, 10) == 0 % Plot every 10 time steps
            plot(x, real(Psi));
            title(sprintf('Time t = %0.2f', t));
            drawnow;
        end
    end
end