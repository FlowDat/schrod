% Parameters
L = 130;                   % Length of the domain
N = 100;                   % Number of spatial points
M = 200;                   % Number of time steps

% Check that N and M are positive integers
if ~isscalar(N) || ~isscalar(M) || N <= 0 || M <= 0
    error('N and M must be positive scalars.');
end

X = linspace(0, L, N);    % Spatial grid
h = X(2) - X(1);          % Grid spacing
T = 10;                   % Total time
tau = T / M;              % Time step size
alpha = 1;                % Diffusion coefficient

% Initialize matrices
A = zeros(N, N);          % Stiffness matrix
B = zeros(N, N);          % Mass matrix
M_matrix = zeros(N, N);  % Modified mass matrix

% Define the finite element basis functions and their derivatives
for i = 1:N-1
    % Local basis functions (linear)
    xi = X(i:i+1);
    h_local = xi(2) - xi(1);
    B(i, i) = B(i, i) + 1/h_local;
    B(i, i+1) = B(i, i+1) - 1/h_local;
    B(i+1, i) = B(i+1, i) - 1/h_local;
    B(i+1, i+1) = B(i+1, i+1) + 1/h_local;
    
    A(i, i) = A(i, i) + 1/(2*h_local);
    A(i, i+1) = A(i, i+1) - 1/(2*h_local);
    A(i+1, i) = A(i+1, i) - 1/(2*h_local);
    A(i+1, i+1) = A(i+1, i+1) + 1/(2*h_local);
end

% Initial condition
Psi = zeros(N, M+1);   % Initialize wave function matrix
Psi(:,1) = exp(-0.5 * ((X - L/2) / 5).^2);  % Gaussian initial condition

% Define the TBC kernel (example: using a simple exponential decay kernel)
n_kernel = 5;  % Length of the kernel
K = exp(-linspace(0, 2, n_kernel));  % Example kernel
K = K / sum(K);  % Normalize kernel

% Time-stepping loop
for m = 1:M
    % Construct the modified mass matrix
    M_matrix = eye(N) + (tau / 2) * alpha * A;
    
    % Construct the right-hand side
    b = (eye(N) - (tau / 2) * alpha * A) * Psi(:,m);
    
    % Apply the TBC at the boundary
    % Convolve the last few elements of Psi with the kernel K
    if N >= n_kernel
        boundary_vals = Psi(end-n_kernel+1:end, m);  % Get the last values
        convolved_vals = K * boundary_vals;  % Convolution
        
        % Update the boundary value
        b(end) = convolved_vals(1);  % Example: use first value of convolution for boundary condition
    end
    
    % Solve the linear system
    Psi(:, m+1) = M_matrix \ b;
end

% Plot the results
figure;
mesh(linspace(0, T, M+1), X, Psi');
xlabel('Time');
ylabel('Position');
zlabel('Wave Function');
title('1D FEM with Discrete TBCs for the Schrödinger Equation');
