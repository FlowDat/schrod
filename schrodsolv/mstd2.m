clear all
clc

% Constants and parameters
L = 10;            % Length of the spatial domain
N = 100;           % Number of spatial points
dx = L/N;           % Spatial step size
x = linspace(-L/2, L/2, N);  % Spatial grid

Tmax = 1;           % Maximum time
dt = 0.1;         % Time step size
timeSteps = Tmax/dt; % Number of time steps

% Initial wave packet parameters
x0 = -20;           % Initial position
k0 = 5;             % Wave number
sigma = 1;          % Width of the packet

% Potential (assume a simple potential barrier for example)
V = zeros(1, N);
V(abs(x) < 5) = 1;  % Potential barrier

% Initial wave packet
psi0 = exp(-(x - x0).^2 / (2 * sigma^2)) .* exp(1i * k0 * x);
psi = psi0;         % Initial wave function

% Absorbing boundary condition parameters
alpha = 0.1;       % Strength of the absorbing boundary

% Time evolution using Crank-Nicolson scheme
A = diag(1 + 2i*dt/(2*dx^2) - i*dt*V) + ...
    diag(-i*dt/(2*dx^2) * ones(1, N-1), 1) + ...
    diag(-i*dt/(2*dx^2) * ones(1, N-1), -1);

B = diag(1 - 2i*dt/(2*dx^2) + i*dt*V) + ...
    diag(i*dt/(2*dx^2) * ones(1, N-1), 1) + ...
    diag(i*dt/(2*dx^2) * ones(1, N-1), -1);

% Absorbing boundary condition (complex potential)
absorber = exp(-alpha * (abs(x) - L/2).^2);
absorber(abs(x) < L/4) = 1;  % Only apply near boundaries

% Time evolution loop
for t = 1:timeSteps
    psi = absorber .* (A \ (B * psi.'));
    %if mod(t, 10) == 0
    plot(abs(psi))% Crank-Nicolson update
    %end
end
plot(abs(psi))
% Calculate reflection and transmission coefficients
reflectionRegion = x < -L/4;  % Define regions
transmissionRegion = x > L/4;

R = trapz(x(reflectionRegion), abs(psi(reflectionRegion)).^2);
T = trapz(x(transmissionRegion), abs(psi(transmissionRegion)).^2);

reflectionCoefficient = R / trapz(x, abs(psi0).^2);
transmissionCoefficient = T / trapz(x, abs(psi0).^2);

% Display results
fprintf('Reflection Coefficient: %.4f\n', reflectionCoefficient);
fprintf('Transmission Coefficient: %.4f\n', transmissionCoefficient);
