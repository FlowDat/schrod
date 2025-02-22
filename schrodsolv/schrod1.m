clear all; close all;
a = 260;  %% Length
M = 1057.27504; %% Mass
N = 1535;
x = linspace(-a/2,a/2,N); x = x';
k = (N*linspace(0,1/2,N)); k = k';
dt = 10; %% Time step
dx=x(2)-x(1);
%V0=0.1;
x0=75;
E0=10;
delta=1;
%V = -V0./cosh((x-30)).^2;
V = 0;
%V=0;
Phi0=sqrt(dx/(delta*sqrt(pi)))*exp(-0.5*((x-x0)/(delta)).^2).*exp(-1i*sqrt(2*M*E0).*x);
%Phi0=init(N, x, M, x0, E0, delta)';
%Phi0 = Phi0/sum(abs(Phi0).^2);
%Phi0c = conj(Phi0); %% real(Phi0)- i*imag(Phi0);
GK = (exp((i*dt/(4*M))*((2*pi/a)^2)*(k.^2))); %% dt/2 kinetic energy propagator
GK2 = (exp((i*dt/(2*M))*((2*pi/a)^2)*(k.^2))); %% dt kinetic energy propagator
GV = exp(i*dt*V); %% Potential spatial interaction
iPhi = fft(Phi0);
Phi = ifft(iPhi.*GK);
Phi = GV.*Phi;
NPt = 40000;
ch = -0.1;
lh = 100;
hw = 1./(1+exp(ch*((-lh:lh))))';
R=0;T=0;
kt=0;
fileID = fopen('RT.txt', 'w');
for nrn = 1:NPt
    iPhi = fft(Phi);    
    Phi = ifft(iPhi.*GK2);
    Phi = GV.*Phi;  
    if mod(nrn,2000) == 0    
    amp=sum(abs(Phi).^2);
    T = trapz(x(x<0), abs(Phi(x<0)).^2);
    R= trapz(x(x>0), abs(Phi(x>0)).^2);
    fprintf(fileID, '%f\t%f\n', T, R);
        fprintf('Norm %f \n',amp);
    %fprintf('Norm %f T %f R %f R+T %f R-T %f\n',amp,T,R, R+T, R-T);
        plot(x,abs(Phi))        
        pause(0.05)
        hold on
     
    end
   Phi(1:2*lh) = Phi(1:2*lh).*hw(1:2*lh);
end
fclose(fileID);
hold off

function W2=init(NQ1, RQ1, RMS, R0, ESK, DELTA)
    % Initialize variables and constants
    WI = [0.0, 1.0] / pi;
    
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
    FAC0 = sqrt(2.0 * RMS / PI);
    
    % Loop over NQ1
    for I = 1:NQ1
        % Calculate FACE and W2
        FACE = exp(-((RQ1(I) - R0) / DELTA)^2 / 2.0);
        W2(I) = FAC * exp(-1i * AK * RQ1(I)) * FACE;    
    end
end
