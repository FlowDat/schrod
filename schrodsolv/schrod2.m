clear all; close all;

m = 1057.27504;
L = 260;
N = 1535;
x0 = 75;
lm = 90;
dx = L/(N-1);
E0 = 10;
delta = 1;
dt = 10;
NPt = 40000;
k = ((1:N)*pi/L).^2/2/m;
x = linspace(-L/2,L/2,N);
%k=fftshift(k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mk = sqrt(2/L)*sin(((pi/L*(1:N)'*(x+L/2))));
%D2 = diff(Mk,2)/dx;
D2=gradient((gradient(Mk,dx)),dx);
%D2=D2/norm(D2,2);
%plot(D2)
[eigenvecs, eigenvals] = eig(D2/2/m);
energies = (sort(abs(diag(eigenvals))))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pk = exp(1i*dt*k);
%pk=fftshift(exp(1i*dt*k));
%pk1 = exp(1i*dt*energies);

V = 0;
psi0 = sqrt(dx/(delta*sqrt(pi)))*exp(-0.5*((x-x0)/(delta)).^2).*exp(-1i*sqrt(2*m*E0).*x);
psi = psi0;
%%%%hw=1./(1+exp(-0.16*(-40:59)));%1
hw = 1./(1+exp(-0.1*((-lm:lm))));%2

for t = 1:NPt
psi = exp(-1i*(V*dt)/2).*psi;
psi = ifft(fft(psi).*pk);
psi = exp(-1i*(V*dt)/2).*psi;
psi(1:2*lm) = psi(1:2*lm).*hw(1:2*lm);%2
%%%%psi(1:100)=psi(1:100).*hw(1:100);%1

    if  mod(t,1000)==0
    amp = max(abs(psi).^2);
    fprintf('Norm %f\n',amp);
    plot(x,(abs((psi))));%axis([-150 150 0 0.15]);
    pause(0.1)
    hold on
    end    
    
end

    hold off