function [y,freqLast] = tutorial6(nt,nfs,L)
% Description: This file implements the analytical solution to the wave equation (guitar string).
%
% Input parameters:
% - nt      - number of timesteps
% - nfs     - number of fourier terms
% - L       - Size of the Domain
%
% Output parameter:
% - y        - guitar string deflection (m)
% - freqLast - last fundamental frequency of system (Hz)
% Spatial (x) parameters
nx     = 101;                       % Number of grid points
xstart = 0;                         % Start of computational domain (m)
xend   = xstart + L;                % End of computational domain (m)
dx     = (xend - xstart)/(nx-1);                        % Spatial step size, delta x (m)
x      = linspace(xstart,xend,nx);                        % Vector of grid points
 
% Temporal (t) parameters
tstart = 0;                         % Starting time (s)
dt     = 0.00001;                   % Time step (s)
timeend   = dt*(nt-1);
time   = tstart:dt:timeend;                        % Time vector 

% Phyiscal parameters
T      = 71;                        % Tension in string (N)
rho    = 0.401e-3;                  % Mass per unit length (kg/m)
c      = sqrt(T/rho);               % PDand E constant 
k      = 0.005;                     % Initial displacement (m) 

%Initialise the arrays B and lambda for all Fourier coeffs
B=zeros(1,nfs);
lambda=zeros(1,nfs);
for n=1:nfs 
    B(n)=8*k/(n^2*pi^2)*sin((n*pi)/2);
    lambda(n)=(c*n*pi/L);
end
freqLast = lambda(nfs)/(2*pi);

%% Main solution loop

% Loop through time
for i=1:length(time)
    
    t = time(i);        % time at current point in loop

    y=zeros(1,nx);               % Initialise the solution vector for the analytical solution

    % Loop through space
    for j=1:nx         % For each grid point

        for n=1:nfs    % Sum each individual term in the series
            y(j)=y(j)+B(n)*cos(lambda(n)*t)*sin(n*pi*x(j)/L);
        end
    end
end