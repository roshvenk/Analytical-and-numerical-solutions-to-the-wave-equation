function [meanError,dt] = tutorial7(nx,tfinal)
% Description: This file implements the analytical vs numerical solution to the wave equation.
%
% Input parameters:
% - nx       - number of gridpoints in x
% - tfinal   - final simulation time (s)
%
% Output parameters:
% - meanError - mean absolute percentage error of numerical solution
% - dt        - maximum temporal step size for stability as per CFL condition
% Spatial (x) parameters
L      = 0.035;                     % size of the domain
xstart = 0;                         % start of computational domain (m)
xend   = xstart + L;                % end of computational domain (m)
dx     = (xend - xstart)/(nx - 1);  % spatial step size, delta x (m)
x      = xstart: dx : xend;         % vector of grid points from p to q

% Physical parameters.
l = 0.001;                          % height of fluid chamber (m)
K = 1;                             % membrane stiffness (N/m)
rho = 10.^3;                        % fluid density (kg/m^3)
c = sqrt(l*K/(2*rho));              % constant c in the governing equation

% Temporal (t) parameters
timestart = 0;                      % starting time (0 seconds);
dt = dx/c;
nt = (tfinal/dt)+1 ; 
nt = ceil(nt);
dt = tfinal/(nt-1);
time      = linspace(timestart, tfinal, nt);   % time vector


% Fourier parameters (coeffs and lambdas)
nfs = 3;                            % number of terms in the Fourier Series Expansion for analytical (since we only have 3 modes initially, the full solution will only have 3 modes)
A = 8.e-9;                          % coefficient of the first mode in the initial condition (m)
B = 5.e-9;                          % coefficient of the second mode in the initial condition (m)
C = 1.e-9;                          % coefficient of the third mode in the initial condition (m)
Bn = [A, B, C];                     % vector of fourier coefficients
lambda2 = c*2*pi/L;                     % lambda of the first mode
lambda5 = c*5*pi/L;                     % lambda of the second mode
lambda19 = c*19*pi/L;                   % lambda of the third mode
lambdan = [lambda2, lambda5, lambda19]; % vector of lambdas

% Initialise the solution arrays for the numerical solution: hn (h^n), hnp1 (h^(n+1)), hnm1 (h^(n-1)) as vectors of zeros with length nx
hn     = zeros(1, nx);    % current value of h (deflection) at time n
hnp1   = zeros(1, nx);    % next value of h at time n+1
hnm1   = zeros(1, nx);    % previous value at time n-1

% Initialise the solution arrays for the analytical solution: hExact as a vector of zeros with length nx 
hExact = zeros(1, nx);    % exact fourier series solution for h

% Specify the initial conditions for the numerical solution
for j=2:nx-1
    hn(j) = Bn(1)*sin(2*pi*x(j)/L) + Bn(2)*sin(5*pi*x(j)/L) + Bn(3)*sin(19*pi*x(j)/L);    % initial conditions
end


% Enforce the boundary conditions - these will never change in the
% simulation
hn(1)      = 0;
hn(nx)     = 0;
hnp1(1)    = 0;
hnp1(nx)   = 0;

% Set h^(n-1) to h^n to get things started
hnm1 = hn;


%% Main solution loop

% Loop through time
for i=1:nt

    t=time(i);                   % current time for outputted solution

    % Analytical solution (same as Week 5)
    hExact=zeros(1,nx);          % re-initialise the solution array
    for j=1:nx;                  % loop through grid points
       hExact(j)= hExact(j) + Bn(1)*cos(lambdan(1)*t)*sin(lambdan(1)*x(j)/c) + Bn(2)*cos(lambdan(2)*t)*sin(lambdan(2)*x(j)/c) + Bn(3)*cos(lambdan(3)*t)*sin(lambdan(3)*x(j)/c);
    end

    % Numerical solution - note that we only need to solve for i=2, nx-1 as
    % the first and last points are fixed by the boundary condition
    sigma=c^2*(dt^2)/(dx^2);
    for j=2:nx-1
        hnp1(j) = 2*hn(j) - hnm1(j)+sigma*(hn(j + 1) - 2*hn(j) + hn(j - 1));
    end

end
hnm1(2:nx-1) = hn(2:nx-1);
    hn(2:nx-1) = hnp1(2:nx-1);
meanError= mean(abs((hn(2: nx - 1) - hExact(2:nx-1))./hExact(2:nx-1)*100));

end