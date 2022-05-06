%% Main file modelling the temperature distribution of a pipe
% define the constants
    % dimensions
R = 0.006; % m, outside radius of the tube
d = 0.001; % m, thickness of the tube wall
N = 100; % number of gridpoints in the r direction
P = 100; % number of gridpoints in the theta direction

dx = d/N;
dy = 2*pi*R/P;

    % simulation
runtime = 10; % s, the time the simulation will run
dt = 0.01; % s, the timestep
    
    % material properties
rho = 295.5; % kg/m^3, density of the pipe material
k = 0.13; % W/(m K), thermal conductivity of the pipe material
cp = 1470; % J/(kg K), specific heat of the pipe material
c = 5000; % W/(m^2 K), heat transfer coefficient pipe -> water

% define the matrix
alpha = k/(rho*cp); % constant to simplify the formulas
beta = c/(rho*cp); % constant to simplify the formulas

M = GenerateM(N, alpha); % generate the M matrix
A = GenerateA(N, alpha, beta); % generate the A matrix

% generate the 3D array (time x X x Y) storing all the the temperature distribution
% of the pipe at every timepoint
T = zeros(runtime/dt, N, P);

% generate the 1D array(time) storing the water temperature at every
% time-point
Tw = zeros(runtime/dt, 1);

% initial conditions
    % set initial temperature
Tinitial = 295; % K, initial temperature of the water and the pipe

    % set initial temperature distribution to be the initial temperature
    % everywhere
T(1, :, :) = Tinitial * ones(N, P);
    % set initial water temperature to be the initial temperature
Tw(1) = Tinitial;

% loop
    % solve matrix equation for the temperature inside the pipe
    
    % derive the new temperture of the water
    
    % store the new timestep in array
    
% store the simulation in as a file