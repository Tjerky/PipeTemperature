%% Main file modelling the temperature distribution of a pipe
% define the constants
    % dimensions
R = 0.006; % m, outside radius of the tube
d = 0.001; % m, thickness of the tube wall
N = 100; % number of gridpoints in the r direction
P = 100; % number of gridpoints in the theta direction
    
    % material properties
rho = 295.5; % kg/m^3, density of the pipe material
k = 0.13; % W/(m K), thermal conductivity of the pipe material
cp = 1470; % J/(kg K), specific heat of the pipe material
c = 5000; % W/(m^2 K), heat transfer coefficient pipe -> water

% define the matrix
alpha = k/(rho*cp); % constant to simplify the formulas
beta = c/(rho*cp); % constant to simplify the formulas

% initial conditions

% loop
    % solve matrix equation for the temperature inside the pipe
    
    % derive the new temperture of the water
    
    % store the new timestep in array
    
% store the simulation in as a file