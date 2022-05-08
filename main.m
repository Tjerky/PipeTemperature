%% Main file modelling the temperature distribution of a pipe
% define the constants
    % dimensions
R = 0.006; % m, outside radius of the tube
d = 0.001; % m, thickness of the tube wall
N = 100; % number of gridpoints in the r direction
P = 100; % number of gridpoints in the theta direction
L = 1; % m, length of the pipe in the z direction
m = 1; % kg, mass of the water in the system

dx = d/N;
dy = 2*pi*R/P;

    % simulation
runtime = 100; % s, the time the simulation will run
dt = 0.01; % s, the timestep
    
    % material properties
        %pipe
rho = 295.5; % kg/m^3, density of the pipe material
k = 0.13; % W/(m K), thermal conductivity of the pipe material
cp = 1470; % J/(kg K), specific heat of the pipe material
c = 500; % W/(m^2 K), heat transfer coefficient pipe -> water
        %water
cw = 4180; % J/(kg K)
    % sunlight
I = 1000; % W/m^2, intensity of the sunlight

% define the matrix
alpha = k*dt/(rho*cp*dx^2); % constant to simplify the formulas
beta = c*dt/(rho*cp*dx); % constant to simplify the formulas
gamma = I*dt/(rho*cp*dx);

M = GenerateM(N, alpha); % generate the M matrix
Mi = inv(M); % generate the inverse of M
A = GenerateA(N, alpha, beta); % generate the A matrix

% generate the 3D array (time x Y x X) storing all the the temperature distribution
% of the pipe at every timepoint
T = zeros(runtime/dt, P, N);

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
Tw(1) = Tinitial-10;

% loop over all timepoints
for t = 2:(runtime/dt)
    % Generate the Q matrix
    Q = GenerateQ(squeeze(T(t-1, :, :)), Tw(t-1), beta, gamma, P, N);
    
    % loop over all pipe segments
    for i = 1:P       
        % get the specific q vector for this segment. This represents the
        % heat generation at every point.
        q = Q(:, i);
        
        % solve the matrix equation Mx = Ay + q
        T(t, i, :) = Mi*(A*squeeze(T(t-1, i, :)) + q);
    end
        
    % derive the new temperture of the water
        % get the temperature difference between all the water touching segments of
        % the pipes and the water.
    dT = Tw(t-1) - squeeze(T(t-1, :, 1));
    
        % derive the heat flux between the all the segments and the water
        % with newton's law of cooling
    qw = sum(-1 * c * dT * dy * L);
    
        % change the water temperate with euler's method and store it in
        % the Tw array
    Tw(t) = Tw(t-1) + qw/(cw * m);
end
    
% store the simulation in as a file