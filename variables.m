%% define the constants
    % dimensions of the storage vessel
Rmax = 0.055; % m, outside radius of the tube
Rmin = 0.025; % m, inside radius of the tube
R1 = 0.0261; % m, outside radius of the first tube
R2 = 0.0539; % m, inside radius of the second tube
L_pvc = 0.55; % m, length of the pvc pipes

    % dimensions of the polyurethane tubing
L_pol = 1; % m, length of the polyurethane tube
R_pol_o = 0.006; % m, outside radius of the polyurethane tube
R_pol_i = 0.004; % m, inside radius of the polyurethane tube

    % dimensions of the solar collector
A = 1.1; % m^2, area of the aluminium plate
d = 0.002; % m, thickness of the aluminium plate
L_cop = 5; % m, length of the copper pipe
R_cop_o = 0.006; % m, outside radius of the copper pipe
R_cop_i = 0.00545; % m, inside radius of the copper pipe
H = 0.007; % m, distance between aluminium plate and the glass plate
    
    % material properties
        % pvc pipe
k_pvc = 0.19; % W/(m K), thermal conductivity of the pipe material (pvc)
hw_pvc = 3000; % W/(m^2 K), heat transfer coefficient pipe -> water
hv_pvc = 0.1; % W/(m^2 K), heat transfer coefficient pipe -> partial vacuum
e_pvc = 0.93; % unitless, emissivity of pvc
        % vacuum
k_v = 1; % W/(m K), thermal conductivity of the vaccuum
        %water
cw = 4180; % J/(kg K)
mw = 1; % kg, mass of the water in the system
hw_cop = 300; % W/(m^2 K), heat transfer coefficient copper -> water
hw_pol = 10; % W/(m^2 K), heat transfer coefficient polyurethane -> water
        % air
rho_air = 1.225; % kg/m^3, density of air
c_air = 718; % J/(kg K), specific heat of air
ha_pvc = 10; % W/(m^2 K), heat transfer coefficient pvc -> air
ha_cop = 10; % W/(m^2 K), heat transfer coefficient copper -> air
ha_al = 10; % W/(m^2 K) heat transfer coefficient aluminium -> air
ha_pol = 10; % W/(m^2 K) heat transfer coefficient polyurethane -> air
e_air = 0.8; % unitless, emissivity of the air
Tair = 295; % K, temperature of the air
        % aluminium
rho_al = 2702; % kg/m^3, density of aluminium
c_al = 880; % J/(kg K), specific heat of aluminium 
k_al = 237; % W/(m K) thermal conductivity of aluminium
e_al = 0.1; % unitless, emissivity of the aluminium plate
        % copper
k_cop = 402; % W/(m K), thermal conductivity of copper
e_cop = 0.03; % unitless, emissivity of the copper tube
        % polyurethane
k_pol = 0.13; % W/(m K), thermal conductivity of polyurethane
    
    % sun
I = 1000; % W/m^2, intensity of the sunlight
    % Univeral constants
sigma = 5.670374419 * 10^-8; % W(m^2 K^4), Stefan-Boltzmann constant

    % simulation
runtime = 20*60; % s, duration of the simulated period
dt = 1; % s, duration of one timestep

P = ceil(runtime/dt); % unitless, amount of timesteps

    % initial condition
T_initial = 295; % K, initial temperature of the water

%% Calculate thermal conductivities of each layer in W/K (storage tank)
    % 1. Interface between the water and the inner surface of the inner pvc
    % pipe
k1 = hw_pvc * 2*pi*Rmin*L_pvc;
    % 2. Conduction through the inner pvc pipe
k2 = 2*pi*k_pvc*L_pvc/log(R1/Rmin);
    % 3. Convection between the partial vacuum and the outer surface of the
    % inner pvc tube
k3 = hv_pvc * 2*pi*R1*L_pvc;
    % 4. Conduction through the partial vacuum
k4 = 2*pi*k_v*L_pvc/log(R2/R1);
    % 5. Convention between the partial vacuum and the inner surface of the
    % outer pvc tube
k5 = hv_pvc * 2*pi*R2*L_pvc;
    % 6. Conduction through the outer pvc pipe
k6 = 2*pi*k_pvc*L_pvc/log(Rmax/R2);
    % 7. Convection between the air and the outer surface of the outer pvc
    % pipe
k7 = ha_pvc * 2*pi*Rmax*L_pvc;

%% Reduce the conductivities by taking together heat fluxes in serie and in parallel (storage tank)
    % k1 and k2 in series
K1 = (k1*k2)/(k1+k2);
    % k3, k4 and k5 in series
K2 = (k3*k4*k5)/(k3*k4+k4*k5+k3*k5);

%% Reduced constants in front of Stefan-Boltzmann's Law (storage tank)
    % Radiation between the outer surface of the inner pvc pipe and the
    % inner surface of the outer pvc pipe
S1 = e_pvc * sigma * 2*pi*R1*L_pvc;
    % Radiation from the outer surface of the outer pvc pipe
S2 = e_pvc * sigma * 2*pi*Rmax*L_pvc;
    % Radiation from the air to the outer surface of the outer pvc pipe
S3 = e_air * sigma * 2*pi*Rmax*L_pvc;

%% Define the system for the heat flow out the storage vessel
% system definition
system = @(Q, T1, T2, T3, Tw) [K1*(Tw-T1) - Q, K2*(T1-T2) + S1*(T1^4-T2^4) - Q, k6*(T2-T3) - Q, k7*(T3-Tair) + S2*T3^4 - S3*Tair^4 - Q];
% initial guess for the solution
F0 = [10, 300, 300, 300];

%% Calculate thermal conductivities of each layer in W/K (solar collector)
    % 1. Convection between the aluminium plate and the air
k8 = ha_al * A;
    % 2. Convection between the air and the copper tube. I assume the
    % contact area covers 9/10 of the outer surface of the copper tube.
k9 = (9/10)*ha_cop * 2*pi*R_cop_o*L_cop;
    % 3. Conduction between the outer surface of the copper tube and the
    % inner surface for the part in contact with the air. I assume
    % the connection covers 9/10 of the outer surface of the copper tube
k10 =  (9/10)*2*pi*k_cop*L_cop/log(R_cop_o/R_cop_i);
    % 4. Conduction between the aluminium plate and the inner surface of
    % the copper tube. I assume the connection covers 1/10 of the outer
    % surface of the copper tube.
k11 = (1/10)*2*pi*k_cop*L_cop/log(R_cop_o/R_cop_i);
    % 5. Convection between the inner surface of the copper tube and the
    % water from the connection with the aluminium plate
k12 = (1/10)*hw_cop * 2*pi*R_cop_i*L_cop;
    % 6. Convection between the inner surface of the copper tube and the
    % water from the connection with the air
k13 = (9/10)*hw_cop * 2*pi*R_cop_i*L_cop;

%% Reduce the thermal conductivities by taking together heat fluxes in serie and parallel (solar collector)
    % k10 and k13 in series, 
K4 = (k10*k13)/(k10+k13);
    % k11 and k12 in series
K5 = (k11*k12)/(k11+k12);

%% Calculate thermal conductivities of each layer in W/K (polyurethane tube)
    % 1. Convection between the air and the polyurethane tube
k13 = ha_pol * 2*pi*R_pol_o*L_pol;
    % 2. Conduction through the polyurethane tube
k14 = 2*pi*k_pol*L_pol/log(R_pol_o/R_pol_i);
    % 3. Convection between polyurethane tube and the water
k15 = hw_pol * 2*pi*R_pol_i*L_pol;

%% Calculate the constants in front of Stefan-Boltzmann's Law (solar collector)
    % 1. Radiation from the metal plate
S4 = sigma * e_al * A;
    % 2. Radiation from the air inside the collector to the metal plate and
    % the outside air
S5 = sigma * e_air * A;
    % 3. Radiation from the air inside the collector to the copper tube
S6 = sigma * e_air * (9/10) * 2*pi*R_cop_o*L_cop;
    % 4. Radiation from the copper tube to the air inside the collector

%% Define the system for the heat from the air to the water (solar collector)
% system definition
system = @(Q, T1, Ta, Tw) [];
% initial guess for the solution
F0 = [10, 300, 300, 300];

%% Reduce the thermal conductivities by taking together heat fluxes in serie and parallel (polyurethane tube)
% k13, k14 and k15 in series
K6 = (k13*k14*k15)/(k13*k14+k14*k15+k13*k15);

save('variables.mat');
