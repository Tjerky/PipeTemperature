%% define the constants
    % dimensions of the storage vessel
R1 = 0.025; % m, inside radius of the inner pvc tube
R2 = 0.0261; % m, outside radius of the inner pvc tube
R3 = 0.0539; % m, inside radius of the outer pvc tube
R4 = 0.055; % m, outside radius of the outer pvc tube
R5 = 0.1; % m, outside radius of the polyethylene foam
L_pvc = 0.55; % m, length of the pvc pipes
d_caps = 2*0.0022; % m, thickness of the two endcaps glued together

    % dimensions of the polyurethane tubing
L_pol = 6; % m, length of the polyurethane tube
R_pol_o = 0.006; % m, outside radius of the polyurethane tube
R_pol_i = 0.004; % m, inside radius of the polyurethane tube

    % dimensions of the solar collector
        % aluminium plate
A = 1.1; % m^2, area of the aluminium plate
d_al = 0.002; % m, thickness of the aluminium plate
m_al = 5.883; % kg, mass of the aluminium plate
        % glass
d_glass = 0.004; % m, thickness of the glass panel
        % copper pipe
L_cop = 5; % m, length of the copper pipe
R_cop_o = 0.006; % m, outside radius of the copper pipe
R_cop_i = 0.00545; % m, inside radius of the copper pipe
c_contact = 0.01; % unitless, ratio between the area of the pipe in contact with the aluminium plate and the total outside area of the pipe.
        % tempex isolation plate
d_tempex = 0.05; % m, thickness of the tempex isolation
        % trespa plates
d_trespa = 0.006; % m, thickness of the trespa plates
        % pinewood enclosure
d_pin = 0.027; % m, thickness of the pinewoord enclosure
A_pin = 0.017; % m^2, area of the pipewood in contact with the air
        % air inside the solar collector
m_air = 0.1225; % kg, mass of the air in the solar collector
        % rest
H = 0.007; % m, distance between aluminium plate and the glass plate
    
    % material properties
        % pvc pipe
k_pvc = 0.19; % W/(m K), thermal conductivity of the pipe material (pvc), source: Project manual
hw_pvc = 30; % W/(m^2 K), heat transfer coefficient pipe -> water, source: 
hv_pvc = 0.1; % W/(m^2 K), heat transfer coefficient pipe -> partial vacuum, source: 
e_pvc = 0.05; % unitless, emissivity of pvc, source: https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
        % vacuum
k_v = 1; % W/(m K), thermal conductivity of the vaccuum
        %water
cw = 4180; % J/(kg K)
mw = 1; % kg, mass of the water in the system
hw_cop = 300; % W/(m^2 K), heat transfer coefficient copper -> water
hw_pol = 300; % W/(m^2 K), heat transfer coefficient polyurethane -> water
        % air
rho_air = 1.225; % kg/m^3, density of air
c_air = 718; % J/(kg K), specific heat of air
ha_foam = 10; % W/(m^2 K), heat transfer coefficient polyethylene foam -> air
ha_cop = 10; % W/(m^2 K), heat transfer coefficient copper -> air
ha_al = 10; % W/(m^2 K), heat transfer coefficient aluminium -> air
ha_pol = 10; % W/(m^2 K), heat transfer coefficient polyurethane -> air
ha_tres = 10; % W/(m^2 K), heat transfer coefficient trespa -> air
ha_glass = 10; % W/(m^2 K), heat transfer coefficient glass -> air
ha_pin = 10; % W/(m^2 K), heat transfer coefficient pinewood -> air
ha_tempex = 10; % W/(m^2 K), heat transfer coefficient tempex -> air
e_air = 0.8; % unitless, emissivity of the air
Tair = 295; % K, temperature of the air
        % aluminium
c_al = 880; % J/(kg K), specific heat of aluminium 
k_al = 237; % W/(m K) thermal conductivity of aluminium
e_al = 0.1; % unitless, emissivity of the aluminium plate
        % copper
k_cop = 402; % W/(m K), thermal conductivity of copper
e_cop = 0.03; % unitless, emissivity of the copper tube
c_cop = 377; % J/(kg K)
        % polyurethane
k_pol = 0.13; % W/(m K), thermal conductivity of polyurethane
c_pol = 1800; % J/(kg K), specific heat polyurethane
        % tempex
k_tempex = 0.03; % W/(m K), thermal conductivity of tempex
c_tempex = 1100; % J(kg K), specific heat of tempex
e_glass = 0.88; % unitless, thermal emissivity of glass
        % glass
k_glass = 0.96; % W/(m K), thermal conductivity of glass
c_glass = 780; % J/(kg K), specific heat glass
t_glass = 0.9; % unitless, transmittance of the glass
        % pinewood
c_pin = 2300; % J/(kg K), specific heat of pinewood
k_pin = 0.12; % W/(m K), thermal conductivity of pinewood
        % trespa
k_trespa = 0.3; % W/(m K), thermal conductivity of trespa
        % polyethylene foam
k_foam = 0.04; % W/(m K), thermal conductivity of the polyethylene foam
    
    % sun
I = 1000; % W/m^2, intensity of the sunlight
    % Univeral constants
sigma = 5.670374419 * 10^-8; % W(m^2 K^4), Stefan-Boltzmann constant

    % simulation
runtime = 120*60; % s, duration of the simulated period
dt = 1; % s, duration of one timestep

P = ceil(runtime/dt); % unitless, amount of timesteps

    % initial condition
T_initial = 295; % K, initial temperature of the water

%% Calculate thermal conductivities of each layer in W/K (storage tank)
    % 1. Interface between the water and the inner surface of the inner pvc
    % pipe
k1 = hw_pvc * 2*pi*R1*L_pvc;
    % 2. Conduction through the inner pvc pipe
k2 = 2*pi*k_pvc*L_pvc/log(R2/R1);
    % 3. Convection between the partial vacuum and the outer surface of the
    % inner pvc tube
k3 = hv_pvc * 2*pi*R2*L_pvc;
    % 4. Conduction through the partial vacuum
k4 = 2*pi*k_v*L_pvc/log(R3/R2);
    % 5. Convention between the partial vacuum and the inner surface of the
    % outer pvc tube
k5 = hv_pvc * 2*pi*R3*L_pvc;
    % 6. Conduction through the outer pvc pipe
k6 = 2*pi*k_pvc*L_pvc/log(R4/R3);
    % 7. Conduction through the polyethylene foam wrapping
k7 = 2*pi*k_foam*L_pvc/log(R5/R4);
    % 8. Convection between the air and the outer surface of the
    % polyethylene foam wrapping
k8 = ha_foam * 2*pi*R5*L_pvc;
    % 9. Convection from the water to the pvc endcaps
k9 = hw_pvc * pi*R1^2;
    % 10. Conduction through both pvc endcaps
k10 = k_pvc * pi*R1^2/d_caps;
    % 11. Conduction through the tempex plate
k11 = k_tempex * pi*R1^2/d_tempex;
    % 12. Convection from the tempex to the air
k12 = ha_tempex * pi*R1^2;

%% Reduce the conductivities by taking together heat fluxes in serie and in parallel (storage tank)
    % k1 and k2 in series, heat transfer to through the inner pvc tube
K1 = (k1*k2)/(k1+k2);
    % k3, k4 and k5 in series, heat transfer through the vacuum
K2 = (k3*k4*k5)/(k3*k4+k4*k5+k3*k5);
    % k6 and k7 in series, heat transfer trough the outer pvc tube and the
    % polyethylene wrapping.
K3 = (k6*k7)/(k6+k7);
    % k9, k10, k11 and k12 in series
K4 = (k9*k10*k11*k12)/(k9*k10*k11 + k9*k10*k12 + k9*k11*k12 + k10*k11*k12);
%% Reduced constants in front of Stefan-Boltzmann's Law (storage tank)
    % Radiation between the outer surface of the inner pvc pipe and the
    % inner surface of the outer pvc pipe
S1 = (e_pvc/(2-e_pvc)) * sigma * 2*pi*R2*L_pvc;
    % Radiation from the outer surface of the outer pvc pipe
S2 = (e_pvc/(2-e_pvc)) * sigma * 2*pi*R4*L_pvc;
    % Radiation from the air to the outer surface of the outer pvc pipe
S3 = e_air * sigma * 2*pi*R4*L_pvc;

%% Define the system for the heat flow out the storage vessel
% system definition
system_storage = @(Q, T1, T2, T3, Tw) [K1*(Tw-T1) - Q, K2*(T1-T2) + S1*(T1^4-T2^4) - Q, K3*(T2-T3) - Q, k8*(T3-Tair) + S2*T3^4 - S3*Tair^4 - Q];
% initial guess for the solution
F0 = [10, 300, 300, 300];

%% Calculate thermal conductivities of each layer in W/K (solar collector)
    % 1. Convection between the aluminium plate and the air
k13 = ha_al * A;
    % 2. Convection between the air and the copper tube.
k14 = (1-c_contact)*ha_cop * 2*pi*R_cop_o*L_cop;
    % 3. Conduction between the outer surface of the copper tube and the
    % inner surface for the part in contact with the air. I assume
    % the connection covers 9/10 of the outer surface of the copper tube
k15 =  (1-c_contact)*2*pi*k_cop*L_cop/log(R_cop_o/R_cop_i);
    % 4. Conduction between the aluminium plate and the inner surface of
    % the copper tube. I assume the connection covers 1/10 of the outer
    % surface of the copper tube.
k16 = c_contact*2*pi*k_cop*L_cop/log(R_cop_o/R_cop_i);
    % 5. Convection between the inner surface of the copper tube and the
    % water from the connection with the aluminium plate
k17 = c_contact*hw_cop * 2*pi*R_cop_i*L_cop;
    % 6. Convection between the inner surface of the copper tube and the
    % water from the connection with the air
k18 = (1-c_contact)*hw_cop * 2*pi*R_cop_i*L_cop;
    % 7. Conduction through the tempex isolation plate
k19 = k_tempex * A/d_tempex;
    % 8. Conduction through the trespa plate
k20 = k_trespa*A/d_trespa;
    % 9. Convection from the trespa plate to the air
k21 = ha_tres * A;
    % 10. Conduction through the glass
k22 = k_glass * A /d_glass;
    % 11. Convection from the air to glass
k23 = ha_glass * A;
    % 12. Convection from the air to the pinewood enclosure
k24 = ha_pin * A_pin;
    % 13. Conduction through the pinewood
k25 = k_pin * A_pin/d_pin;

%% Reduce the thermal conductivities by taking together heat fluxes in serie and parallel (solar collector)
    % k14, k15 and k18 in series, heat transfer from the air to the water
    % through the copper pipe
K5 = (k14*k15*k18)/(k14*k15+k14*k18+k15*k18);
    % k16 and k17 in series, heat transfer from the aluminium plate to the
    % water through the copper pipe
K6 = (k16*k17)/(k16+k17);
    % k14, k15 and k16 in series, heat transfer through the tempex and the
    % trespa
K7 = (k19*k20*k21)/(k19*k20+k20*k21+k19*k21);
    % k22, k23 in series, heat transfer through the glass panel to the
    % outside air
K8 = (k22*k23)/(k22+k23);
    % k24, k25 and k24 in series, heat transfer through the pinewood
    % enclosure
K9 = (k25*k24^2)/(2*k25*k24+k24^2);

%% Calculate the constants in front of Stefan-Boltzmann's Law (solar collector)
    % 1. Radiation from the metal plate to the glass
S4 = sigma * e_al * e_glass * A;
    % 2. Radiation from the air inside the collector to the metal plate and
    % the outside air
S5 = sigma * e_air * A;
    % 3. Radiation from the air inside the collector to the copper tube
S6 = sigma * e_air * (9/10) * 2*pi*R_cop_o*L_cop;
    % 4. Radiation from the copper tube
S7 = sigma * e_cop * (9/10) * 2*pi*R_cop_o*L_cop;

%% Define the system for the heat transfer through the glass panel
% system definition
system_glass = @(Q1, Q2, T1, Ta, Tp) [S4*(Tp^4-T1^4) - Q1, k23*(Ta-T1) - Q2, K8*(T1-Tair) - Q1-Q2];
% initial guess for the solution
F1 = [10, 10, 300];

%% Calculate thermal conductivities of each layer in W/K (polyurethane tube)
    % 1. Convection between the air and the polyurethane tube
k26 = ha_pol * 2*pi*R_pol_o*L_pol;
    % 2. Conduction through the polyurethane tube
k27 = 2*pi*k_pol*L_pol/log(R_pol_o/R_pol_i);
    % 3. Convection between polyurethane tube and the water
k28 = hw_pol * 2*pi*R_pol_i*L_pol;

%% Reduce the thermal conductivities by taking together heat fluxes in serie and parallel (polyurethane tube)
% k26, k27 and k28 in series
K10 = (k26*k27*k28)/(k26*k27+k26*k28+k27*k28);

save('variables.mat');
