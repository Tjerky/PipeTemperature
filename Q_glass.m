%% Function to find the heat fluxes to and through the glass panel
% 1) heat flux from the metal plate to the glass panel
% 2) heat flux from the air inside the collector to the glass panel

function [Q_plate_glass, Q_air_glass] = Q_glass(Tp, Ta)
% load the necessary variables from the variables.mat file
load('variables.mat', 'system_glass', 'F1');

% adjust the system of equations into the form needed by fsolve
F = @(x) system_glass(x(1), x(2), x(3), x(4), Ta, Tp);
% use fsolve to solve the system
x = fsolve(F, F1);

% Extract the heat flows from the solution of the system
Q_plate_glass = x(1);
Q_air_glass = x(2);

end