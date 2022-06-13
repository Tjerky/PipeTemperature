%% Function to find the heat flux out of the storage tank

function Q = Q_water_storage(Tw)
% load the necessary variables from the variables.mat file
load('variables.mat', 'system_storage', 'F0', 'K4', 'Tair');

% Calculate the total heat flow out of the storage tank
    % 1. Calculate the heat flow through the cylindrical side of the storage tank
        % adjust the system of equations into the form needed by fsolve
F = @(x) system_storage(x(1), x(2), x(3), x(4), Tw);
        % use fsolve to solve the system for Q, T1, T2 and T3
x = fsolve(F, F0);
        % extract Q (the heat flux) out of the answer array
Q_cyl = x(1);
    % 2. Calculate the heat flow out of the endcaps of the storage tank
Q_caps = K4 * (Tw-Tair);
    % 3. Add the two heat fluxes to obtain the total heat flux out of the storage tank
Q = Q_cyl + 2 * Q_caps;

end