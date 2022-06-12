%% Function to find the heat flux out of the storage tank

function Q = Q_water_storage(Tw)
load('variables.mat', 'system_storage', 'F0');

F = @(x) system_storage(x(1), x(2), x(3), x(4), Tw);

x = fsolve(F, F0);

Q = x(1);

end