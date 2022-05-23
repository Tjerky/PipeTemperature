%% Function to find the heat flux out of the storage tank

function Q = Q_storage_tank(Tw)
load('variables.mat', 'system', 'F0');

F = @(x) system(x(1), x(2), x(3), x(4), Tw);

x = fsolve(F, F0);

Q = x(1);

end