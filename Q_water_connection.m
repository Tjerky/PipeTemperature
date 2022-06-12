%% Function to find the heat flux out of the polyurethane tubes

function Q = Q_water_connection(Tw)
load('variables.mat', 'K9', 'Tair')

Q = K9 * (Tw-Tair);
end