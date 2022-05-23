%% Function to find the heat flux out of the polyurethane tubes

function Q = Q_polyurethane_tube(Tw)
load('variables.mat', 'K6', 'Tair')

Q = K6 * (Tw-Tair);
end