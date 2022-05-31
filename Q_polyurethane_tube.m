%% Function to find the heat flux out of the polyurethane tubes

function Q = Q_polyurethane_tube(Tw)
load('variables.mat', 'K8', 'Tair')

Q = K8 * (Tw-Tair);
end