%% Function to calculate the net heat flowing into the aluminum plate

function Q = Q_plate(Tp, Ta, Tw)
load('variables.mat', 'k8', 'K5', 'K6', 'Tair', 'I', 'A')

% I*A: adsorption of sunlight, K6*(Tp-Tair): conduction through the tempex
% and trespa, k8*(Tp-Ta): convection to the air, K5*(Tp-Tw): conduction to
% the water through the copper pipes.
Q = I*A - K6*(Tp - Tair) - k8*(Tp-Ta) - K5*(Tp-Tw);

end