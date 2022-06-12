%% Function to calculate the net heat flowing into the aluminum plate

function Q = Q_plate(Tp, Ta, Tw)
load('variables.mat', 'K6', 'K7', 'k13', 'Tair', 'I', 'A')

% I*A: adsorption of sunlight, K6*(Tp-Tair): conduction through the tempex
% and trespa, k8*(Tp-Ta): convection to the air, K5*(Tp-Tw): conduction to
% the water through the copper pipes.
Q = I*A - K7*(Tp - Tair) - k13*(Tp-Ta) - K6*(Tp-Tw);

end