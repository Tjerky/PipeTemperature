%% Function to calculate the net heat flowing into the water inside the solar collector

function Q = Q_water_collector(Ta, Tp, Tw)
load('variables.mat', 'K4', 'K5')

Q = K4*(Ta-Tw) + K5*(Tp-Tw);

end