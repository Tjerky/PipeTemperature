%% Function to calculate the net heat flowing into the air inside the solar collector

function Q = Q_air(Ta, Tp, Tw)
load('variables.mat', 'k8','K9', 'K5', 'Tair')

Q = k8*(Tp-Ta) - K9*(Ta-Tair) - K5*(Ta-Tw);

end