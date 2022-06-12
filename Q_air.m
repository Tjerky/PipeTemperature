%% Function to calculate the net heat flowing into the air inside the solar collector

function Q = Q_air(Ta, Tp, Tw)
load('variables.mat', 'k8','K4', 'K7', 'K8', 'Tair')

Q = k8*(Tp-Ta) - (K7+K8)*(Ta-Tair) - K4*(Ta-Tw);

end