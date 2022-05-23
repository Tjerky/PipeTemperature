%% Function to find the heat flux from the metal plate to the water

function Q = Q_plate_to_water(Tw, Tp)
load('variables.mat', 'K5')

Q = K5*(Tp - Tw);
end