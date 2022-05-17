%% function to animate the outcome of a simulation

% function definition
function plotTemperatureDistribution(T, R, Phi)
    % convert the polar coordinates
    [X, Y] = pol2cart(Phi(2:end-1, :),R(2:end-1, :));

    mesh(X, Y, T(2:end-1, :))
end