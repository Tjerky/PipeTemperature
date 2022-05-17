%% Function to generate a matrix Q with as columns the q vectors in the matrix equations

% define the function
function Q = GenerateQ(T, Tw, beta, gamma, P, N)
    % transpose T, since it is somehow flipped (this should be fixed in a
    % better way)
    T = transpose(T);
    
    % generate a NxP matrix filled with zeros
    Q = 0*T;
    
    % fill the matrix
        % upper row with the heat flux from the water calculated with
        % newton's law of cooling.
        Q(2, :) = 2*beta*(Tw - T(2, :));
        % lower row with the heating term due to the adsorption of sunlight
    for j = 1:P
        if j <= P/2
            Q(N-1, j) = 2*gamma*sin(2*pi*(j-1)/P);
        end
    end