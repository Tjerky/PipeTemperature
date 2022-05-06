%% Function to generate a matrix Q with as columns the q vectors in the matrix equations

% define the function
function Q = GenerateQ(P, N, T, beta, I)
    % generate a NxP matrix filled with zeros
    Q = zeros(N, P);
    
    % fill the matrix
    for j = 1:P
        % upper row with the cooling term due to the flow of water through
        % the pipe
        Q(1, j) = 2*beta*T;
        % lower row with the heating term due to the adsorption of sunlight
        if j <= P/2
            Q(N, j) = I * sin(2*pi*(j-1)/P);
        end
    end
    
    % fil