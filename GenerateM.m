%% Function to generate the M matrix in the matrix equation

%function definition
function M = GenerateM(N, alpha)
    % generate a NxN matrix filled with zeros
    M = zeros(N, N);
    
    % fill middle rows
    for j = 2:N-1
        % 2 + 2 * alpha on the diagonal
        M(j, j) = 2*alpha + 2;
        % -alpha on the left diagonal
        M(j, j-1) = -1*alpha;
        % -alpha on the right diagonal
        M(j, j+1) = -1*alpha;
    end
    
    %fill upper row
    M(1, 1) = 2-alpha;
    M(1, 2) = alpha;
    % fill lower row
    M(N,N) = 2+alpha;
    M(N, N-1) = -1*alpha;
end