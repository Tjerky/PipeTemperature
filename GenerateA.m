%% Function to generate the A matrix in the matrix equation

% function definition
function A = GenerateA(N, alpha, beta)
    % generate a NxN matrix filled with zeros
    A = zeros(N, N);
    
    % fill middle rows
    for j = 2:N-1
        % 1-alpha on the diagonal
        A(j, j) = 1-alpha;
        % alpha on the left diagonal
        A(j, j-1) = alpha;
        % alpha on the right diagonal
        A(j, j+1) = alpha;
    
    % fill upper row
    A(1,1) = 2+alpha-2*beta;
    A(1,2) = alpha;
    
    % fill lower row
    A(N,N) = 2-alpha;
    A(N, N-1) = alpha;
end