%% finding the steady state of a pipe

% define the constants
    % dimensions
Rmax = 0.006; % m, outside radius of the tube
Rmin = 0.0049; % m, inside radius of the tube
N = 100; % number of gridpoints in the r direction
P = 100; % number of gridpoints in the theta direction

dr = (Rmax-Rmin)/(N-2);
dphi = 2*pi/P;
    
    % material properties
        %pipe
rho = 295.5; % kg/m^3, density of the pipe material
k = 398; % W/(m K), thermal conductivity of the pipe material
cp = 1470; % J/(kg K), specific heat of the pipe material
c = 5000; % W/(m^2 K), heat transfer coefficient pipe -> water
        %water
cw = 4180; % J/(kg K)
Tw = 295; % K
    % sunlight
I = 1000; % W/m^2, intensity of the sunlight

% define the constants in front
r = cat(2, [0], linspace(Rmin, Rmax, N-2), [0]); % vector of all r positions
phi = linspace(0, 2*pi, P+1); % vector of all phi positions +1 extra
phi = phi(1:end-1); % vector of all phi positions

[Phi, R] = meshgrid(phi, r); % meshgrid of all positions

    % define the constanst in a grid
alpha = R*dr*dphi^2./(2*R.^2*dphi^2+2*dr^2); % constant to simplify the formulas
beta = R.^2*dphi^2./(2*R.^2*dphi^2+2*dr^2); % constant to simplify the formulas
gamma = dr^2./(2*R.^2*dphi^2+2*dr^2); % constant to simplify the formulas
delta = R*dr^2*dphi^2./(2*k*(R.^2*dphi^2+dr^2)); % constant to simplify the formulas

% initialization of the temperature grid
    % generate a matrix filled with zeros
T = 300*ones(N, P);

% Gauss-seidel relaxation of the temperature distribution. This essentially
% means that we fill in the linear system a hundred times or so and hope that it
% converges.
for b = 1:10
    for s = 1:100000
        % loop the grid in the r-direction
        for j = 3:N-2
            % loop over the grid in the phi-direction
            for i = 2:P-1
                % calculate the temperature at T(j, i) given the current guess
                % of the temperatures around it.
                T(j, i) = (alpha(j, i) + beta(j, i))*T(j+1, i) + (beta(j, i) -alpha(j, i))*T(j-1, i) + gamma(j, i)*(T(j, i+1)+T(j, i-1));

            % boundary conditions
                % at i = 1, it loops around
            T(j, 1) = (alpha(j,1) + beta(j,1))*T(j+1, 1) + (beta(j,1) -alpha(j,1))*T(j-1, 1) + gamma(j, 1)*(T(j, 2)+T(j, P));
                % at i = P, it loops around
            T(j, P) = (alpha(j,P) + beta(j,P))*T(j+1, P) + (beta(j,P) -alpha(j,P))*T(j-1, P) + gamma(j, P)*(T(j, 1)+T(j, P-1));
            end
        end

        for i = 1:P
            % at j = 1, no-flux boundary
            T(1, i) = T(2, i);
            % at j = N, no-flux boundary
            T(N, i) = T(N-1, i);
        end

        for i = 2:P-1
            % at j = 2, extra heat generation term due to heat exchange with
            % the water
            T(2, i) = (1/((delta(2, i)*c)/dr+1))*((alpha(2, i) + beta(2, i))*T(3, i) + (beta(2, i) -alpha(2, i))*T(2, i) + gamma(2, i)*(T(2, i+1)+T(2, i-1))+delta(2, i)*c*Tw/dr);

            % at j = N-1, extra heat generation term due to adsorption of
            % sunlight
                % calculate the intensity of the sun at this point 
            I_ = -I * cos(phi(i));
            if I_ < 0
                I_ = 0;
            end
            T(N-1, i) = ((alpha(N-1, i) + beta(N-1, i))*T(N, i) + (beta(N-1, i) -alpha(N-1, i))*T(N-2, i) + gamma(N-1, i)*(T(N-1, i+1)+T(N-1, i-1))+delta(N-1, i)*I_/dr); 
        end

        % four super special cases, where i = 1 or P and j = 2 or N-1
        % j = 2 and i = 1
        T(2, 1) = (1/((delta(2, 1)*c)/dr+1))*((alpha(2, 1) + beta(2, 1))*T(3, 1) + (beta(2, 1) -alpha(2, 1))*T(2, 1) + gamma(2, 1)*(T(2, 2)+T(2, P))+delta(2, 1)*c*Tw/dr);
        % j = 2 and i = P
        T(2, P) = (1/((delta(2, P)*c)/dr+1))*((alpha(2, P) + beta(2, P))*T(3, P) + (beta(2, P) -alpha(2, P))*T(1, P) + gamma(2, P)*(T(2, 1)+T(2, P-1))+delta(2, 1)*c*Tw/dr);
        % j = N-1 and i = 1
        T(N-1, 1) = (alpha(N-1,1) + beta(N-1,1))*T(N, 1) + (beta(N-1,1) -alpha(N-1,1))*T(N-2, 1) + gamma(N-1, 1)*(T(N-1, 2)+T(N-1, P));
        % j = N-1 and i = 1
        T(N-1, P) = (alpha(N-1,P) + beta(N-1,P))*T(N, P) + (beta(N-1,P) -alpha(N-1,P))*T(N-2, P) + gamma(N-1, P)*(T(N-1, 1)+T(N-1, P-1));
    end
    fprintf('Completed %i of the 10000000 \n', b*100000) 
end
