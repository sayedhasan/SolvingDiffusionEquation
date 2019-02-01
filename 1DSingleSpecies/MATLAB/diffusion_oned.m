%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Single Species 1D diffusion equation using forward marching in time.
%  A simple TOY problem.
%
%  Equation:
% 
%       dP/dt = D*d^2P/dx^2
%          With Neumann boundary condition at x=0 and x=x1
%          And with initial profile: P(x,0) = P0;  for  0<= x <= 0.2
%             
%  Use: Just run the program by typing it's name
%  You can change the input parameters. All values are in SI unit.
%
%  Sayed Hasan 1/31/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;   % clear the workspace



%--------------------------------------------------------
%--- INPUT PARAMETERS -----------------------------------
%--------------------------------------------------------
T0 = 0;      % initial time [sec]
T1 = 2;      % final time   [sec]
dt = 1e-2;   % time step    [sec]
D  = 0.05;   % Diffusion constant [m2/sec]
x0 = 0;      % initial x-location [m]
x1 = 1;      % final x-location   [m]
NX = 100;    % number of x-points
P0 = 10;     % initial concentration [#/m3]
X0 = 0.2;    % initial concentrention extent [m]
%--- Note ---
%   P(x,0) = P0;  for    0<= x <= x0
%---------------------------------------------------------




%------------------ Program Implementation --------------
clf;     % clear figure
clc;     % clear the screen


%--- calculate number of time steps
NT = (T1 - T0)/dt; 
t  = T0:dt:T1;  

%--- position vector
x  = linspace(x0, x1, NX);
dx = x(2) - x(1);

%--- discretize the diffusion operator using center derivative
%    M = D * d2/dx2
%
alpha = dt*D/(dx*dx);
M = diag(2*alpha*ones(NX,1)) - diag(alpha*ones(NX-1,1),1) - ...
    diag(alpha*ones(NX-1,1),-1);

%--- apply Neumann B.C. at both x0 and x1, meaning we are assuming Conc.
%    gradient must be zero at the both ends -- a reasonable thing to do.
M(1,1)   = alpha;  % this is for node x0
M(NX,NX) = alpha;  % this is for node x1

%--- now add the contribution from previous time step. This basically adds
%    1 to the diagonal of the Diffusion operator M.
M = M + eye(NX);

%--- define the initial concentration.
P = zeros(NX, NT);          %-- row is position, column is time
P(1:round(X0/dx), 1) = P0;  %-- P(x0, 0) = P0;    0<= x <= x0

%--- start the time iteration now
figure(1);
set(gcf, 'color', 'w');
set(gca, 'fontsize', 14);

for it = 1:NT-1
    %--- solve the system
    Pini = P(:,it);    %-- load solution at t as Pini
    Pfin = M \ Pini;   %-- solve conc. for time t+dt
    P(:,it+1) = Pfin;  %-- save the solution at column it+1
    
    %--- plot the solution
    plot(x, P(:,it), 'linew', 2);  %-- basic line plot
    set(gca, 'ylim', [0, 1.1*P0]); %-- set fixed y-scale limit
    xlabel('Position X (m)', 'fontsize', 14); %-- xlabel
    ylabel('Concentration (#/m^3)', 'fontsize', 14); %-- ylabel
    title(sprintf('Profile at time = %2.1e Sec', t(it+1))); %-- title
    drawnow;
end
