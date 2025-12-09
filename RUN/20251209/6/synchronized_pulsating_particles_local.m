% 2025-11-07 by Xinyun Liu
% Sync Pulsating Particles.pptx
% Computer "Experiments" on Classical Fluids. I. Thermodynamical Properties
% of Lennard-Jones Molecules Physical Review 1967
% Consider a system of N particles, enclosed in a cube of side L. 
% with periodic boundary conditions interacting through a two-body potential of the Lennard-Jones type  
% x0, y0, x1, y1, x2, y2, are the unwrapped coordinate sets, used in update
% x0_w, y0_w, x1_w, x2_w  are the wrapped coordinates for interaction, and data saving 

clc; clear all; close all;
folder = "RES";
mkdir(folder)
addpath("../../../")

batch_omega = 2.1823*2;
batch_epsilon = 0.00;

[batch_omega, batch_epsilon] = meshgrid(batch_omega, batch_epsilon);

numsimulations = size(batch_omega(:),1);

for i = 1:numsimulations

% Geometry Parameters
N = 20;              % number of particles
rho = 2;             % density of particles rho = N/L^2, number of particles per area
L = 2^(1/6)*N;       % domain width
l0 = 2^(1/6);        % approximate distance between two particles

% Driving and Damping Parameters
beta = 0.;           % damping coefficient
omega   = batch_omega(i);           % forcing frequency
epsilon = batch_epsilon(i);         % forcing strength

% Celllist  Structures
rc = 1.5;                   % cut off distance of Lennard Jones potential
nc = floor(L/rc);           % number of cells in x or y direction
lc = L/nc;                  % length of the cell

% Time resolution
h = min(0.001, 1/(2*pi*omega)/50)/20; 
h2 = h^2;
betah_2 = beta*h/2;
M = 100000;

% Parameter Container p
p = struct('beta', {beta}, 'epsilon', {epsilon}, 'omega', {omega}, 'N', {N}, 'L', {L}, 'rc', {rc}, 'nc', {nc}, 'lc', {lc}, 'rho', {rho}, 'h', {h}, 'M', {M});


% Result Container
numStepPerSave = 10;
T = zeros(floor(M/numStepPerSave) + 1, 1);
X = zeros(floor(M/numStepPerSave) + 1, N);
Y = zeros(floor(M/numStepPerSave) + 1, N);
U = zeros(floor(M/numStepPerSave) + 1, 1);
K = zeros(floor(M/numStepPerSave) + 1, 1);
E = zeros(floor(M/numStepPerSave) + 1, 1);

%% Initial Condition j = 0
% The positions are initially taken, in general, at the Nodes of the square lattice which has the desire density
p.t = 0; 
perb = 4E-1*ones(N,1); perb(1:2:end) = -perb(1:2:end);
x0 =  (0:N-1).'*2^(1/6) + 2^(1/6)/2 + perb; 
y0 = zeros(N,1) + L/2; 
u0 = 0E-3*randn(size(x0));
v0 = 0E-3*randn(size(x0));

% CellList for Interaction 
x0_w = mod(x0, L); y0_w = mod(y0, L);
[cellList, particlesPerCell] = buildCellList(x0_w, y0_w, p);
[ax0, ay0, U0] = accel_LJ(x0_w, y0_w, cellList, particlesPerCell, p);      % acceleration of Lennard-Jones part
ax0 = ax0 - beta*u0;  ay0 = ay0 - beta*v0;                                 % acceleration including damping

% Update: j = 1
p.t = p.t + h;
x1 = x0 + u0*h + 0.5*ax0*h2;                      
y1 = y0 + v0*h + 0.5*ay0*h2; 

% First Data Saving
T(1) = p.t;
X(1,:) = x0_w.';
Y(1,:) = y0_w.';
U(1) = U0;
K(1) = 1/2*sum(u0.^2 + v0.^2);
E(1) = U(1) + K(1);

% j = 2 to M
saveCounter = 1; 
for j = 2:M

    % CellList for Interaction 
    x1_w = mod(x1, L); y1_w = mod(y1, L);
    [cellList, particlesPerCell] = buildCellList(x1_w, y1_w, p);
    [ax1, ay1, U1] = accel_LJ(x1_w, y1_w, cellList, particlesPerCell, p);
    
    % Update
    p.t = p.t + h;
    x2 = ((-1 + betah_2)*x0 + 2*x1 + ax1*h2)/(1 + betah_2); 
    y2 = ((-1 + betah_2)*y0 + 2*y1 + ay1*h2)/(1 + betah_2); 

    % Save every numStepPerSave steps
    if mod(j, numStepPerSave) == 0
        saveCounter = saveCounter + 1;
        T(saveCounter) = p.t;
        X(saveCounter, :) = x1_w.';
        Y(saveCounter, :) = y1_w.';
        U(saveCounter) = U1;
        K(saveCounter) = 1/2*sum( ((x2-x0)/(2*h)).^2 + ((y2-y0)/(2*h)).^2 );
        E(saveCounter) = U(saveCounter) + K(saveCounter);
    end

    % updata particle positions 
    x0 = x1;    y0 = y1;
    x1 = x2;    y1 = y2;

end

fileName = fullfile(folder, sprintf('res_N=%d_rho=%.2f_beta=%.1f_epsilon=%.2f_omega=%.2f.mat', N, rho, beta, epsilon, omega))
parsave(fileName, T, X, Y, U, K, E, p);
disp(['Finish:' fileName])

end


%% Build CellList Func
function [cellList, particlesPerCell] = buildCellList(x, y, p)
    ix = floor(x/p.lc) + 1;
    iy = floor(y/p.lc) + 1; iy = p.nc - iy + 1;
    cellList = cell(p.nc,p.nc);         % creat cell structure 
    particlesPerCell = zeros(p.nc,p.nc);  
    for n = 1:p.N
        cellList{iy(n), ix(n)}(end+1) = n;
        particlesPerCell(iy(n), ix(n)) = particlesPerCell(iy(n), ix(n)) + 1;  % note the position exchange of ix iy
    end
end

function parsave(fileName, T, X, Y, U, K, E, p)
    save(fileName, 'T', 'X', 'Y', 'U', 'K', 'E', 'p');
end

