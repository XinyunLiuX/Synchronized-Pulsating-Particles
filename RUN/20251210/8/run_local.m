% 2025-12-10 by Xinyun Liu
% ODE45 solver is included
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

batch_omega = 2.6*2;
batch_epsilon = 0.07;

[batch_omega, batch_epsilon] = meshgrid(batch_omega, batch_epsilon);

numsimulations = size(batch_omega(:),1);

for i = 1:numsimulations

% Geometry Parameters
N = 20;              % number of particles
rho = 2;             % density of particles rho = N/L^2, number of particles per area
L = 2^(1/6)*N;       % domain width
l0 = 2^(1/6);        % approximate distance between two particles

% Driving and Damping Parameters
beta = 1.4;           % damping coefficient
omega   = batch_omega(i);           % forcing frequency
epsilon = batch_epsilon(i);         % forcing strength

% Celllist  Structures
rc = 1.5;                   % cut off distance of Lennard Jones potential
nc = floor(L/rc);           % number of cells in x or y direction
lc = L/nc;                  % length of the cell

% Time resolution
dt = 0.1;
tspan = 0:dt:3000;
M = length(tspan);

% Parameter Container p
p = struct('beta', {beta}, 'epsilon', {epsilon}, 'omega', {omega}, 'N', {N}, 'M', {M}, 'L', {L}, 'rc', {rc}, 'nc', {nc}, 'lc', {lc}, 'rho', {rho}, 'dt', {dt});

%% Initial Condition 

% The positions are initially taken, in general, at the Nodes of the square lattice which has the desire density
perb = 1E-5*ones(N,1); perb(1:2:end) = -perb(1:2:end);
x0 =  (0:N-1).'*2^(1/6) + 2^(1/6)/2 + perb; 
y0 = zeros(N,1) + L/2; 
u0 = 0E-3*randn(size(x0));
v0 = 0E-3*randn(size(x0));

X0 = [x0; y0; u0; v0];

% Numerical tolerances
RelTol = 1E-6;
AbsTol = 1E-6;

opts = odeset('RelTol',RelTol, 'AbsTol',AbsTol);

[T, X] = ode45(@(t,X) synchronized_pulsating_particles(t, X ,p), tspan, X0, opts);

x = X(:,1:N);
y = X(:,N+1:2*N);
u = X(:,2*N+1:3*N);
v = X(:,3*N+1:4*N);

fileName = fullfile(folder, sprintf('res_N=%d_rho=%.2f_beta=%.1f_epsilon=%.2f_omega=%.2f.mat', N, rho, beta, epsilon, omega))
parsave(fileName, T, x, y, u, v, p);
disp(['Finish:' fileName])

end


function parsave(fileName, T, x, y, u, v, p)
    save(fileName, 'T', 'x', 'y', 'u', 'v','p');
end

