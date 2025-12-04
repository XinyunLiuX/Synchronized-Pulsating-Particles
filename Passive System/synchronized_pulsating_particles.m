% 2025-11-07 by Xinyun Liu
% Sync Pulsating Particles.pptx
% Computer "Experiments" on Classical Fluids. I. Thermodynamical Properties
% of Lennard-Jones Molecules Physical Review 1967
% Consider a system of N particles, enclosed in a cube of side L. 
% with periodic boundary conditions interacting through a two-body potential of the Lennard-Jones type  

folder = "RES";
mkdir(folder)

N = 15^2;           % number of particles
rho = 1;          % density of particles rho = N/L^2, number of particles per area
L = sqrt(N/rho);    % domain width
l0 = sqrt(1/rho);   % approximate distance between two particles

% Cell Structures
rc = 2.5;                   % cut off distance of Lennard Jones potential
nc = floor(L/rc);           % number of cells in x or y direction
lc = L/nc;                  % length of the cell
p = struct('N', {N}, 'L', {L}, 'rc', {rc}, 'nc', {nc}, 'lc', {lc}, 'rho', {rho});

% t = 0. The positions are initially taken, in general, at the Nodes of the square lattice which has the desire density.
x0 = 0.5*l0:l0:(sqrt(N)-0.5)*l0; x0 = mod(x0, L);
y0 = 0.5*l0:l0:(sqrt(N)-0.5)*l0; y0 = mod(y0, L);
[X0, Y0] = meshgrid(x0, y0);
x0 = X0(:);
y0 = Y0(:);
[cellList, particlesPerCell] = buildCellList(x0, y0, p);

% Time Matching
h = 0.016; 
p.h = h;
h2 = h^2;
M = 5000;

% res.mat
numStepPerSave = 2;
X0 = zeros(floor(M/numStepPerSave) + 1, N);
Y0 = zeros(floor(M/numStepPerSave) + 1, N);
X0(1,:) = x0.';
Y0(1,:) = y0.';

% j = 0
u0 = 1E-6*randn(size(x0));
v0 = 1E-6*randn(size(x0));

% j = 1
[ax0, ay0] = accel(x0, y0, cellList, particlesPerCell, p);
x1 = x0 + u0*h + 0.5*ax0*h2; x1 = mod(x1, L);
y1 = y0 + v0*h + 0.5*ay0*h2; y1 = mod(y1, L);
[cellList, particlesPerCell] = buildCellList(x1, y1, p);

% j = 2 to M
saveCounter = 1; 
for j = 2:M
    [ax1, ay1] = accel(x1, y1, cellList, particlesPerCell, p);
    x2 = -x0 + 2*x1 + ax1*h2; x2 = mod(x2, L);
    y2 = -y0 + 2*y1 + ay1*h2; y2 = mod(y2, L);
    % updata particle positions and the cell list
    x0 = x1; y0 = y1;
    x1 = x2; y1 = y2;
    [cellList, particlesPerCell] = buildCellList(x1, y1, p);
    % Save every numStepPerSave steps
    if mod(j, numStepPerSave) == 0
        saveCounter = saveCounter + 1;
        X0(saveCounter, :) = x1.';
        Y0(saveCounter, :) = y1.';
    end
end

save(sprintf('RES/res_N=%d_rho=%.2f.mat', N, rho), 'X0', 'Y0', 'p');

%%
function [cellList, particlesPerCell] = buildCellList(x, y, p)
    ix = floor(x/p.lc) + 1;
    iy = floor((p.L - y)/p.lc) + 1;
    cellList = cell(p.nc,p.nc);         % creat cell structure 
    particlesPerCell = zeros(p.nc,p.nc);  
    for n = 1:p.N
        cellList{ix(n), iy(n)}(end+1) = n;
        particlesPerCell(ix(n), iy(n)) = particlesPerCell(ix(n), iy(n)) + 1;
    end
end