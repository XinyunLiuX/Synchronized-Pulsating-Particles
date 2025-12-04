function [AX, AY, U] = accel_LJ(x, y, cellList, numParCell, p)
% Sync Pulsating Particles.pptx
% Calculate Acceleration based on current position x and y
% Domain is divided by ncxnc cells for acceleration

% input: 
% x, y:           vectors (N x 1) of coordinates for particle i = 1,2, ..., N
% CELL:           cell (nc x nc), each contains the particle id within the cell 
% numParCell:     cell (nc x nc), each contains the number of particle within the cell
% p:              parameter container

% output;
% AX, AY: vectors (N x 1) of x (y)-acceleration for particle i = 1,2, ..., N


N = p.N;            % num of particles 
rc = p.rc;          % cutoff interaction distance
nc = p.nc;          % num of cells in x,y direction
AX = zeros(N,1);    % x-acceleration 
AY = zeros(N,1);    % y-acceleration
U  = 0;
for i = 1:nc
    for j = 1:nc    % loop over all cells
        n0 = numParCell(i,j);               % number of particles in the current cell
        for k = 1:n0                        % loop over particles in the current cell
            id0 = cellList{i,j}(k);             % focus on the current particle id0

            % Intracellular Interaction
            c1 = i; c2 = j;
            for l = k+1:n0                  % loop over interacting particles in the same cell that are futher down the list
               id1 = cellList{c1,c2}(l);        % for interacting particle id1
               [ax, ay, Uij] = LennardJones(x(id0), x(id1), y(id0), y(id1), p);     
               AX(id0) = AX(id0) + ax; 
               AY(id0) = AY(id0) + ay;
               % Newton's third law
               AX(id1) = AX(id1) - ax;
               AY(id1) = AY(id1) - ay;
               % Potential Energy
               U = U + Uij;
           end

           % Intercellular Interaction: 1: East 2: Southwest 3. South 4. Southeast
           for neighbor = 1:4                   % loop over 4 neighboring cells
               switch neighbor
                   case 1                       
                       c1 = i; c2 = j+1;        % east
                       if j == nc
                           c2 = 1;      
                       end
                   case 2
                       c1 = i+1; c2 = j-1;      % southwest
                       if i == nc
                           c1 = 1;
                       end
                       if j == 1
                           c2 = nc;
                       end
                   case 3
                       c1 = i+1; c2 = j;        % south
                       if i == nc
                           c1 = 1; 
                       end
                   case 4
                       c1 = i+1; c2 = j+1;      % southeast
                       if i == nc
                           c1 = 1;
                       end
                       if j == nc
                           c2 = 1;
                       end
               end
               n1 = numParCell(c1,c2);          % number of particles in that neighboring cell
               for l = 1:n1                         % loop over all particles in that neighboring cell
                   id1 = cellList{c1,c2}(l);            % for interacting particle id1
                   [ax, ay, Uij] = LennardJones(x(id0), x(id1), y(id0), y(id1), p);     
                   AX(id0) = AX(id0) + ax; 
                   AY(id0) = AY(id0) + ay;
                   % Newton's third law
                   AX(id1) = AX(id1) - ax;
                   AY(id1) = AY(id1) - ay;
                   % Potential Energy
                   U = U + Uij;
               end
           end


        end

    end
end




         
end