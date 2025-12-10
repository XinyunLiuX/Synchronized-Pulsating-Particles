function dXdt = synchronized_pulsating_particles(t, X ,p)
    p.t = t;
    N = p.N;
    
    dXdt = zeros(4*N,1);
    
    x = X(1:N);
    y = X(N+1:2*N);
    u = X(2*N+1:3*N);
    v = X(3*N+1:4*N);
    
    % CellList for Interaction 
    x_w = mod(x, p.L); y_w = mod(y, p.L);
    [cellList, particlesPerCell] = buildCellList(x_w, y_w, p);
    [ax, ay, ~] = accel_LJ(x_w, y_w, cellList, particlesPerCell, p);      % acceleration of Lennard-Jones part
    
    dXdt(1:N) = u;
    dXdt(N+1:2*N) = v;
    dXdt(2*N+1:3*N) = -p.beta*u + ax;
    dXdt(3*N+1:4*N) = -p.beta*v + ay;

end