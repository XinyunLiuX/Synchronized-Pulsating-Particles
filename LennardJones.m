function [ax, ay, Uij] = LennardJones(t, x0, x1, y0, y1, p)
    rc = p.rc;
    L = p.L;

    dx = x0 - x1;     dx = dx - L*round(dx/L);
    dy = y0 - y1;     dy = dy - L*round(dy/L);

    e = 1 + p.epsilon*sin(p.omega*t);
    e2 = e^2;

    r2 = dx^2 + dy^2;
    r2 = r2/e2;

    if r2 > rc^2
        ax = 0; 
        ay = 0;
        Uij = 0;
        return
    end
    ax = e2^(-1)*(r2^(-7) - 0.5*r2^(-4))*dx;
    ay = e2^(-1)*(r2^(-7) - 0.5*r2^(-4))*dy;
    Uij = 1/12*(r2^(-6) - r2^(-3));

end