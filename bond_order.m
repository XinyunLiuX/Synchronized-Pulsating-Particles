function [Psi4, Psi6] = bond_order(x, y)
% BOND_ORDER  Compute bond-orientational order parameters in 2D

% ------------------------------------------------------------
% Setup
% ------------------------------------------------------------
pos = [x(:), y(:)];
N   = size(pos,1);


psi4 = zeros(N,1);
psi6 = zeros(N,1);
Nb   = zeros(N,1);

a_nn = zeros(N,1);

for i = 1:N
    dx = pos(:,1) - pos(i,1);
    dy = pos(:,2) - pos(i,2);
    r  = hypot(dx,dy);
    r(i) = inf;                 % exclude self
    a_nn(i) = min(r);
end

a = median(a_nn);               % robust against boundaries
rc = 1.3 * a;                   % nearest-neighbor cutoff

% ------------------------------------------------------------
% Neighbor loop (metric, symmetry-faithful)
% ------------------------------------------------------------
for i = 1:N
    for j = i+1:N
        dx = pos(j,1) - pos(i,1);
        dy = pos(j,2) - pos(i,2);
        r  = hypot(dx, dy);

        if r < rc && r > 0
            theta = atan2(dy, dx);

            % add bond contribution

            psi4(i) = psi4(i) + exp(1i*4*theta);
            psi4(j) = psi4(j) + exp(1i*4*(theta+pi));

            psi6(i) = psi6(i) + exp(1i*6*theta);
            psi6(j) = psi6(j) + exp(1i*6*(theta+pi));

            Nb(i) = Nb(i) + 1;
            Nb(j) = Nb(j) + 1;
        end
    end
end

% ------------------------------------------------------------
% Normalize local order parameters
% ------------------------------------------------------------

psi4 = psi4 ./ Nb;
psi6 = psi6 ./ Nb;

% ------------------------------------------------------------
% Global, rotation-invariant order parameters
% ------------------------------------------------------------
Psi4 = sqrt(mean(abs(psi4).^2));   % square (C4)
Psi6 = sqrt(mean(abs(psi6).^2));   % hexagonal (C6)

end
