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