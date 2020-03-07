function lap_vec = lap_operator(lap_input,grid)
if size(lap_input,2) == 1
    lap_input = reshape(lap_input,grid.Nx,grid.Ny);
end
lap_extend = zeros(grid.Nx+2,grid.Ny+2);
lap_extend(2:grid.Nx+1,2:grid.Ny+1) = lap_input;
lap_vec = (lap_extend(3:grid.Nx+2,2:grid.Ny+1) - 2*lap_extend(2:grid.Nx+1,2:grid.Ny+1) + lap_extend(1:grid.Nx,2:grid.Ny+1))/grid.dx^2 ...
    + (lap_extend(2:grid.Nx+1,3:grid.Ny+2) - 2*lap_extend(2:grid.Nx+1,2:grid.Ny+1) + lap_extend(2:grid.Nx+1,1:grid.Ny))/grid.dy^2;
lap_vec = reshape(lap_vec,grid.Nx*grid.Ny,1);