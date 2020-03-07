function lap_mat = laplace_mat(grid)

index = reshape((1:grid.Nx*grid.Ny),grid.Nx,grid.Ny);

lap_left_row = reshape(index(:,2:grid.Ny),grid.Nx*(grid.Ny-1),1);
lap_left_col = reshape(index(:,1:grid.Ny-1),grid.Nx*(grid.Ny-1),1);
lap_left_val = 1/grid.dy^2*ones(grid.Nx*(grid.Ny-1),1);

lap_right_row = reshape(index(:,1:grid.Ny-1),grid.Nx*(grid.Ny-1),1);
lap_right_col = reshape(index(:,2:grid.Ny),grid.Nx*(grid.Ny-1),1);
lap_right_val = 1/grid.dy^2*ones(grid.Nx*(grid.Ny-1),1);

lap_up_row = reshape(index(2:grid.Nx,:),(grid.Nx-1)*grid.Ny,1);
lap_up_col = reshape(index(1:grid.Nx-1,:),(grid.Nx-1)*grid.Ny,1);
lap_up_val = 1/grid.dx^2*ones((grid.Nx-1)*grid.Ny,1);

lap_down_row = reshape(index(1:grid.Nx-1,:),(grid.Nx-1)*grid.Ny,1);
lap_down_col = reshape(index(2:grid.Nx,:),(grid.Nx-1)*grid.Ny,1);
lap_down_val = 1/grid.dx^2*ones((grid.Nx-1)*grid.Ny,1);

lap_mid_row = reshape(index,grid.Nx*grid.Ny,1);
lap_mid_col = reshape(index,grid.Nx*grid.Ny,1);
lap_mid_val = (-2/grid.dx^2-2/grid.dy^2)*ones(grid.Nx*grid.Ny,1);

lap_mat = sparse(lap_left_row,lap_left_col,lap_left_val,grid.Nx*grid.Ny,grid.Nx*grid.Ny);
lap_mat = lap_mat + sparse(lap_right_row,lap_right_col,lap_right_val,grid.Nx*grid.Ny,grid.Nx*grid.Ny);
lap_mat = lap_mat + sparse(lap_up_row,lap_up_col,lap_up_val,grid.Nx*grid.Ny,grid.Nx*grid.Ny);
lap_mat = lap_mat + sparse(lap_down_row,lap_down_col,lap_down_val,grid.Nx*grid.Ny,grid.Nx*grid.Ny);
lap_mat = lap_mat + sparse(lap_mid_row,lap_mid_col,lap_mid_val,grid.Nx*grid.Ny,grid.Nx*grid.Ny);


