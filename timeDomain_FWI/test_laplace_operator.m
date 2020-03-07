clear
grid.x = (0:0.005:1)'*ones(1,201);
grid.y = ones(201,1)*(0:0.005:1);
grid.Nx = 201;
grid.Ny = 201;
grid.dx = 0.005;
grid.dy = 0.005;
U_ext = exp(grid.x).*sin(grid.y);
lap_mat = laplace_mat(grid);
lap_U = lap_mat*reshape(U_ext,201^2,1);
lap_U = reshape(lap_U,201,201);
imagesc(lap_U(2:200,2:200))