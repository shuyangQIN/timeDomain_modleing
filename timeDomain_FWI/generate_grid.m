function grid = generate_grid(Lx_comp, Ly_comp, dx, dy)
grid.Lx_comp = Lx_comp;
grid.Ly_comp = Ly_comp;
grid.dx = dx;
grid.dy = dy;

grid.Nx = round(2*Lx_comp/dx)+1;
grid.Ny = round(2*Ly_comp/dy)+1;
grid.x = (-Lx_comp:dx:Lx_comp)'*ones(1,grid.Ny);
grid.y = ones(grid.Nx,1)*(-Ly_comp:dy:Ly_comp);

grid.Nx_pml = 20;
grid.Ny_pml = 20;
grid.Lx = 1.000;
grid.Ly = 1.000;

grid.x_start = round((Lx_comp - grid.Lx)/dx) + 1;
grid.x_end = grid.Nx - round((Lx_comp - grid.Lx)/dx);
grid.y_start = round((Ly_comp - grid.Ly)/dy) + 1;
grid.y_end = grid.Ny - round((Ly_comp - grid.Ly)/dy);