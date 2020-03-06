dx = 0.01;
Nx = 201;
grid.x = (-1:dx:1)'*ones(1,201);
grid.y = ones(201,1)*(-1:dx:1);
U_ext = exp(grid.x).*sin(grid.y);
V_ext = zeros(201,201);

U1_ext = -exp(grid.x).*sin(grid.y);
U2_ext = exp(grid.x).*cos(grid.y);
U3_ext = exp(grid.x).*sin(grid.y);
U4_ext = -exp(grid.x).*cos(grid.y);

Dirichlet = [U_ext(1,1:200)';U_ext(1:200,201);U_ext(201,201:-1:2)';U_ext(201:-1:2,1)];
Neumann = [U1_ext(1,:)';U2_ext(:,201);U3_ext(201,201:-1:1)';U4_ext(201:-1:1,1)];

[U,V] = biharmonic_solver(Dirichlet,Neumann,dx,Nx);

figure
imagesc(U-U_ext(2:200,2:200))
colormap jet
colorbar

figure
imagesc(V-V_ext)
colormap jet
colorbar