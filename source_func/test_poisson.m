%% test poisson equation solver
N = [20;40;80;160];

for j = 1:size(N,1)
    Nx = N(j)+1;
    Ny = N(j)+1;
    dx = 1/N(j);
    dy = 1/N(j);

    grid.x = (0:dx:1)'*ones(1,Ny);
    grid.y = ones(Nx,1)*(0:dy:1);

    f = @(x,y) exp(x.*y);
    U_ext = f(grid.x,grid.y);
%     figure
%     imagesc(U_ext)
%     colormap jet
%     colorbar

    measurement = [U_ext(1,1:Ny-1)';U_ext(1:Nx-1,Ny);U_ext(Nx,Ny:-1:2)';U_ext(Nx:-1:2,1)];
    H = -(grid.x.^2 + grid.y.^2).*exp(grid.x.*grid.y);
    F = H(2:Nx-1,2:Ny-1);
    U = poisson_solver(measurement,Nx,Ny,dx,dy,F);
    
    figure
    imagesc(abs(U-U_ext)./abs(U_ext))
    colormap jet
    colorbar
end
% figure
% imagesc(U)
% colormap jet
% colorbar
% 
% figure
% imagesc(abs(U-U_ext)./abs(U_ext))
% colormap jet
% colorbar

%% test the convergence order