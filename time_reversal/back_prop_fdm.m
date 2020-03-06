function U = back_prop_fdm(U0, Ut, sensor_data, speed, dt, dx, dy)
cmap = customize_colormap(2);

% Ut = smooth_source(Ut,'Gaussian');
[Nx,Ny] = size(U0);

% U_old = smooth_source(U0,'Gaussian'); % pressure field at final time
U_old = U0;
U_now = U0 - dt*Ut;
U_now = smooth_source(U_now,'Gaussian');

% Apply boundary condition
Nt = size(sensor_data,2);
U_now(1,1:Ny-1) = sensor_data(1:Ny-1,Nt-1)';
U_now(1:Nx-1,Ny) = sensor_data(Ny:Nx+Ny-2,Nt-1);
U_now(Nx,Ny:-1:2) = sensor_data(Nx+Ny-1:Nx+2*Ny-3,Nt-1)';
U_now(Nx:-1:2,1) = sensor_data(Nx+2*Ny-2:2*Nx+2*Ny-4,Nt-1);
%U_now = smooth_source(U_now,'Gaussian');

for iter = Nt-2:-1:1
    U = zeros(Nx,Ny);
    U(1,1:Ny-1) = sensor_data(1:Ny-1,iter)';
    U(1:Nx-1,Ny) = sensor_data(Ny:Nx+Ny-2,iter);
    U(Nx,Ny:-1:2) = sensor_data(Nx+Ny-1:Nx+2*Ny-3,iter)';
    U(Nx:-1:2,1) = sensor_data(Nx+2*Ny-2:2*Nx+2*Ny-4,iter);
    U(2:Nx-1,2:Ny-1) = -U_old(2:Nx-1,2:Ny-1)+2*U_now(2:Nx-1,2:Ny-1)+(dt^2/dx^2)*speed(2:Nx-1,2:Ny-1).^2.*...
        (U_now(1:Nx-2,2:Ny-1)+U_now(3:Nx,2:Ny-1)-2*U_now(2:Nx-1,2:Ny-1))+...
        (dt^2/dy^2)*speed(2:Nx-1,2:Ny-1).^2.*(U_now(2:Nx-1,1:Ny-2)+U_now(2:Nx-1,3:Ny)-2*U_now(2:Nx-1,2:Ny-1));
%     if iter > Nt-100
%         U = smooth_source(U,'Gaussian');
%     end
    U_old = U_now;
    U_now = U;
    if mod(iter,40) == 0
        imagesc(U)
        title(sprintf('time =%e s',(iter-1)*dt));
        caxis([-1,1])
        colormap(cmap)
        colorbar
        axis equal
        drawnow
    end
end






