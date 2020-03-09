function V = neumann_adjoint_wave(neumann_data, speed, Nx, Ny, Nt, dx, dy, dt)
if size(neumann_data,1) ~= 2*(Nx+Ny-2)
    error('Data size not compatible!')
end
cmap = customize_colormap(2);
V = zeros(Nx,Ny,Nt);
% at time step t = T
V(:,:,Nt) = zeros(Nx,Ny);
Vt_extend  = apply_neumann_condition(V(:,:,Nt),neumann_data(:,Nt),dx,dy);
% at time step t = T - dt
V(:,:,Nt-1) = V(:,:,Nt) + (dt^2/2)*speed.^2.*...
    ((Vt_extend(3:Nx+2,2:Ny+1)-2*Vt_extend(2:Nx+1,2:Ny+1)+Vt_extend(1:Nx,2:Ny+1))/dx^2 ...
    + (Vt_extend(2:Nx+1,3:Ny+2)-2*Vt_extend(2:Nx+1,2:Ny+1)+Vt_extend(2:Nx+1,1:Ny))/dy^2);
Vt_extend  = apply_neumann_condition(V(:,:,Nt-1),neumann_data(:,Nt-1),dx,dy);
for j = Nt-2:-1:1
    V(:,:,j) = 2*V(:,:,j+1) - V(:,:,j+2) + dt^2*speed.^2.*...
        ((Vt_extend(3:Nx+2,2:Ny+1)-2*Vt_extend(2:Nx+1,2:Ny+1)+Vt_extend(1:Nx,2:Ny+1))/dx^2 ...
        + (Vt_extend(2:Nx+1,3:Ny+2)-2*Vt_extend(2:Nx+1,2:Ny+1)+Vt_extend(2:Nx+1,1:Ny))/dy^2);
    Vt_extend  = apply_neumann_condition(V(:,:,j),neumann_data(:,j),dx,dy);
    if mod(j,20) == 0
        imagesc(V(:,:,j))
        %imagesc(p0)
        title(sprintf('time =%e s',(j-1)*dt));
        caxis([-1e-3,1e-3])
        colormap(cmap)
        colorbar
        %axis equal
        drawnow
    end
end
end

function Vt_extend = apply_neumann_condition(Vt,data,dx,dy)
[Nx,Ny] = size(Vt);
Vt_extend = zeros(Nx+2,Ny+2);
Vt_extend(2:Nx+1,2:Ny+1) = Vt;
data_extend = [data(1:Ny);data(Ny:Nx+Ny-1);data(Nx+Ny-1:Nx+2*Ny-2);data(Nx+2*Ny-2:2*Nx+2*Ny-4);data(1)];
Vt_extend(1,2:Ny+1) = Vt_extend(3,2:Ny+1) + 2*dx*data_extend(1:Ny)';
Vt_extend(2:Nx+1,Ny+2) = Vt_extend(2:Nx+1,Ny) + 2*dy*data_extend(Ny+1:Nx+Ny);
Vt_extend(Nx+2,Ny+1:-1:2) = Vt_extend(Nx,Ny+1:-1:2) + 2*dx*data_extend(Nx+Ny+1:Nx+2*Ny)';
Vt_extend(Nx+1:-1:2,1) = Vt_extend(Nx+1:-1:2,3) + 2*dy*data_extend(Nx+2*Ny+1:2*Nx+2*Ny);
end

