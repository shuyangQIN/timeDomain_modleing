%%  reproducing dissipative time reversal algorithm
%   reference -- A dissipative time reversal technique for photoacoustic
%       tomography in a cavity
%                Author: Linh V. Nguyen and Leonid A. Kunyansky
clear
%clc
close all
addpath('../source_func')

%%  Computational Domain
% specs: computational domain centered at (0,0)
%     Lx_comp & Ly_comp: one-sided length of comp domain
%     dx & dy: spatial spacings
%     Nx & Ny: total spatial grids
%     grid: struct variable, location of grids
%     Nx_pml & Ny_pml: number of PML grids on each side
%     Lx & Ly: interior domain length
%     x_start & y_start: the element of left upper most grid in interior domain
%     x_end & y_end: the element of right lower most grid in interior domain
Lx_comp = 1.000;
Ly_comp = 1.000;
dx = 0.010;
dy = 0.010;
Nx = round(2*Lx_comp/dx)+1;
Ny = round(2*Ly_comp/dy)+1;
grid.x = (-Lx_comp:dx:Lx_comp)'*ones(1,Ny);
grid.y = ones(Nx,1)*(-Ly_comp:dy:Ly_comp);

Nx_pml = 20;
Ny_pml = 20;
Lx = 1.000;
Ly = 1.000;

x_start = round((Lx_comp - Lx)/dx) + 1;
x_end = Nx - round((Lx_comp - Lx)/dx);
y_start = round((Ly_comp - Ly)/dy) + 1;
y_end = Ny - round((Ly_comp - Ly)/dy);

%%  Colormap for visualization
cmap = customize_colormap(2);

%%  Medium: background velocity + inhomogeneous part
medium.density = 1;
medium.speed = 1.000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
%medium.speed = medium.speed + 0.500 * (dist < 0.500);
medium.speed = medium.speed+(0.200*sin(2*pi*grid.x)+0.100*cos(2*pi*grid.y)).*(abs(grid.x)<=Lx).*(abs(grid.y)<=Ly);

% smoothing the medium
window_type = 'Gaussian';
medium.speed = smooth_source(medium.speed,window_type);
figure
imagesc(medium.speed)
colormap(cmap)
title('True Velocity Model')
colorbar

% CFL number: 0.1
grid.dt = sqrt(0.04*dx^2/max(max(medium.speed))^2);
grid.Nt = ceil(1.2*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.speed))));
grid.Nt = round(grid.Nt*2);

% choose Nt to be bigger than 2 times the longest geodesic

%%  Source: either time dependent or initial time
phantom_image = phantom('Modified Shepp-Logan',x_end-x_start+1);
source.p0 = zeros(Nx,Ny);
source.p0(x_start:x_end,y_start:y_end) = phantom_image;
%imshow(source.p0)
window_type = 'Gaussian';
% source.p0 = smooth_source(source.p0,window_type);
% imagesc(source.p0)
% colormap(cmap)
% colorbar

%%  Sensor: location of receivers
sensor.mask = rect_sensor([x_start;y_start], [x_end;y_end], grid.x, grid.y);
sensor_data_fdm_neumann = zeros(size(sensor.mask,2),grid.Nt);

%%  Initialize Pressure field -- extension as neumann boundary condition
U_old = zeros(Nx+2,Ny+2);
U_old(2:Nx+1,2:Ny+1) = source.p0;
U_old(1,2:Nx+1) = U_old(3,2:Nx+1);
U_old(Nx+2,2:Nx+1) = U_old(Nx,2:Nx+1);
U_old(2:Ny+1,1) = U_old(2:Ny+1,3);
U_old(2:Ny+1,Ny+2) = U_old(2:Ny+1,Ny);
sensor_data_fdm_neumann(:,1) = record_data(U_old, [x_start+1;y_start+1], [x_end+1;y_end+1]);

U_now = zeros(Nx+2,Ny+2);
U_now(2:Nx+1,2:Ny+1) = U_old(2:Nx+1,2:Ny+1) + (grid.dt^2/2)*medium.speed.^2.*...
    ((U_old(3:Nx+2,2:Ny+1)-2*U_old(2:Nx+1,2:Ny+1)+U_old(1:Nx,2:Ny+1))/dx^2 + ...
    (U_old(2:Nx+1,3:Ny+2)-2*U_old(2:Nx+1,2:Ny+1)+U_old(2:Nx+1,1:Ny))/dy^2);
U_now(1,2:Nx+1) = U_now(3,2:Nx+1);
U_now(Nx+2,2:Nx+1) = U_now(Nx,2:Nx+1);
U_now(2:Ny+1,1) = U_now(2:Ny+1,3);
U_now(2:Ny+1,Ny+2) = U_now(2:Ny+1,Ny);
sensor_data_fdm_neumann(:,2) = record_data(U_now, [x_start+1;y_start+1], [x_end+1;y_end+1]);

U = zeros(Nx+2,Ny+2);

%% Forward
for j = 3:grid.Nt
    U(2:Nx+1,2:Ny+1) = 2*U_now(2:Nx+1,2:Ny+1) - U_old(2:Nx+1,2:Ny+1) + (grid.dt^2)*medium.speed.^2.*...
        ((U_now(3:Nx+2,2:Ny+1)-2*U_now(2:Nx+1,2:Ny+1)+U_now(1:Nx,2:Ny+1))/dx^2 + ...
        (U_now(2:Nx+1,3:Ny+2)-2*U_now(2:Nx+1,2:Ny+1)+U_now(2:Nx+1,1:Ny))/dy^2);
    U(1,2:Nx+1) = U(3,2:Nx+1);
    U(Nx+2,2:Nx+1) = U(Nx,2:Nx+1);
    U(2:Ny+1,1) = U(2:Ny+1,3);
    U(2:Ny+1,Ny+2) = U(2:Ny+1,Ny);
    sensor_data_fdm_neumann(:,j) = record_data(U, [x_start+1;y_start+1], [x_end+1;y_end+1]);
    U_old = U_now;
    U_now = U;
    if mod(j,20) == 0
        imagesc(U(2:Nx+1,2:Ny+1))
        % imagesc(p0)
        title(sprintf('time =%e s',(j-1)*grid.dt));
        caxis([-1,1])
        colormap(cmap)
        colorbar
        % axis equal
        drawnow
    end
end

%% time reversal -- idea: compute wavefield in interior and boundary using wave equation
%   then update wavefield values in ghost layer, using lambda == 1

% at t = T, u_t = 0 and u = 0
U_old = zeros(Nx+2,Ny+2);

% Harmonic Extension
F = zeros(Nx-2,Ny-2);
measurement = sensor_data_fdm_neumann(:,grid.Nt);
U0 = poisson_solver(measurement,Nx,Ny,dx,dy,F);
U_old(2:Nx+1,2:Ny+1) = U0;

% second order approximation of gt
gt = (1.5*sensor_data_fdm_neumann(:,grid.Nt) - 2*sensor_data_fdm_neumann(:,grid.Nt-1)...
    + 0.5*sensor_data_fdm_neumann(:,grid.Nt-2))/grid.dt;
gt = [gt(1:Ny);gt(Ny:Nx+Ny-1);gt(Nx+Ny-1:Nx+2*Ny-2);gt(Nx+2*Ny-2:2*Nx+2*Ny-4);gt(1)];
% extend to the boundary
U_old(1,2:Ny+1) = U_old(3,2:Ny+1) -2*dx*gt(1:Ny)';
U_old(2:Nx+1,Ny+2) = U_old(2:Nx+1,Ny) - 2*dy*gt(Ny+1:Nx+Ny);
U_old(Nx+2,Ny+1:-1:2) = U_old(Nx,Ny+1:-1:2) - 2*dx*gt(Nx+Ny+1:Nx+2*Ny)';
U_old(Nx+1:-1:2,1) = U_old(Nx+1:-1:2,3) - 2*dy*gt(Nx+2*Ny+1:2*Nx+2*Ny);

% at t = T-dt
U_now = zeros(Nx+2,Ny+2);
U_now(2:Nx+1,2:Ny+1) = U_old(2:Nx+1,2:Ny+1) + (grid.dt^2/2)*medium.speed.^2.*...
    ((U_old(3:Nx+2,2:Ny+1)-2*U_old(2:Nx+1,2:Ny+1)+U_old(1:Nx,2:Ny+1))/dx^2 + ...
    (U_old(2:Nx+1,3:Ny+2)-2*U_old(2:Nx+1,2:Ny+1)+U_old(2:Nx+1,1:Ny))/dy^2);
% second order approximation of gt
gt = (sensor_data_fdm_neumann(:,grid.Nt) - sensor_data_fdm_neumann(:,grid.Nt-2))/(2*grid.dt);
gt = [gt(1:Ny);gt(Ny:Nx+Ny-1);gt(Nx+Ny-1:Nx+2*Ny-2);gt(Nx+2*Ny-2:2*Nx+2*Ny-4);gt(1)];
% extend to the boundary
U_now(1,2:Ny+1) = U_now(3,2:Ny+1) -2*dx*gt(1:Ny)';
U_now(2:Nx+1,Ny+2) = U_now(2:Nx+1,Ny) - 2*dy*gt(Ny+1:Nx+Ny);
U_now(Nx+2,Ny+1:-1:2) = U_now(Nx,Ny+1:-1:2) - 2*dx*gt(Nx+Ny+1:Nx+2*Ny)';
U_now(Nx+1:-1:2,1) = U_now(Nx+1:-1:2,3) - 2*dy*gt(Nx+2*Ny+1:2*Nx+2*Ny);

U = zeros(Nx+2,Ny+2);
Ut = zeros(2*Nx+2*Ny,1);

%% backward update in time
for j = grid.Nt-2:-1:1
    % update interior and boundary points
    U(2:Nx+1,2:Ny+1) = 2*U_now(2:Nx+1,2:Ny+1) - U_old(2:Nx+1,2:Ny+1) + (grid.dt^2)*medium.speed.^2.*...
        ((U_now(3:Nx+2,2:Ny+1)-2*U_now(2:Nx+1,2:Ny+1)+U_now(1:Nx,2:Ny+1))/dx^2 + ...
        (U_now(2:Nx+1,3:Ny+2)-2*U_now(2:Nx+1,2:Ny+1)+U_now(2:Nx+1,1:Ny))/dy^2);
    % second order approximation of gt
    gt = (-1.5*sensor_data_fdm_neumann(:,j) + 2*sensor_data_fdm_neumann(:,j+1)...
        - 0.5*sensor_data_fdm_neumann(:,j+2))/grid.dt;
    gt = [gt(1:Ny);gt(Ny:Nx+Ny-1);gt(Nx+Ny-1:Nx+2*Ny-2);gt(Nx+2*Ny-2:2*Nx+2*Ny-4);gt(1)];
    Ut(1:Ny) = (-1.5*U(2,2:Ny+1) + 2*U_now(2,2:Ny+1) - 0.5*U_old(2,2:Ny+1))'/grid.dt;
    Ut(Ny+1:Nx+Ny) = (-1.5*U(2:Nx+1,Ny+1) + 2*U_now(2:Nx+1,Ny+1) - 0.5*U_old(2:Nx+1,Ny+1))/grid.dt;
    Ut(Nx+Ny+1:Nx+2*Ny) = (-1.5*U(Nx+1,Ny+1:-1:2) + 2*U_now(Nx+1,Ny+1:-1:2) - 0.5*U_old(Nx+1,Ny+1:-1:2))'/grid.dt;
    Ut(Nx+2*Ny+1:2*Nx+2*Ny) = (-1.5*U(Nx+1:-1:2,2) + 2*U_now(Nx+1:-1:2,2) - 0.5*U_old(Nx+1:-1:2,2))/grid.dt;
    U(1,2:Ny+1) = U(3,2:Ny+1) + 2*dx*(Ut(1:Ny)'-gt(1:Ny)');
    U(2:Nx+1,Ny+2) = U(2:Nx+1,Ny) + 2*dy*(Ut(Ny+1:Nx+Ny)-gt(Ny+1:Nx+Ny));
    U(Nx+2,Ny+1:-1:2) = U(Nx,Ny+1:-1:2) + 2*dx*(Ut(Nx+Ny+1:Nx+2*Ny)'-gt(Nx+Ny+1:Nx+2*Ny)');
    U(Nx+1:-1:2,1) = U(Nx+1:-1:2,3) + 2*dy*(Ut(Nx+2*Ny+1:2*Nx+2*Ny)-gt(Nx+2*Ny+1:2*Nx+2*Ny));
    if mod(j,20) == 0
        imagesc(U(2:Nx+1,2:Ny+1))
        % imagesc(p0)
        title(sprintf('time =%e s',(j-1)*grid.dt));
        caxis([-1,1])
        colormap(cmap)
        colorbar
        % axis equal
        drawnow
    end
    U_old = U_now;
    U_now = U;
end

%% visualization
p0_recon = U(2:Nx+1,2:Ny+1);

shift = 0;
%shift = mean(mean(source.p0-p0_recon));

figure
imagesc(source.p0)
colormap(cmap)
caxis([-1,1])
colorbar
title('True image')

figure
imagesc(p0_recon+shift)
colormap(cmap)
caxis([-1,1])
colorbar
title('Reconstructed image')

figure
imagesc(p0_recon - source.p0 + shift)
colormap(cmap)
%caxis([-1,1])
colorbar
title('Pixel-wise error')

fprintf('Frobenius norm of reconstruction: %f\n',norm(p0_recon - source.p0 + shift,'fro'));



