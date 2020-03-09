%%  Main Script --- forward modeling wave propagation in 2D with perfectly matched layer
%   reference -- The perfectly matched layer for acoustic waves
%                in absorptive media, Author: Qing-Huo Liu & Jianping Tao
clear
clc
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
Lx_comp = 1.300;
Ly_comp = 1.300;
dx = 0.01;
dy = 0.01;
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
medium.speed = medium.speed + 0.500 * (dist < 0.500);
%medium.speed = medium.speed+(0.200*sin(2*pi*grid.x)+0.100*cos(2*pi*grid.y)).*(abs(grid.x)<=0.9*Lx).*(abs(grid.y)<=0.9*Ly);

% smoothing the medium
window_type = 'Gaussian';
medium.speed = smooth_source(medium.speed,window_type);
figure
imagesc(medium.speed)
colormap(cmap)
title('True Velocity Model')
colorbar

% CFL number: 0.5
grid.dt = sqrt(0.04*dx^2/max(max(medium.speed))^2);
grid.Nt = ceil(1.2*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.speed))));
grid.Nt = round(grid.Nt*1.5);

% choose Nt to be bigger than 2 times the longest geodesic

%%  Source: either time dependent or initial time
phantom_image = phantom('Modified Shepp-Logan',x_end-x_start+1);
source.p0 = zeros(Nx,Ny);
source.p0(x_start:x_end,y_start:y_end) = phantom_image;
%imshow(source.p0)
window_type = 'Gaussian';
source.p0 = smooth_source(source.p0,window_type);
% imagesc(source.p0)
% colormap(cmap)
% colorbar

%%  Sensor: location of receivers
sensor.mask = rect_sensor([x_start;y_start], [x_end;y_end], grid.x, grid.y);
sensor_data_fdm11 = zeros(size(sensor.mask,2),grid.Nt);

%%  Finite Difference Time Domain + Absorbing Boundary Conditions
%   Use first order systems to to marching forward in time

%%  PML absorption coefficients
grid.PML_size_x = 20;   % number of grids in PML region along one side of x-direction
grid.PML_size_y = 20;   % number of grids in PML region along one side of y-direction
const = 1e2;
grid.PML_length_x = grid.PML_size_x * dx;
grid.PML_length_y = grid.PML_size_y * dy;
PML_size_x = Lx_comp - grid.PML_length_x;
PML_size_y = Ly_comp - grid.PML_length_y;

% compute PML coefficients on standard computational grid
grid.xPML_std_grid = grid.x(2:Nx-1,2:Ny-1);
grid.yPML_std_grid = grid.y(2:Nx-1,2:Ny-1);
grid.xPML_std_x = (const/grid.PML_length_x)*((abs(grid.xPML_std_grid)>PML_size_x).*(abs(grid.xPML_std_grid)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_std_y = (const/grid.PML_length_y)*((abs(grid.yPML_std_grid)>PML_size_y).*(abs(grid.yPML_std_grid)-PML_size_y).^4/grid.PML_length_y^4);

% compute PML coefficients on staggered grid in x-direction
grid.xPML_stag_grid = (grid.x(2:Nx,2:Ny-1) + grid.x(1:Nx-1,2:Ny-1))/2;
grid.xPML_stag_x = (const/grid.PML_length_x)*((abs(grid.xPML_stag_grid)>PML_size_x).*(abs(grid.xPML_stag_grid)-PML_size_x).^4/grid.PML_length_x^4);

% compute PML coefficients on staggered grid in y-direction
grid.yPML_stag_grid = (grid.y(2:Nx-1,2:Ny) + grid.y(2:Nx-1,1:Ny-1))/2;
grid.yPML_stag_y = (const/grid.PML_length_y)*((abs(grid.yPML_stag_grid)>PML_size_y).*(abs(grid.yPML_stag_grid)-PML_size_y).^4/grid.PML_length_y^4);

%% Initialize Pressure field
%   pressure field initialized in the interior
Up = source.p0(2:Nx-1,2:Nx-1);
%   wave speed in the interior
speed = medium.speed(2:Nx-1,2:Ny-1);
Up_extend = source.p0;
sensor_data_fdm11(:,1) = record_data(Up_extend, [x_start;y_start], [x_end;y_end]);

%% Initialize Particle Velocity field
%   assume the particle velocity is small enough to be ignored
%   initialize vx and vy at t=0 to be 0

Vx = -(Up_extend(2:Nx,2:Ny-1) - Up_extend(1:Nx-1,2:Ny-1))*grid.dt/(2*dx*medium.density);
Vy = -(Up_extend(2:Nx-1,2:Ny) - Up_extend(2:Nx-1,1:Ny-1))*grid.dt/(2*dy*medium.density);

% initialize pressure field components
Up_x = Up/2;
Up_y = Up/2;

%   assembling matrix for updating pressure
B1_x = 1/grid.dt + grid.xPML_std_x/2;
B2_x = 1/grid.dt - grid.xPML_std_x/2;
B1_y = 1/grid.dt + grid.yPML_std_y/2;
B2_y = 1/grid.dt - grid.yPML_std_y/2;
%   assembling matrix for updating particle velocity
C1_x = medium.density/grid.dt + medium.density*grid.xPML_stag_x/2;
C2_x = medium.density/grid.dt - medium.density*grid.xPML_stag_x/2;
C1_y = medium.density/grid.dt + medium.density*grid.yPML_stag_y/2;
C2_y = medium.density/grid.dt - medium.density*grid.yPML_stag_y/2;

%% Marching forward in time
for j = 2:grid.Nt
    % compute pressure field based on particle velocity field
    Up_x = (-medium.density*speed.^2.*(Vx(2:Nx-1,:)-Vx(1:Nx-2,:))/dx + B2_x.*Up_x)./B1_x;
    Up_y = (-medium.density*speed.^2.*(Vy(:,2:Ny-1)-Vy(:,1:Ny-2))/dy + B2_y.*Up_y)./B1_y;
    Up = Up_x + Up_y;
    % extend the wavefield to compute velocity field
    Up_extend(2:Nx-1,2:Ny-1) = Up;
    p0 = zeros(size(source.p0));
%     for k = 1:source_num
%         p0(source_grid(k,1),source_grid(k,2)) = amplitude(k,j);
%     end
%     if j <= 1
%         p0 = smooth_source(p0,window_type);
%     end
    Up_extend = Up_extend + p0;
    % recording wave field at sensor locations
    sensor_data_fdm11(:,j) = record_data(Up_extend, [x_start;y_start], [x_end;y_end]);
    % compute particle velocity field based on pressure field
    Vx = (-(Up_extend(2:Nx,2:Ny-1) - Up_extend(1:Nx-1,2:Ny-1))/dx + C2_x.*Vx)./C1_x;
    Vy = (-(Up_extend(2:Nx-1,2:Ny) - Up_extend(2:Nx-1,1:Ny-1))/dy + C2_y.*Vy)./C1_y;
    if mod(j,40) == 0
        imagesc(Up_extend)
        %imagesc(p0)
        title(sprintf('time =%e s',(j-1)*grid.dt));
        caxis([-1,1])
        colormap(cmap)
        colorbar
        %axis equal
        drawnow
    end
end

%%
plot(sensor_data_fdm11(1,:))
save('data/sensor_data_fdm11')

%%  Down-sampling recorded data
sampling_rate = 1;
sensor_data_back_prop = sensor_data_fdm11(1:sampling_rate:end,:);
x_temp = x_end - x_start + 1; y_temp = y_end - y_start + 1;
neumann1 = [Up_extend(x_start-1,y_start:y_end)';Up_extend(x_start:x_end,y_end+1);...
    Up_extend(x_end+1,y_end:-1:y_start)';Up_extend(x_end:-1:x_start,y_start-1)];
neumann2 = [Up_extend(x_start+1,y_start:y_end)';Up_extend(x_start:x_end,y_end-1);...
    Up_extend(x_end-1,y_end:-1:y_start)';Up_extend(x_end:-1:x_start,y_start+1)];
neumann_data_back_prop = (neumann1 - neumann2)/(2*dx);

%% time reversal -- 0 initial data for back propagation
lengthx = x_end-x_start+1; lengthy = y_end-y_start+1;
U0 = zeros(lengthx,lengthy);
Ut = zeros(lengthx,lengthy);
speed_interior = medium.speed(x_start:x_end,y_start:y_end);
p0_recon = back_prop_fdm(U0, Ut, sensor_data_back_prop, speed_interior, grid.dt, dx, dy);
err0 = p0_recon - source.p0(x_start:x_end,y_start:y_end);
fprintf('Frobenius norm via 0 initial data: %f\n',norm(err0,'fro'));

% imagesc(p0_recon - phantom_image)
% norm(p0_recon - phantom_image,'fro')
% figure
% imagesc(p0_recon - source.p0(x_start:x_end,y_start:y_end))
% norm(p0_recon - source.p0(x_start:x_end,y_start:y_end),'fro')
% colormap(cmap)
% colorbar

%% time reversal -- harmonic extension
F = zeros(lengthx-2,lengthy-2);
measurement = sensor_data_back_prop(:,grid.Nt);
U0 = poisson_solver(measurement,lengthx,lengthy,dx,dy,F);
Ut = zeros(lengthx,lengthy);
speed_interior = medium.speed(x_start:x_end,y_start:y_end);
p1_recon = back_prop_fdm(U0, Ut, sensor_data_back_prop, speed_interior, grid.dt, dx, dy);
err1 = p1_recon - source.p0(x_start:x_end,y_start:y_end);
fprintf('Frobenius norm via harmonic extension: %f\n',norm(err1,'fro'))

% imagesc(p1_recon - phantom_image)
% norm(p1_recon - phantom_image,'fro')
% figure
% imagesc(p1_recon - source.p0(x_start:x_end,y_start:y_end))
% norm(p1_recon - source.p0(x_start:x_end,y_start:y_end),'fro')
% colormap(cmap)
% colorbar

%% time reversal -- biharmonic extension
F = zeros(lengthx-2,lengthy-2);
measurement = sensor_data_back_prop(:,grid.Nt);
U0 = poisson_solver(measurement,lengthx,lengthy,dx,dy,F);
dirichlet_data = sensor_data_back_prop(:,grid.Nt);
neumann_data = neumann_data_back_prop;
[U,V] = biharmonic_solver(dirichlet_data,neumann_data,dx,x_end-x_start+1);
Ut = -(V.*(V<0));
speed_interior = medium.speed(x_start:x_end,y_start:y_end);
p2_recon = back_prop_fdm(U0, Ut, sensor_data_back_prop, speed_interior, grid.dt, dx, dy);
err2 = p2_recon - source.p0(x_start:x_end,y_start:y_end);
fprintf('Frobenius norm via biharmonic extension: %f\n',norm(err2,'fro'))

%% visualization
figure
imagesc(err0)
colormap(cmap)
colorbar

figure
imagesc(err1)
colormap(cmap)
colorbar

figure
imagesc(err2)
colormap(cmap)
colorbar