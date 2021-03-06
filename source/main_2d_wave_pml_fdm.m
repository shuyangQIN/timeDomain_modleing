%%  Main Script --- forward modeling wave propagation in 2D with perfectly matched layer
%   reference -- The perfectly matched layer for acoustic waves
%                in absorptive media, Author: Qing-Huo Liu & Jianping Tao
clear
clc
close all
addpath('source_func')

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
Lx_comp = 1300;
Ly_comp = 1300;
dx = 10;
dy = 10;
Nx = round(2*Lx_comp/dx)+1;
Ny = round(2*Ly_comp/dy)+1;
grid.x = (-Lx_comp:dx:Lx_comp)'*ones(1,Ny);
grid.y = ones(Nx,1)*(-Ly_comp:dy:Ly_comp);

Nx_pml = 20;
Ny_pml = 20;
Lx = 1000;
Ly = 1000;

x_start = round((Lx_comp - Lx)/dx) + 1;
x_end = Nx - round((Lx_comp - Lx)/dx);
y_start = round((Ly_comp - Ly)/dy) + 1;
y_end = Ny - round((Ly_comp - Ly)/dy);

%%  Colormap for visualization
d1 = [255 255-128 0]/(255*63);
d2 = [0 0 255]/(255*63);

c1 = zeros(64,3);
c1(1,:) = [0 128 255]/255;
c2 = zeros(64,3);
c2(1,:) = [255 255 0]/255;
for j = 2:64
    c1(j,:) = c1(1,:)+(j-1)*d1;
    c2(j,:) = c2(1,:)+(j-1)*d2;
end

c2 = flipud(c2);
cmap = [c1;c2(2:64,:)];

%%  Medium: background velocity + inhomogeneous part
medium.density = 1200;
medium.speed = 1000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
%medium.speed = medium.speed + 500 * (dist < 500);
%medium.speed = medium.speed+(200*sin(0.002*pi*grid.x)+100*cos(0.002*pi*grid.y)).*(abs(grid.x)<=Lx).*(abs(grid.y)<=Ly);

% smoothing the medium
window_type = 'Gaussian';
medium.speed = smooth_source(medium.speed,window_type);
figure
imagesc(medium.speed)
colormap(cmap)
title('True Velocity Model')
colorbar

% CFL number: 0.5
grid.dt = sqrt(0.25*dx^2/max(max(medium.speed))^2);
grid.Nt = ceil(sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.speed))));
% choose Nt to be bigger than 2 times the longest geodesic

%%  Source: either time dependent or initial time
% specs:
%   source_loc: point source locations on computational grid
%   source_mag: magnitude of point sources
%   source_grid: row/column number of sources in the matrix
%   source.time_dependent: time dependent 0 & independent 1
source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; 0 0]*1e3;
source_mag = [1;1;1;1;1]*1e6;
source_index = [1;2];
source_num = size(source_index,1);
source_grid = [round(source_loc(source_index,1)/dx)+(Nx+1)/2 round(source_loc(source_index,2)/dx)+(Ny+1)/2];
source.time_dependent = 0;
if source.time_dependent == 0 % time dependent
    % due to limited memory, only save source-time function
    source.firing_time = [0.2;0.5]; % point source firing time
    if size(source.firing_time,1) ~= source_num
        error('Wrong firing time input! Timing should be equal to %d\n', num_source);
    end
    source.peak_freq = 10;  % peak frequency of source wavelet
    a = pi^2*source.peak_freq^2;
    t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),2,1);
    amplitude = source_mag(source_index) .* exp(-a*(t-source.firing_time).^2); % magnitude?
    source.p0 = zeros(Nx,Ny);
elseif source.time_dependent == 1 % initial displacement
    source.p0 = zeros(Nx,Ny);
    for j = 1:source_num
        source.p0(source_grid(j,1),source_grid(j,2)) = source_mag(source_index(j));
    end
    window_type = 'Hanning';
    source.p0 = smooth_source(source.p0,window_type);
    figure
    imagesc(source.p0)
    colormap(cmap)
    title('Initial Displacement')
    colorbar
else
    error('Unsupported source function!')
end

%%  Sensor: location of receivers
sensor.mask = rect_sensor([x_start;y_start], [x_end;y_end], grid.x, grid.y);
sensor_data_fdm = zeros(size(sensor.mask,2),grid.Nt);

%%  Finite Difference Time Domain + Absorbing Boundary Conditions
%   Use first order systems to to marching forward in time

%%  PML absorption coefficients
grid.PML_size_x = 20;   % number of grids in PML region along one side of x-direction
grid.PML_size_y = 20;   % number of grids in PML region along one side of y-direction
const = 1e8;
grid.PML_length_x = grid.PML_size_x * dx;
grid.PML_length_y = grid.PML_size_y * dy;
PML_size_x = Lx_comp - grid.PML_length_x;
PML_size_y = Ly_comp - grid.PML_length_y;
% compute PML coefficients on staggered grid
grid.xPML_std_grid = (const/grid.PML_length_x)*((abs(grid.x)>PML_size_x).*(abs(grid.x)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_std_grid = (const/grid.PML_length_y)*((abs(grid.y)>PML_size_y).*(abs(grid.y)-PML_size_y).^4/grid.PML_length_y^4);
grid.xPML_std_grid = grid.xPML_std_grid(2:Nx-1,2:Ny-1);
grid.yPML_std_grid = grid.yPML_std_grid(2:Nx-1,2:Ny-1);

% compute PML coefficients on x-direction staggered grid
grid.stag_x_x = (grid.x(1:Nx-1,:) + grid.x(2:Nx,:))/2;
grid.stag_x_y = (grid.y(1:Nx-1,:) + grid.y(2:Nx,:))/2;
grid.xPML_stag_x = (const/grid.PML_length_x)*((abs(grid.stag_x_x)>PML_size_x).*(abs(grid.stag_x_x)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_stag_x = (const/grid.PML_length_y)*((abs(grid.stag_x_y)>PML_size_y).*(abs(grid.stag_x_y)-PML_size_y).^4/grid.PML_length_y^4);
clear grid.stag_x_x; clear grid.stag_x_y

% compute PML coefficients on y-direction staggered grid
grid.stag_y_x = (grid.x(:,1:Ny-1) + grid.x(:,2:Ny))/2;
grid.stag_y_y = (grid.y(:,1:Ny-1) + grid.y(:,2:Ny))/2;
grid.xPML_stag_y = (const/grid.PML_length_x)*((abs(grid.stag_y_x)>PML_size_x).*(abs(grid.stag_y_x)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_stag_y = (const/grid.PML_length_y)*((abs(grid.stag_y_y)>PML_size_y).*(abs(grid.stag_y_y)-PML_size_y).^4/grid.PML_length_y^4);
clear grid.stag_y_x; clear grid.stag_y_y

% compute PML coefficients on xy-direction staggered grid
grid.stag_xy_x = (grid.x(1:Nx-1,1:Ny-1) + grid.x(1:Nx-1,2:Ny) + grid.x(2:Nx,1:Ny-1) + grid.x(2:Nx,2:Ny))/4;
grid.stag_xy_y = (grid.y(1:Nx-1,1:Ny-1) + grid.y(1:Nx-1,2:Ny) + grid.y(2:Nx,1:Ny-1) + grid.y(2:Nx,2:Ny))/4;
grid.xPML_stag_xy = (const/grid.PML_length_x)*((abs(grid.stag_xy_x)>PML_size_x).*(abs(grid.stag_xy_x)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_stag_xy = (const/grid.PML_length_y)*((abs(grid.stag_xy_y)>PML_size_y).*(abs(grid.stag_xy_y)-PML_size_y).^4/grid.PML_length_y^4);
clear grid.stag_xy_x; clear grid.stag_xy_y

%% Initialize Pressure field
%   interpolate pressure field onto staggered grid
Up_interp = (source.p0(1:Nx-1,1:Ny-1) + source.p0(1:Nx-1,2:Ny) + source.p0(2:Nx,1:Ny-1) + source.p0(2:Nx,2:Ny))/4;
%   interpolate pressure field onto staggered grid
speed_interp = (medium.speed(1:Nx-1,1:Ny-1) + medium.speed(1:Nx-1,2:Ny) + medium.speed(2:Nx,1:Ny-1) + medium.speed(2:Nx,2:Ny))/4;
Up_interp_extend = zeros(Nx+1,Ny+1);
Up_interp_extend(2:Nx,2:Ny) = Up_interp;
sensor_data_fdm(:,1) = record_data(source.p0, [x_start;y_start], [x_end;y_end]);

%% Initialize Particle Velocity field
%   assume the particle velocity is small enough to be ignored
%   initialize vx and vy at t=0 to be 0

Vx = -(Up_interp_extend(2:Nx+1,2:Ny) - Up_interp_extend(1:Nx,2:Ny))*grid.dt/(2*dx*medium.density);
Vy = -(Up_interp_extend(2:Nx,2:Ny+1) - Up_interp_extend(2:Nx,1:Ny))*grid.dt/(2*dy*medium.density);

Up_x = Up_interp/2;
Up_y = Up_interp/2;

%   assembling matrix for updating pressure
B1_x = 1/grid.dt + grid.xPML_stag_xy/2;
B2_x = 1/grid.dt - grid.xPML_stag_xy/2;
B1_y = 1/grid.dt + grid.yPML_stag_xy/2;
B2_y = 1/grid.dt - grid.yPML_stag_xy/2;
%   assembling matrix for updating particle velocity
C1_x = medium.density/grid.dt + medium.density*grid.xPML_stag_y/2;
C2_x = medium.density/grid.dt - medium.density*grid.xPML_stag_y/2;
C1_y = medium.density/grid.dt + medium.density*grid.yPML_stag_x/2;
C2_y = medium.density/grid.dt - medium.density*grid.yPML_stag_x/2;

%% Marching forward in time
for j = 2:grid.Nt
    % compute pressure field based on particle velocity field
    Up_x = (-medium.density*speed_interp.^2.*(Vx(2:Nx,:)-Vx(1:Nx-1,:))/dx + B2_x.*Up_x)./B1_x;
    Up_y = (-medium.density*speed_interp.^2.*(Vy(:,2:Ny)-Vy(:,1:Ny-1))/dy + B2_y.*Up_y)./B1_y;
    Up_interp = Up_x + Up_y;
    % extend the wavefield to compute velocity field
    Up_interp_extend(2:Nx,2:Ny) = Up_interp;
    % recording wave field at sensor locations
    Up = (Up_interp_extend(1:Nx,1:Ny) + Up_interp_extend(1:Nx,2:Ny+1) + Up_interp_extend(2:Nx+1,1:Ny) + Up_interp_extend(2:Nx+1,2:Ny+1))/4;
    if source.time_dependent == 0
        p0 = zeros(size(source.p0));
        for k = 1:source_num
            p0(source_grid(k,1),source_grid(k,2)) = amplitude(k,j);
        end
        %p0 = smooth_source(p0,window_type);
        Up = Up + grid.dt*p0;
        Up_interp = (Up(1:Nx-1,1:Ny-1) + Up(1:Nx-1,2:Ny) + Up(2:Nx,1:Ny-1) + Up(2:Nx,2:Ny))/4;
        Up_interp_extend(2:Nx,2:Ny) = Up_interp;
    end
    sensor_data_fdm(:,j) = record_data(Up, [x_start;y_start], [x_end;y_end]);
    % compute particle velocity field based on pressure field
    Vx = (-(Up_interp_extend(2:Nx+1,2:Ny) - Up_interp_extend(1:Nx,2:Ny))/dx + C2_x.*Vx)./C1_x;
    Vy = (-(Up_interp_extend(2:Nx,2:Ny+1) - Up_interp_extend(2:Nx,1:Ny))/dy + C2_y.*Vy)./C1_y;
    if mod(j,5) == 0
        imagesc(Up_interp_extend)
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
plot(sensor_data_fdm(1,:))
save('data/sensor_data_fdm')