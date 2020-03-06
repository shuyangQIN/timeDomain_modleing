%%  Main Script --- forward modeling wave propagation in 2D with dirichlet boundary condition
%   second order finite difference marching forward scheme
%   time dependent & time independent
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
cmap = customize_colormap(2);

%%  Medium: background velocity + inhomogeneous part
medium.density = 1200;
medium.speed = 1000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
%medium.speed = medium.speed + 500 * (dist < 500);
%medium.speed = medium.speed+(200*sin(0.002*pi*grid.x)+100*cos(0.002*pi*grid.y)).*(abs(grid.x)<=Lx).*(abs(grid.y)<=Ly);
window_type = 'Gaussian';
medium.speed = smooth_source(medium.speed,window_type); % do we need to smooth medium?
figure
imagesc(medium.speed)
colormap(cmap)
title('True Velocity Model')
colorbar
grid.dt = sqrt(0.3*dx^2/max(max(medium.speed))^2);
grid.Nt = ceil(1.2*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.speed))));

%%  Source: either time dependent or initial time
source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; -0.5 -0.5]*1e3;
source_mag = [1;1;1;1;-1]*1e6;
source_index = [3;5];
source_num = size(source_index,1);
source_grid = [round(source_loc(source_index,1)/dx)+(Nx+1)/2 round(source_loc(source_index,2)/dx)+(Ny+1)/2];
source.time_dependent = 0;
if source.time_dependent == 0 % time dependent
    % due to limited memory, only save source-time function
    source.firing_time = [0.0;0.0];
    if size(source.firing_time,1) ~= source_num
        error('Wrong firing time input! Timing should be equal to %d\n', num_source);
    end
    source.peak_freq = 10;
    a = pi^2*source.peak_freq^2;
    t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),source_num,1);
    amplitude = source_mag(source_index) .* exp(-a*(t-source.firing_time).^2); % magnitude?
    p0 = zeros(Nx,Ny);
    for k = 1:source_num
        p0(source_grid(k,1),source_grid(k,2)) = amplitude(k,1);
    end
    source.p0 = smooth_source(p0,window_type);
elseif source.time_dependent == 1 % initial displacement
    source.p0 = zeros(Nx,Ny);
    for j = 1:source_num
        source.p0(source_grid(j,1),source_grid(j,2)) = source_mag(source_index(j));
    end
    window_type = 'Hanning';
    %source.p0 = smooth_source(source.p0,window_type);
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
sensor_data_dirichlet = zeros(size(sensor.mask,2),grid.Nt);

%% time t = dt
Up_past = source.p0;
sensor_data_dirichlet(:,1) = record_data(Up_past, [x_start;y_start], [x_end;y_end]);
Up_now = zeros(Nx,Ny);
Up_now(2:Nx-1,2:Ny-1) = Up_past(2:Nx-1,2:Ny-1) + grid.dt^2*medium.speed(2:Nx-1,2:Ny-1).^2.*...
    ((Up_past(1:Nx-2,2:Ny-1)-2*Up_past(2:Nx-1,2:Ny-1)+Up_past(3:Nx,2:Ny-1))/dx^2 + ...
    (Up_past(2:Nx-1,1:Ny-2)-2*Up_past(2:Nx-1,2:Ny-1)+Up_past(2:Nx-1,3:Ny))/dy^2)/2;
p0 = zeros(size(source.p0));
for k = 1:source_num
    p0(source_grid(k,1),source_grid(k,2)) = amplitude(k,2);
end
p0 = smooth_source(p0,window_type);
Up_now = Up_now + p0;
sensor_data_dirichlet(:,2) = record_data(Up_now, [x_start;y_start], [x_end;y_end]);
Up_future = zeros(Nx,Ny);

%%  start marching forward in time
for j = 3:grid.Nt
    Up_future(2:Nx-1,2:Ny-1) = 2*Up_now(2:Nx-1,2:Ny-1) - Up_past(2:Nx-1,2:Ny-1) +...
        grid.dt^2*medium.speed(2:Nx-1,2:Ny-1).^2.*...
        ((Up_now(1:Nx-2,2:Ny-1)-2*Up_now(2:Nx-1,2:Ny-1)+Up_now(3:Nx,2:Ny-1))/dx^2 + ...
        (Up_now(2:Nx-1,1:Ny-2)-2*Up_now(2:Nx-1,2:Ny-1)+Up_now(2:Nx-1,3:Ny))/dy^2);
    if source.time_dependent == 0
        p0 = zeros(size(source.p0));
        for k = 1:source_num
            p0(source_grid(k,1),source_grid(k,2)) = amplitude(k,j);
        end
        p0 = smooth_source(p0,window_type);
        Up_future = Up_future + p0;
    end
    if mod(j,5) == 0
        imagesc(Up_future)
        title(sprintf('time =%e s',(j-1)*grid.dt));
        caxis([-1e5,1e5])
        colormap(cmap)
        colorbar
        % axis equal
        drawnow
    end
    Up_past = Up_now;
    Up_now = Up_future;
    sensor_data_dirichlet(:,j) = record_data(Up_future, [x_start;y_start], [x_end;y_end]);
end

%%
plot(sensor_data_dirichlet(1,:))
save('data/sensor_data_dirichlet')