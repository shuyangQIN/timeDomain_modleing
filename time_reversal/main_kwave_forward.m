%%  Main Script --- forward modeling wave propagation in 2D with perfectly matched layer
%   reference -- The perfectly matched layer for acoustic waves
%                in absorptive media, Author: Qing-Huo Liu & Jianping Tao
clear
clc
close all
addpath('../source_func')
addpath('../../k-wave-toolbox-version-1.2.1/k-Wave/')

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
dx = 0.010;
dy = 0.010;
Nx = round(2*Lx_comp/dx)+1;
Ny = round(2*Ly_comp/dy)+1;

grid = kWaveGrid(Nx,dx,Ny,dy);

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
medium.sound_speed = 1.000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
%medium.speed = medium.speed + 500 * (dist < 500);
%medium.sound_speed = medium.sound_speed+(0.200*sin(2*pi*grid.x)+0.100*cos(2*pi*grid.y)).*(abs(grid.x)<=Lx).*(abs(grid.y)<=Ly);

% smoothing the medium
window_type = 'Gaussian';
medium.sound_speed = smooth_source(medium.sound_speed,window_type);
figure
imagesc(medium.sound_speed)
colormap(cmap)
title('True Velocity Model')
colorbar

% CFL number: 0.5
grid.dt = sqrt(0.25*dx^2/max(max(medium.sound_speed))^2);
grid.Nt = ceil(1.2*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.sound_speed))));
grid.Nt = round(grid.Nt);

% choose Nt to be bigger than 2 times the longest geodesic

%%  Source: either time dependent or initial time
% specs:
%   source_loc: point source locations on computational grid
%   source_mag: magnitude of point sources
%   source_grid: row/column number of sources in the matrix
%   source.time_dependent: time dependent 0 & independent 1
% source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; 0 0]*1e3;
% source_mag = [1;1;1;1;1]*1e6;
% source_index = [1;2];
% source_num = size(source_index,1);
% source_grid = [round(source_loc(source_index,1)/dx)+(Nx+1)/2 round(source_loc(source_index,2)/dx)+(Ny+1)/2];
% source.time_dependent = 1;
% if source.time_dependent == 0 % time dependent
%     % due to limited memory, only save source-time function
%     source.firing_time = [0.2;0.5]; % point source firing time
%     if size(source.firing_time,1) ~= source_num
%         error('Wrong firing time input! Timing should be equal to %d\n', num_source);
%     end
%     source.peak_freq = 10;  % peak frequency of source wavelet
%     a = pi^2*source.peak_freq^2;
%     t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),source_num,1);
%     amplitude = source_mag(source_index) .* exp(-a*(t-source.firing_time).^2); % magnitude?
%     source.p0 = zeros(Nx,Ny);
% elseif source.time_dependent == 1 % initial displacement
%     source.p0 = zeros(Nx,Ny);
%     for j = 1:source_num
%         source.p0(source_grid(j,1),source_grid(j,2)) = source_mag(source_index(j));
%     end
%     window_type = 'Gaussian';
%     source.p0 = smooth_source(source.p0,window_type);
%     figure
%     imagesc(source.p0)
%     colormap(cmap)
%     title('Initial Displacement')
%     colorbar
%     % due to limited memory, only save source-time function
%     source.firing_time = [0.0;0.0]; % point source firing time
%     if size(source.firing_time,1) ~= source_num
%         error('Wrong firing time input! Timing should be equal to %d\n', num_source);
%     end
%     source.peak_freq = 10;  % peak frequency of source wavelet
%     a = pi^2*source.peak_freq^2;
%     t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),source_num,1);
%     amplitude = source_mag(source_index) .* exp(-a*(t-source.firing_time).^2); % magnitude?
% else
%     error('Unsupported source function!')
% end
phantom_image = phantom('Modified Shepp-Logan',x_end-x_start+1);
source.p0 = zeros(Nx,Ny);
source.p0(x_start:x_end,y_start:y_end) = phantom_image;
imshow(source.p0)
mat_sm = smooth(grid,source.p0);
% window_type = 'Gaussian';
% source.p0 = smooth_source(source.p0,window_type);
% imagesc(source.p0)
% colormap(cmap)
% colorbar

%%  Sensor: location of receivers
sensor.mask = rect_sensor([x_start;y_start], [x_end;y_end], grid.x, grid.y);
%sensor_data_fdm = zeros(size(sensor.mask,2),grid.Nt);

%%
sensor_data_kwave1 = kspaceFirstOrder2D(grid,medium,source,sensor);

%%
plot(sensor_data_kwave1(1,:))
save('data/sensor_data_kwave1')

%%  Down-sampling recorded data
sampling_rate = 1;
sensor_back.mask = sensor.mask(:,1:sampling_rate:end);
sensor_back.time_reversal_boundary_data = sensor_data_kwave1(1:sampling_rate:end,:);
source_back.p0 = zeros(Nx,Ny);
p0_recon = kspaceFirstOrder2D(grid,medium,source_back,sensor_back);
p0_recon = p0_recon(x_start:x_end,y_start:y_end);

%% harmonic extension
lengthx = x_end-x_start+1; lengthy = y_end-y_start+1;
sensor_back.mask = sensor.mask(:,1:sampling_rate:end);
sensor_back.time_reversal_boundary_data = sensor_data_kwave1(1:sampling_rate:end,:);
F = zeros(lengthx-2,lengthy-2);
measurement = sensor_data_kwave1(:,grid.Nt);
U0 = poisson_solver(measurement,lengthx,lengthy,dx,dy,F);
source_back.p0 = zeros(Nx,Ny);
source_back.p0(x_start:x_end,y_start:y_end) = U0;
p1_recon = kspaceFirstOrder2D(grid,medium,source_back,sensor_back);
p1_recon = p1_recon(x_start:x_end,y_start:y_end);

%%
phantom_image = mat_sm(x_start:x_end,y_start:y_end);
figure
imagesc(p0_recon-phantom_image)
colormap(cmap)
colorbar

norm(p0_recon-phantom_image,'fro')

figure
imagesc(p1_recon-phantom_image)
colormap(cmap)
colorbar

norm(p1_recon-phantom_image,'fro')
