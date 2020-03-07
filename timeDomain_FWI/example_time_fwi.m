%%  forward modeling wave propagation in 2D with perfectly matched layer
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
Lx_comp = 1300;
Ly_comp = 1300;
dx = 10;
dy = 10;
grid = generate_grid(Lx_comp,Ly_comp,dx,dy);

%%  Colormap for visualization
cmap = customize_colormap(2);

%%  Medium: background velocity + inhomogeneous part
medium.density = 1200;
medium.speed = 1000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
%medium.speed = medium.speed + 500 * (dist < 500);
medium.speed = medium.speed+(200*sin(0.002*pi*grid.x)+100*cos(0.002*pi*grid.y)).*(abs(grid.x)<=grid.Lx).*(abs(grid.y)<=grid.Ly);

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
grid.Nt = ceil(1.5*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.speed))));
% choose total simulation time to be bigger than 2 times the longest geodesic

%%  Source: either time dependent or initial time
% specs:
%   source_loc: point source locations on computational grid
%   source_mag: magnitude of point sources
%   source_grid: row/column number of sources in the matrix
%   source.time_dependent: time dependent 0 & independent 1
source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; 0 0]*1e3;
source_mag = [1;1;1;1;1]*1e6;
source_index = [1;2];
time_dependent = 1;
firing_time = [0.0;0.0];
source  = generate_source(source_loc(source_index,:),source_mag(source_index),firing_time,time_dependent,grid);

%%  Sensor: location of receivers
sensor.mask = rect_sensor([grid.x_start;grid.y_start], [grid.x_end;grid.y_end], grid.x, grid.y);

%%  2D Wave Forward Propagation
sensor_data_fdm = forward_prop(grid,medium,source,sensor);

%%  Down-sampling recorded data
sampling_rate = 1;
d_obs = sensor_data_fdm(1:sampling_rate:end,:);
P_interp = interp_mat(grid,sensor);

%%  Recovering medium via full waveform inversion
medium_recover.density = 1200;
medium_recover.speed = 1000*ones(grid.Nx,grid.Ny);
d_syn = forward_prop(grid,medium_recover,source,sensor);
d_res = d_syn-d_obs;
f_val = 0.5*norm(d_res,'fro')^2;
df = dF(grid,medium_recover,source,sensor,d_res);
df_mat = reshape(df,grid.Nx,grid.Ny);
imagesc(df_mat)
colormap(cmap)
colorbar
