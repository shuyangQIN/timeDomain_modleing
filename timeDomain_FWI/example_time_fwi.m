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
Lx_comp = 1.300;
Ly_comp = 1.300;
dx = 0.010;
dy = 0.010;
grid = generate_grid(Lx_comp,Ly_comp,dx,dy);

%%  Colormap for visualization
cmap = customize_colormap(2);

%%  Medium: background velocity + inhomogeneous part
medium.density = 1.200;
medium.speed = 1.000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
medium.speed = medium.speed + 0.500 * (dist < 0.500);
%medium.speed = medium.speed+(200*sin(0.002*pi*grid.x)+100*cos(0.002*pi*grid.y)).*(abs(grid.x)<=grid.Lx).*(abs(grid.y)<=grid.Ly);

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
source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; 0 0];%*1e3;
source_mag = [1;1;1;1;1];%*1e6;
source_index = [5];
time_dependent = 0;
firing_time = [0.2];
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
medium_recover.density = 1.200;
medium_recover.speed = 1.000*ones(grid.Nx,grid.Ny);
% % time independent update trial -- failed
Nx = grid.x_end - grid.x_start + 1;
Ny = grid.y_end - grid.y_start + 1;
m0 = reshape(1./medium_recover.speed(grid.x_start:grid.x_end,grid.y_start:grid.y_end).^2,Nx*Ny,1);
fm = @(x) fh(x,grid,medium_recover,source,sensor,d_obs);
mn_r = lbfgs_wWolfe_modified(fm,m0);
medium_recover.speed = reshape(sqrt(1./mn_r),Nx,Ny);
m0 = mn_r;

% [d_syn,lap_U] = forward_prop(grid,medium_recover,source,sensor);
% d_res = d_syn - d_obs;
% V = adjoint_wave(grid,medium,source,sensor,d_res);
% weight = ones(1,1,grid.Nt); weight(:,:,1) = 0.5; weight(:,:,grid.Nt) = 0.5;
% df = sum(weight.*V.*lap_U,3);
% imagesc(df)
% colormap(cmap)
% colorbar
% g = reshape(df,size(V,1)*size(V,2),1);

%% 
cmap = customize_colormap(1);
imagesc(medium_recover.speed)
colormap(cmap)
colorbar
