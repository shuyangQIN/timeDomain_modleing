%%  Main Script --- forward modeling wave propagation in 2D with perfectly matched layer
%   reference -- k-wave matlab toolbox
%                using pseudospectral method
clear
clc
close all
addpath('source_func')
addpath('../k-wave-toolbox-version-1.2.1/k-Wave')

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
% grid.x = (-Lx_comp:dx:Lx_comp)'*ones(1,Ny);
% grid.y = ones(Nx,1)*(-Ly_comp:dy:Ly_comp);
% 
% Nx_pml = 20;
% Ny_pml = 20;
Lx = 1000;
Ly = 1000;
grid = kWaveGrid(Nx,dx,Ny,dy);

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
medium.density = 1000;
medium.sound_speed = 1000*ones(size(grid.x));
dist = sqrt(grid.x.^2 + grid.y.^2);
medium.sound_speed = medium.sound_speed + 500 * (dist < 500);
% medium.sound_speed = medium.sound_speed+(200*sin(0.002*pi*grid.x)+100*cos(0.002*pi*grid.y)).*(abs(grid.x)<=Lx).*(abs(grid.y)<=Ly);
% window_type = 'Gaussian';
% medium.sound_speed = smooth_source(medium.sound_speed,window_type); % do we need to smooth medium?
figure
imagesc(medium.sound_speed)
colormap(cmap)
title('True Velocity Model')
colorbar
grid.dt = sqrt(0.25*dx^2/max(max(medium.sound_speed))^2);
grid.Nt = ceil(1.2*sqrt(2)*2*Lx_comp/(grid.dt*min(min(medium.sound_speed))));

%%  Source: either time dependent or initial time
source_loc = [-0.7 0.3; 0.6 -0.3; 0.5 0.5; -0.6 -0.4; 0 0]*1e3;
source_mag = [1;1;1;1;1]*1e6;
source_index = [1;2];
source_num = size(source_index,1);
source_grid = [round(source_loc(source_index,1)/dx)+(Nx+1)/2 round(source_loc(source_index,2)/dx)+(Ny+1)/2];
time_dependent = 1;
if time_dependent == 0 % time dependent
    % due to limited memory, only save source-time function
    source.firing_time = [0.5;1];
    if size(source.firing_time,1) ~= num_source
        error('Wrong firing time input! Timing should be equal to %d\n', num_source);
    end
elseif time_dependent == 1 % initial displacement
    source.p0 = zeros(Nx,Ny);
    for j = 1:source_num
        source.p0(source_grid(j,1),source_grid(j,2)) = source_mag(source_index(j));
    end
    window_type = 'Gaussian';
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

%% forward simulation
sensor_data_kwave = kspaceFirstOrder2D(grid,medium,source,sensor);

%%
save('data/sensor_data_kwave')
plot(sensor_data_kwave(1,:))