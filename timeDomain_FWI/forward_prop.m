function [sensor_data_fdm,lap_U] = forward_prop(grid,medium,source,sensor)
cmap = customize_colormap(2);
sensor_data_fdm = zeros(size(sensor.mask,2),grid.Nt);
lap_U = zeros(grid.x_end-grid.x_start+1,grid.y_end-grid.y_start+1,grid.Nt);

%%  Finite Difference Time Domain + Absorbing Boundary Conditions
%   Use first order systems to to marching forward in time

%%  PML absorption coefficients
grid.PML_size_x = 20;   % number of grids in PML region along one side of x-direction
grid.PML_size_y = 20;   % number of grids in PML region along one side of y-direction
const = 1e2;
grid.PML_length_x = grid.PML_size_x * grid.dx;
grid.PML_length_y = grid.PML_size_y * grid.dy;
PML_size_x = grid.Lx_comp - grid.PML_length_x;
PML_size_y = grid.Ly_comp - grid.PML_length_y;

% compute PML coefficients on standard computational grid
grid.xPML_std_grid = grid.x(2:grid.Nx-1,2:grid.Ny-1);
grid.yPML_std_grid = grid.y(2:grid.Nx-1,2:grid.Ny-1);
grid.xPML_std_x = (const/grid.PML_length_x)*((abs(grid.xPML_std_grid)>PML_size_x).*(abs(grid.xPML_std_grid)-PML_size_x).^4/grid.PML_length_x^4);
grid.yPML_std_y = (const/grid.PML_length_y)*((abs(grid.yPML_std_grid)>PML_size_y).*(abs(grid.yPML_std_grid)-PML_size_y).^4/grid.PML_length_y^4);

% compute PML coefficients on staggered grid in x-direction
grid.xPML_stag_grid = (grid.x(2:grid.Nx,2:grid.Ny-1) + grid.x(1:grid.Nx-1,2:grid.Ny-1))/2;
grid.xPML_stag_x = (const/grid.PML_length_x)*((abs(grid.xPML_stag_grid)>PML_size_x).*(abs(grid.xPML_stag_grid)-PML_size_x).^4/grid.PML_length_x^4);

% compute PML coefficients on staggered grid in y-direction
grid.yPML_stag_grid = (grid.y(2:grid.Nx-1,2:grid.Ny) + grid.y(2:grid.Nx-1,1:grid.Ny-1))/2;
grid.yPML_stag_y = (const/grid.PML_length_y)*((abs(grid.yPML_stag_grid)>PML_size_y).*(abs(grid.yPML_stag_grid)-PML_size_y).^4/grid.PML_length_y^4);

%% Initialize Pressure field
%   pressure field initialized in the interior
Up = source.p0(2:grid.Nx-1,2:grid.Nx-1);
%   wave speed in the interior
speed = medium.speed(2:grid.Nx-1,2:grid.Ny-1);
Up_extend = source.p0;
lap_U(:,:,1) = (Up_extend(grid.x_start+1:grid.x_end+1,grid.y_start:grid.y_end) - 2*Up_extend(grid.x_start:grid.x_end,grid.y_start:grid.y_end) + Up_extend(grid.x_start-1:grid.x_end-1,grid.y_start:grid.y_end))/grid.dx^2 ...
    + (Up_extend(grid.x_start:grid.x_end,grid.y_start+1:grid.y_end+1) - 2*Up_extend(grid.x_start:grid.x_end,grid.y_start:grid.y_end) + Up_extend(grid.x_start:grid.x_end,grid.y_start-1:grid.y_end-1))/grid.dy^2;
lap_U(:,:,1) = lap_U(:,:,1).*medium.speed(grid.x_start:grid.x_end,grid.y_start:grid.y_end).^2;
sensor_data_fdm(:,1) = record_data(Up_extend, [grid.x_start;grid.y_start], [grid.x_end;grid.y_end]);

%% Initialize Particle Velocity field
%   assume the particle velocity is small enough to be ignored
%   initialize vx and vy at t=0 to be 0

Vx = -(Up_extend(2:grid.Nx,2:grid.Ny-1) - Up_extend(1:grid.Nx-1,2:grid.Ny-1))*grid.dt/(2*grid.dx*medium.density);
Vy = -(Up_extend(2:grid.Nx-1,2:grid.Ny) - Up_extend(2:grid.Nx-1,1:grid.Ny-1))*grid.dt/(2*grid.dy*medium.density);

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
window_type = 'Gaussian';
for j = 2:grid.Nt
    % compute pressure field based on particle velocity field
    Up_x = (-medium.density*speed.^2.*(Vx(2:grid.Nx-1,:)-Vx(1:grid.Nx-2,:))/grid.dx + B2_x.*Up_x)./B1_x;
    Up_y = (-medium.density*speed.^2.*(Vy(:,2:grid.Ny-1)-Vy(:,1:grid.Ny-2))/grid.dy + B2_y.*Up_y)./B1_y;
    Up = Up_x + Up_y;
    % extend the wavefield to compute velocity field
    Up_extend(2:grid.Nx-1,2:grid.Ny-1) = Up;
    p0 = zeros(size(source.p0));
    for k = 1:source.num
        p0(source.grid(k,1),source.grid(k,2)) = source.amplitude(k,j);
    end
    if j <= 100 || source.time_dependent == 0
        p0 = smooth_source(p0,window_type);
    end
    Up_extend = Up_extend + p0;
    lap_U(:,:,j) = (Up_extend(grid.x_start+1:grid.x_end+1,grid.y_start:grid.y_end) - 2*Up_extend(grid.x_start:grid.x_end,grid.y_start:grid.y_end) + Up_extend(grid.x_start-1:grid.x_end-1,grid.y_start:grid.y_end))/grid.dx^2 ...
        + (Up_extend(grid.x_start:grid.x_end,grid.y_start+1:grid.y_end+1) - 2*Up_extend(grid.x_start:grid.x_end,grid.y_start:grid.y_end) + Up_extend(grid.x_start:grid.x_end,grid.y_start-1:grid.y_end-1))/grid.dy^2;
    weight = ones(1,j); weight(1) = 0.5; weight(j) = 0.5;
    amplitude_temp = sum(weight.*source.amplitude(:,1:j),2);
    p1 = zeros(grid.Nx,grid.Ny);
    for k = 1:source.num
        p1(source.grid(k,1),source.grid(k,2)) = amplitude_temp(k);
    end
    if j <= 100 || source.time_dependent == 0
        p1 = smooth_source(p1,window_type);
    end
    lap_U(:,:,j) = lap_U(:,:,j) + p1(grid.x_start:grid.x_end,grid.y_start:grid.y_end);
    lap_U(:,:,j) = lap_U(:,:,j).*medium.speed(grid.x_start:grid.x_end,grid.y_start:grid.y_end).^2;
    % recording wave field at sensor locations
    sensor_data_fdm(:,j) = record_data(Up_extend, [grid.x_start;grid.y_start], [grid.x_end;grid.y_end]);
    % compute particle velocity field based on pressure field
    Vx = (-(Up_extend(2:grid.Nx,2:grid.Ny-1) - Up_extend(1:grid.Nx-1,2:grid.Ny-1))/grid.dx + C2_x.*Vx)./C1_x;
    Vy = (-(Up_extend(2:grid.Nx-1,2:grid.Ny) - Up_extend(2:grid.Nx-1,1:grid.Ny-1))/grid.dy + C2_y.*Vy)./C1_y;
    if mod(j,20) == 0
        imagesc(Up_extend)
        %imagesc(p0)
        title(sprintf('time =%e s',(j-1)*grid.dt));
        caxis([-1e-3,1e-3])
        colormap(cmap)
        colorbar
        %axis equal
        drawnow
    end
end