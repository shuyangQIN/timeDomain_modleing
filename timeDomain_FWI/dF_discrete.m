function df = dF_discrete(grid,medium,source,sensor,d_res)
%%  PML absorption coefficients
grid.PML_size_x = 20;   % number of grids in PML region along one side of x-direction
grid.PML_size_y = 20;   % number of grids in PML region along one side of y-direction
const = 1e5;
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

%% assemble matrix for Laplace operator
%lap_mat = laplace_mat(grid);
P_interp = interp_mat(grid,sensor);
adj_source = P_interp' * d_res;
speed_vec = reshape(medium.speed,grid.Nx*grid.Ny,1);
speed_vec2 = speed_vec.^2;

%% Initialize Pressure field
%   first time step (time 0) does not contribute to change of medium
%   pressure field initialized in the interior
Up = source.p0(2:grid.Nx-1,2:grid.Nx-1);
%   wave speed in the interior
speed = medium.speed(2:grid.Nx-1,2:grid.Ny-1);
Up_extend = source.p0;
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

dU_old = zeros(grid.Nx*grid.Ny,1);
% second time step (t=dt) contains information via Taylor expansion
dU_now = grid.dt^2*speed_vec.*lap_operator(Up_extend,grid);
df = dU_now.*adj_source(:,2);
%% Marching forward in time
window_type = 'Gaussian';
for j = 2:grid.Nt-1
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
    % compute particle velocity field based on pressure field
    Vx = (-(Up_extend(2:grid.Nx,2:grid.Ny-1) - Up_extend(1:grid.Nx-1,2:grid.Ny-1))/grid.dx + C2_x.*Vx)./C1_x;
    Vy = (-(Up_extend(2:grid.Nx-1,2:grid.Ny) - Up_extend(2:grid.Nx-1,1:grid.Ny-1))/grid.dy + C2_y.*Vy)./C1_y;
    dU_future  = 2*dU_now - dU_old;
    dU_future = dU_future + 2*speed_vec.*lap_operator(Up_extend,grid)*grid.dt^2;
    temp = lap_operator(dU_now,grid);
    dU_future = dU_future + speed_vec2.*temp*grid.dt^2;
    df = df + dU_future.*adj_source(:,j+1);
    dU_old = dU_now;
    dU_now = dU_future;
end
