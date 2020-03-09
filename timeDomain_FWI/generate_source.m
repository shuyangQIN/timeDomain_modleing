function source  = generate_source(source_loc,source_mag,firing_time,time_dependent,grid)
cmap = customize_colormap(1);
source.loc = source_loc;
source.mag = source_mag;
source.time_dependent = time_dependent;
source.num = size(source.loc,1);
source.grid = [round(source.loc(:,1)/grid.dx)+(grid.Nx+1)/2 round(source_loc(:,2)/grid.dy)+(grid.Ny+1)/2];
if source.time_dependent == 0 % time dependent
    % due to limited memory, only save source-time function
    if size(firing_time,1) ~= source.num
        error('Wrong firing time input! Timing should be equal to %d\n', num_source);
    end
    source.firing_time = firing_time; % point source firing time
    source.peak_freq = 10;  % peak frequency of source wavelet
    a = pi^2*source.peak_freq^2;
    t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),source.num,1);
    amplitude = source.mag .* ricker_wavelet(t,source.firing_time,sqrt(1/a));
    %amplitude = source.mag .* exp(-a*(t-source.firing_time).^2); % magnitude?
    source.amplitude = amplitude;
    source.p0 = zeros(grid.Nx,grid.Ny);
elseif source.time_dependent == 1 % initial displacement
    source.p0 = zeros(grid.Nx,grid.Ny);
    for j = 1:source.num
        source.p0(source.grid(j,1),source.grid(j,2)) = source.mag(j);
    end
    window_type = 'Gaussian';
    source.p0 = smooth_source(source.p0,window_type);
    figure
    imagesc(source.p0)
    colormap(cmap)
    title('Initial Displacement')
    colorbar
    % due to limited memory, only save source-time function
    source.firing_time = [0.0;0.0]; % point source firing time
    if size(source.firing_time,1) ~= source.num
        error('Wrong firing time input! Timing should be equal to %d\n', source.num);
    end
    source.peak_freq = 10;  % peak frequency of source wavelet
    a = pi^2*source.peak_freq^2;
    t = repmat((0:grid.dt:(grid.Nt-1)*grid.dt),source.num,1);
    amplitude = source.mag .* exp(-a*(t-source.firing_time).^2); % magnitude?
    source.amplitude = amplitude;
else
    error('Unsupported source function!')
end