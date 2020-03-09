function [f,g] = fh(m0,grid,medium,source,sensor,d_obs)
Nx = grid.x_end - grid.x_start + 1;
Ny = grid.y_end - grid.y_start + 1;
medium.speed(grid.x_start:grid.x_end,grid.y_start:grid.y_end) = reshape(sqrt(1./m0),Nx,Ny);

[d_syn,lap_U] = forward_prop(grid,medium,source,sensor);
d_res =  d_syn - d_obs;
f = 0.5*norm(d_res,'fro')^2;
V = adjoint_wave(grid,medium,source,sensor,d_res);
% V = neumann_adjoint_wave(d_res, medium.speed(grid.x_start:grid.x_end,grid.y_start:grid.y_end), grid.x_end-grid.x_start+1, ...
%     grid.y_end-grid.y_start+1, grid.Nt, grid.dx, grid.dy, grid.dt);
weight = ones(1,1,grid.Nt); weight(:,:,1) = 0.5; weight(:,:,grid.Nt) = 0.5;
df = -sum(weight.*V.*lap_U,3)*grid.dt;
g = reshape(df,size(V,1)*size(V,2),1);

end