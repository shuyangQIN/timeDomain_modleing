function P_interp = interp_mat(grid,sensor)

LL = (grid.Nx + 1)/2;
UL_index = [floor(sensor.mask(1,:)/grid.dx) + LL; floor(sensor.mask(2,:)/grid.dy) + LL];
UR_index = [floor(sensor.mask(1,:)/grid.dx) + LL; ceil(sensor.mask(2,:)/grid.dy) + LL];
DL_index = [ceil(sensor.mask(1,:)/grid.dx) + LL; floor(sensor.mask(2,:)/grid.dy) + LL];
DR_index = [ceil(sensor.mask(1,:)/grid.dx) + LL; ceil(sensor.mask(2,:)/grid.dy) + LL];

% UL_loc = [floor(sensor.mask(1,:)/grid.dx)*grid.dx; floor(sensor.mask(2,:)/grid.dy)*grid.dy];
% UR_loc = [floor(sensor.mask(1,:)/grid.dx)*grid.dx; ceil(sensor.mask(2,:)/grid.dy)*grid.dy];
% DL_loc = [ceil(sensor.mask(1,:)/grid.dx)*grid.dx; floor(sensor.mask(2,:)/grid.dy)*grid.dy];
% DR_loc = [ceil(sensor.mask(1,:)/grid.dx)*grid.dx; ceil(sensor.mask(2,:)/grid.dy)*grid.dy];

dL = sensor.mask(2,:) - floor(sensor.mask(2,:)/grid.dy)*grid.dy;
dR = ceil(sensor.mask(2,:)/grid.dy)*grid.dy - sensor.mask(2,:);
dU = sensor.mask(1,:) - floor(sensor.mask(1,:)/grid.dx)*grid.dx;
dD = ceil(sensor.mask(1,:)/grid.dx)*grid.dx - sensor.mask(1,:);

for j = 1:size(sensor.mask,2)
    if dU(j) == 0 && dD(j) == 0
        dU(j) = grid.dx/2; dD(j) = grid.dx/2;
    end
    if dL(j) == 0 && dR(j) == 0
        dL(j) = grid.dy/2; dR(j) = grid.dy/2;
    end
end

wUL = dD.*dR/(grid.dx*grid.dy);
wUR = dD.*dL/(grid.dx*grid.dy);
wDL = dU.*dR/(grid.dx*grid.dy);
wDR = dU.*dL/(grid.dx*grid.dy);

UL_index_vec = (UL_index(2,:)-1)*grid.Nx + UL_index(1,:);
UR_index_vec = (UR_index(2,:)-1)*grid.Nx + UR_index(1,:);
DL_index_vec = (DL_index(2,:)-1)*grid.Nx + DL_index(1,:);
DR_index_vec = (DR_index(2,:)-1)*grid.Nx + DR_index(1,:);

P_interp = sparse(repmat(1:size(sensor.mask,2),1,4),[UL_index_vec UR_index_vec DL_index_vec DR_index_vec],...
    [wUL wUR wDL wDR],size(sensor.mask,2),grid.Nx*grid.Ny);

end
