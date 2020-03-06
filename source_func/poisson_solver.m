function U = poisson_solver(measurement,Nx,Ny,dx,dy,F)
%% Checking Inputs
U = zeros(Nx,Ny);
if size(measurement,1) ~= 2*(Nx+Ny-2)
    error('Dimension of measurement not matching!\n')
end

%% Placing the measurement into matrix
U(1,1:Ny-1) = measurement(1:Ny-1)';
U(1:Nx-1,Ny) = measurement(Ny:Nx+Ny-2);
U(Nx,Ny:-1:2) = measurement(Nx+Ny-1:Nx+2*Ny-3);
U(Nx:-1:2,1) = measurement(Nx+2*Ny-2:2*Nx+2*Ny-4);

%% Finding neighbors of each element
index = reshape((1:(Nx-2)*(Ny-2)),Nx-2,Ny-2);
index_center_row = reshape(index,Nx-2,Ny-2);
index_center_col = reshape(index,Nx-2,Ny-2);
coef_center = 2/dx^2 + 2/dy^2;

index_x_pos_row = reshape(index(1:end-1,:),(Nx-3)*(Ny-2),1);
index_x_pos_col = reshape(index(2:end,:),(Nx-3)*(Ny-2),1);
coef_x_pos = -1/dx^2;

index_x_neg_row = reshape(index(2:end,:),(Nx-3)*(Ny-2),1);
index_x_neg_col = reshape(index(1:end-1,:),(Nx-3)*(Ny-2),1);
coef_x_neg = -1/dx^2;

index_y_pos_row = reshape(index(:,1:end-1),(Nx-2)*(Ny-3),1);
index_y_pos_col = reshape(index(:,2:end),(Nx-2)*(Ny-3),1);
coef_y_pos = -1/dy^2;

index_y_neg_row = reshape(index(:,2:end),(Nx-2)*(Ny-3),1);
index_y_neg_col = reshape(index(:,1:end-1),(Nx-2)*(Ny-3),1);
coef_y_neg = -1/dy^2;

%% Assemble matrix on L.H.S
coef_mat = sparse(index_center_row,index_center_col,coef_center,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
coef_mat = coef_mat + sparse(index_x_pos_row,index_x_pos_col,coef_x_pos,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
coef_mat = coef_mat + sparse(index_x_neg_row,index_x_neg_col,coef_x_neg,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
coef_mat = coef_mat + sparse(index_y_pos_row,index_y_pos_col,coef_y_pos,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2));
coef_mat = coef_mat + sparse(index_y_neg_row,index_y_neg_col,coef_y_neg,(Nx-2)*(Ny-2),(Nx-2)*(Ny-2));

%% Constructing R.H.S vector using boundary conditions
if size(F,1) ~= Nx-2 || size(F,2) ~= Ny-2
    error('Input external term F incorrect\n')
end
F = reshape(F,(Nx-2)*(Ny-2),1);
F(index(1,:)) = F(index(1,:)) + measurement(2:Ny-1)/dx^2;
F(index(:,Ny-2)) = F(index(:,Ny-2)) + measurement(Ny+1:Nx+Ny-2)/dy^2;
F(index(Nx-2,:)) = F(index(Nx-2,:)) + flipud(measurement(Nx+Ny:Nx+2*Ny-3))/dx^2;
F(index(:,1)) = F(index(:,1)) + flipud(measurement(Nx+2*Ny-1:2*Nx+2*Ny-4))/dx^2;

%% Solving Laplace equation 
V = coef_mat\F;
U(2:Nx-1,2:Ny-1) = reshape(V,Nx-2,Ny-2);
