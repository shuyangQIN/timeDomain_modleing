%--------------------------------------------------------------------%
%------------ 1D FDTD method ----------------------------------------%
%------------ with PML Boundary conditions---------------------------%
%------------ all units in eV,fs, nm --------------------------------%
%--------------------------------------------------------------------%
clc
close all
clear all
%--------------------------------------------------------------------%
h=0.6582119514 ; % Planck's constant
dt=0.13; % time step in femto seconds
tmax=100; % maximum time limit
t=-20:dt:tmax; % time array
c=299.792458; % speed of light in nm/fs
S=1; % courant factor
dz=c*dt/S;
mu0=2.013354451e-4; % permeability of free space in (V fs^2/e nm)
ep0=55.26349597e-3; % permittivity of free space in (e / V nm)
%!-----------------create E-field array
wlev=1.5d0 ;   % freq of light in eV
wl=(wlev/h) ;  % freq of light in (1/fs)
delta=10.0d0 ; % width of electric field
source=exp(-((t).^2)/(delta^2)).*cos(wl.*t); % Electric field is an Gaussian envelop
%--------------------------------------
a1=dt/dz; % intermediate term for computation
%---Gaussian envelop source----------
Zs=200; % position index of source
M=500; % no. of spatial grid points
Z=(0:M-1).*dz; % space axis
ep(1:M)=ep0; % permitivity array
mu(1:M)=mu0; % pemeability array
%---- PML absorbing boundary condition----
sigma(1:M)=0; % initialize conductivity array
d=55; % width of PML layer
m=3; % polynomial order for grading sigma array (pp 292, Taflove)
neta=sqrt(mu0/ep0); 
R=1e-8; % required reflectivity
sigma_max=-(m+1)*log(R)/(2*neta*d*dz);
Pright=((1:d+1)./d).^m*sigma_max;
sigma(M-d:M)=Pright; % lossy conductivity profile
sigma(1:d+1)=fliplr(Pright);
sigma_star(1:M)=sigma.*mu0./ep0; % Eq 7.8 Taflove, pp 275
%------------- PML constants ----------------------------------------%
A=((mu-0.5*dt*sigma_star)./(mu+0.5*dt*sigma_star)); 
B=(dt/dz)./(mu+0.5*dt*sigma_star);                          
C=((ep-0.5*dt*sigma)./(ep+0.5*dt*sigma)); 
D=(dt/dz)./(ep+0.5*dt*sigma);                     
%---------------------------------------------------------------------%
%------initialize fields in space array
Hy(1:M)=0.0; 
Ex(1:M)=0.0;
%---- begin time loop----
fh = figure(1);
set(fh, 'Color', 'white'); 
for n=1:length(t)
        Ex(Zs)=source(n); % insert source in certain space grid
        Hy(1:M-1)=A(1:M-1).*Hy(1:M-1)-B(1:M-1).*(Ex(2:M)-Ex(1:M-1));
        Ex(2:M-1)=C(2:M-1).*Ex(2:M-1)-D(2:M-1).*(Hy(2:M-1)-Hy(1:M-2));
        Ex(M)=Ex(M-1);
        figure(1)
        plot(Z,Ex)
        xlabel('Z in nm','fontSize',14);
        ylabel('E_x','fontSize',14);  
        title('1D FDTD with PML','fontSize',14);              
        axis([0 M*dz -1 1]);
        getframe();
end
%-----------------------------------------------------------------------%