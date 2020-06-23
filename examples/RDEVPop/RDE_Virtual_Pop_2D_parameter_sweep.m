% File: RDE_Virtual_Pop_2D_parameter_sweep.m
% Author: Nick Henscheid 
% Date: 6-2020
%
% Purpose: Solves the 2D reaction diffusion (Fisher-Kolmogorov) problem
%          u_t = div(D*grad(u)) + rho*u*(1-u/k) using the method of lines 
%          and random field coefficients (D,rho and k).
%          The random field parameters are then swept to obtain estimates
%          of the probability distribution of N(t_f) 
% Notes:   The Laplacian is discretized using the standard 2nd order 
%          5-point stencil. The resulting ODE system is solved using ode45.
%          Both of these discreizations are sub-optimal, but tend to work
%          in reasonable parameter regimes.
%          Units of u(x,t) are cells/cm^2 with t in days (2D)
%          Units of D are thus cm^2/day
%          Units of alpha are cells/day
%          Units of k are cells/cm^2 

% See documentation for full details. 

%%
close all
%% Set up sample storage
% Number of samples to generate
Nsamp = 64;
%% Set up spatial grid
Ngrid = 128;       % Number of grid points in ea. direction
h = 1/(Ngrid+1);
xx = linspace(0,1,Ngrid+2);  
[x,y] = meshgrid(xx);  % Uniform grid on [0,1]^2  (Spatial units are cm)
%% Set up time grid
Nt = 2;        % Number of time points to keep (only need 2 for this study)
tmax = 365;   % Time (days)
t = linspace(0,tmax,Nt);
%%
disp('Generating initial condition and random coefficient fields');
%% Set up initial condition
s  = 0.01;    % Std. dev of initial condition
I0 = 5;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y,xc,yc) I0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(x,y,0.5,0.5);
%% Lumpy isotropic diffusion coefficient D(x,y)
L_D     = LumpyBgnd('N',Ngrid+2,'Kbar',20,'b',1e-7,'cov',0.04);
%% Lumpy growth function rho(x,y)
nrho_samples = 4; 
rhob         = linspace(0.1,1,nrho_samples); 
rhocov       = linspace(0.001,0.003,nrho_samples); 
L_rho   = LumpyBgnd('N',Ngrid+2,'Kbar',200,'b',0.25,'cov',0.002);
%% Lumpy carrying cap function
L_kappa = LumpyBgnd('N',Ngrid+2,'Kbar',100,'b',5e7,'cov',0.1); 
%% Set up RDE Solver
R = RDE; 
R.rho = L_rho; 
R.kappa = L_kappa;
R.D = L_D;
R.u0 = u0;
R.grid = {xx,xx};
%% rho parameter sweep
Nfinal = zeros(nrho_samples,nrho_samples,Nsamp);
for irhob = 1:nrho_samples
    fprintf('irhob = %i/%i\n',irhob,nrho_samples); 
    L_rho.b = rhob(irhob); 
    for irhocov = 1:nrho_samples
        fprintf('irhocov = %i/%i\n',irhocov,nrho_samples); 
        L_rho.cov = rhocov(irhocov);
        R.rho = L_rho; 
        [n,rho,kappa,D] = R.Sample(Nsamp,t);
        for i=1:Nsamp
            N = n{i}.TumorBurden;
            Nfinal(irhob,irhocov,i) = N(end); 
        end
    end
end
%% Save
save('Results.mat',Nfinal); 

fprintf('Done!\n');