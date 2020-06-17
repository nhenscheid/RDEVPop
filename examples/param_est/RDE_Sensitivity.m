% File: RDE_Virtual_Pop_2D.m
% Author: Nick Henscheid 
% Date: 4-2019
%
% Purpose: Solves the 2D reaction diffusion (Fisher-Kolmogorov) problem
%          u_t = div(D*grad(u)) + rho*u*(1-u/k) using the method of lines 
%          and random field parameters (D,rho and k).
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

plotting = 1;  % Display plots during simulation?
%%  Use gpu to compute random fields?  (Marginal speedup for small spatial grid)
if(~isempty(which('gpuDeviceCount')))
    gpu = sign(gpuDeviceCount);
else
    gpu = 0;
end
gpu = 0;  %Or, set gpu manually

%% Set up spatial grid
Ngrid = 256;       % Number of grid points in ea. direction
h     = 1/(Ngrid+1);
xx    = linspace(0,1,Ngrid+2);  % Need extra grid points for ghost cells
[x,y] = meshgrid(xx);  % Uniform grid on [0,1]^2  (Spatial units are cm)

%% Set up time grid
Nt   = 40;        % Number of time points to keep 
tmax = 365;   % Time (days)
t    = linspace(0,tmax,Nt);

%% Set up initial condition
s  = 0.01;    % Std. dev of initial condition
I0 = 5;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y,xc,yc) I0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(x,y,0.5,0.5);
if(plotting)
    thetafig = figure; set(thetafig,'Position',[thetafig.Position(1),thetafig.Position(2),1200,800]);
    subplot(2,2,1);
    imagesc(xx,xx,u0);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.05,0.55,0.4,0.4])
    title(sprintf('Initial Condition (I0 = %f)',I0),'FontSize',14);
end

%% Lumpy isotropic diffusion coefficient D(x,y)
L_D = LumpyBgnd('N',Ngrid+2,'Kbar',3,'b0',1e-7,'cov',0.5,'gpu',gpu);
L_D.israndnumlumps = 0;   % The number of lumps is constant (K == Kbar)
L_D.randomize; 
L_D.b0 = L_D.b0*ones(1,L_D.K); 

if(plotting)
    figure(thetafig);subplot(2,2,2);
    imagesc(xx,xx,L_D.eval);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.55,0.55,0.4,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end

%% Lumpy growth function rho(x,y)
% Construct lumpy background object.  Can set parameters here or later w/
% e.g. L.Kbar, L.cov, etc.
L_rho = LumpyBgnd('N',Ngrid+2,'Kbar',5,'b0',0.25,'cov',0.02,'gpu',gpu);
L_rho.israndnumlumps = 0; 
L_rho.randomize; 
L_rho.b0 = L_rho.b0*ones(1,L_rho.K);
% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    imagesc(xx,xx,L_rho.eval);
    axis image;set(gca,'YDir','normal');colorbar;
    set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Sample growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
end

%% Lumpy carrying cap function
L_kappa = LumpyBgnd('N',Ngrid+2,'Kbar',2,'b0',5e7,'cov',0.5,'gpu',gpu); 
L_kappa.israndnumlumps = 0; 
L_kappa.randomize; 
L_kappa.b0 = L_kappa.b0*ones(1,L_kappa.K); 
K0 = 10;  %NOTE: This isn't the carrying cap, just a scale factor.
if(plotting)
    figure(thetafig); subplot(2,2,4);
    imagesc(xx,xx,L_kappa.eval);axis image;set(gca,'YDir','normal');colorbar;
    kappaaxis = gca;
    set(gca,'Position',[0.55,0.05,0.4,0.4]);hold on;
    title('Sample carrying capacity function $\kappa(x,y)$','FontSize',14);
    xlabel('$x$'); ylabel('$y$'); zlabel('$\kappa(x,y)$ (Units cells/cm$^2$)','FontSize',14);
end


%% Set up RDE Solver
R = RDE; 

R.rho = L_rho; 
R.kappa = L_kappa;
R.D = L_D;
R.u0 = u0;
R.grid = {xx,xx};

%% Solve
n = R.solve(t);

%% Plot one solution path and the estimated tumor burden over time

n.plot(1:10); 
figure; 
plot(t,n.TumorBurden); title('Total tumor burden versus time');


%% Define the initial parameter vector 
%  Each field has 2*K params, namely the center of ea. lump
P = 2*L_D.K + 2*L_rho.K + 2*L_kappa.K; 

beta = [L_D.centers(:);L_rho.centers(:);L_kappa.centers(:)]; 



%% Compute a local sensitivity matrix
[delta_n,delta_N] = local_fin_diff_lumpy(R,t,1:40); 


%% Save
save('results.mat','n','R','rho','kappa','D','N');

fprintf('Done!\n');