% File: RDE_Virtual_Pop_2D_homog.m
% Author: Nick Henscheid 
% Date: 4-2019, 11-2019, 5-2020
%
% Purpose: Solves the 2D reaction diffusion (Fisher-Kolmogorov) problem
%          u_t = div(D*grad(u)) + rho*u*(1-u/k) using the method of lines 
%          and random homogeneous parameters (D,rho and k).
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
gpu = 0;     % Or, set gpu manually
%% Number of samples to generate
Nsamp = 1;   %
%% Set up spatial grid
Ngrid = 256;               % Number of primary grid points in ea. direction
x = linspace(0,1,Ngrid+2); % NOTE: we need Ngrid+2 points for ghost cells
[xg,yg] = meshgrid(x);     % Uniform grid on [0,1]^2  (Spatial units are cm)
%% Set up time grid
Nt = 50;                    % Number of time points to keep 
tmax = 3*365;               % Time to simulate (days)
t = linspace(0,tmax,Nt);
%% Set up initial condition
I0 = 5;       % Integral of intial cell density ("number of initial cells")
s  = 0.01;    % Std. dev of initial condition
initcond = @(x,y,xc,yc) I0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(xg,yg,0.5,0.5);
if(plotting)
    thetafig = figure; 
    set(thetafig,'Position',[thetafig.Position(1),thetafig.Position(2),1200,800]);
    subplot(2,2,1);
    imagesc(x,x,u0);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.05,0.55,0.4,0.4])
    title(sprintf('Initial Condition (I0 = %f)',I0),'FontSize',14);
end
%% Uniform isotropic diffusion coefficient D(x,y)
L_D = UnifRandField('N',Ngrid+2,'a',1e-7,'b',5e-7); 

if(plotting)
    figure(thetafig);subplot(2,2,2);
    plot(L_D); colorbar; set(gca,'CLim',[L_D.a,L_D.b]);
    set(gca,'Position',[0.55,0.55,0.4,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end
%% Uniform growth function rho(x,y)
L_rho = UnifRandField('N',Ngrid+2,'a',0.2,'b',0.5);

% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    plot(L_rho); colorbar; set(gca,'CLim',[L_rho.a,L_rho.b]);
    set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Sample growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
end
%% Lumpy carrying cap function
L_kappa = UnifRandField('N',Ngrid+2,'a',1e7,'b',5e8);
if(plotting)
    figure(thetafig); subplot(2,2,4);
    plot(L_kappa);colorbar;set(gca,'CLim',[L_kappa.a,L_kappa.b]);
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
R.grid = {x,x};
%% Solve
[n,rho,kappa,D] = R.sample(Nsamp,t);

N = zeros(Nt,Nsamp);
for i=1:Nsamp
    N(:,i) = n{i}.TumorBurden;
end

%% Plot one solution path and the estimated tumor burden over time

n{1}.plot(1:10); 
figure; 
plot(N); title('Total tumor burden versus time');

%% Save
save('results_homog.mat','n','R','rho','kappa','D','N');

fprintf('Done!\n');