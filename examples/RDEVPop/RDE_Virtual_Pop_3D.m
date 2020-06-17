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

close all

plotting = 1;  % Display plots during simulation?

%%  Use gpu to compute random fields?
if(~isempty(which('gpuDeviceCount')))
    gpu = sign(gpuDeviceCount);
else
    gpu = 0;
end
%gpu = 0;  %Or, set gpu manually

%% Number of samples to generate
Nsamp  = 1;   % Number of samples to generate
Usamp  = cell(Nsamp,1); % This is where samples will be stored
Ufinal = cell(Nsamp,1); % Storage for the final tumor 

%% Set up spatial grid
N = 64;       % Number of grid points in ea. direction 
S = RectSupport([0,1;0,1;0,1]);  % Unit cube support function
X = unifGrid(S,[N,N,N]);
x = X(:,:,:,1);
y = X(:,:,:,2);
z = X(:,:,:,3);

%% Set up time points
Nt = 8;        % Number of time points to keep 
tmax = 1000;   % Arbitrary time units T
t = linspace(0,tmax,Nt);

%% Set up initial condition
s  = 0.05;    % Std. dev of initial condition
sigma = diag(s^2*ones(3,1));
I0 = 5;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y,z) I0*reshape(mvnpdf([x(:),y(:),z(:)],[0.5,0.5,0.5],sigma),size(x));
u0 = initcond(x,y,z);
if(plotting)
    DisplayImages([8,8],u0,'title','Initial Condition','clim',[0,100]);
end

%%  Diffusion coefficient (constant in 3D for now!!)
D     = 1e-7;   % Diffusion coeff. (Units L^2/T) (assume constant) 

%% Lumpy background growth 
% Construct lumpy background object.  Can set parameters here or later w/
% e.g. L.Kbar, L.cov, etc.

L_rho = LumpyBgnd('S',S,'N',N,'Kbar',200,'b0',5e-3,'cov',0.01,'gpu',gpu);
rho = @(x,y,z) reshape(L_rho.eval([x(:),y(:),z(:)]),size(x)); 
rho_grid = rho(x,y,z); % Evaluate on the grid
% Plot a sample from the growth function
if(plotting)
    figure;
    plot(L_rho);
    title(sprintf('Sample growth rate $\\rho(x,y,z)$ (units 1/T) ($\\rho = %2.2e$ isosurface)',mean(rho_grid(:))),'FontSize',22); 
end

%% Lumpy carrying cap function (growth is contained within a sphere)
L_kappa = LumpyBgnd('S',S,'N',N,'Kbar',20,'b0',1e8,'cov',0.05,'gpu',gpu); 
K0 = 10;
kappa = @(x,y,z) 0.001+K0*((x-0.5).^2 + (y-0.5).^2 + (z-0.5).^2 < 0.4^2).*reshape(L_kappa.eval([x(:),y(:),z(:)]),size(x));
Kgrid = kappa(x,y,z); % Evaluate on the grid
if(plotting)
    figure;
    plot(L_kappa)
    title(sprintf('Sample carrying capacity function $\\kappa(x,y,z)$ ($\\kappa = %2.2e$ isosurface)',mean(Kgrid(:))),'FontSize',22);
    xlabel('x','FontSize',22,'Interpreter','latex'); ylabel('y','FontSize',22,'Interpreter','latex'); 
    zlabel('z','FontSize',22,'Interpreter','latex');
end

%% Set up the ODE right-hand-side
h = 1/(N-1);
odefun = @(T,Y)ODEFun3D(Y,D,rho_grid,Kgrid,h,N);
%% Generate samples 
U = zeros(N,N,N,Nt);
rho_samp = zeros(N,N,N,Nsamp);   % Sampled growth fields
kappa_samp = zeros(N,N,N,Nsamp); % Sampled carrying caps
fprintf('Computing 3D RDE Solutions (may take O(minutes) per sample!) (Nsamp = %i, N_ode = %i)\n',Nsamp,N^3);
for isamp=1:Nsamp
    if(plotting)
        progressbar(isamp/Nsamp);
    end
    fprintf('Randomizing L_rho\n')
    L_rho.randomize;   % Create a new growth function;
    rho_samp(:,:,:,isamp) = rho(x,y,z); % Evaluate on the grid
    fprintf('Randomizing L_kappa\n')
    L_kappa.randomize; % Create a new carrying cap
    kappa_samp(:,:,:,isamp) = kappa(x,y,z);
    odefun = @(T,Y)ODEFun3D(Y,D,rho_samp(:,:,:,isamp),kappa_samp(:,:,:,isamp),h,N);
    tic
    fprintf('Solving ODE system\n')
    [T,Y] = ode45(odefun,t,u0(:));
    fprintf("Time elapsed = %f\n",toc)

    for i=1:Nt
        U(:,:,:,i) = reshape(Y(i,:),[N,N,N]);
    end
    Usamp{isamp} = U;
    Ufinal{isamp} = U(:,:,:,end);
end

%% Plot the time series for a single sample
kplot = 1;
Uplot = Usamp{kplot};
isovalue = 1e3;
if(plotting)  
    figure;
    pat = patch(isosurface(x,y,z,Uplot(:,:,:,end),isovalue));
    pat.FaceColor = [0 0.4471 0.7412];
    pat.EdgeColor = [0,0,0];
    pat.FaceAlpha = 0.75;
    pat.EdgeAlpha = 0.5;
    
end

%% Show a movie
kplot = 1;
Uplot = Usamp{kplot};
fig = figure;
set(fig,'Position',[fig.Position(1),fig.Position(2),1000,1000]);
isovalue = 20;
for i=1:Nt
    isof3 = isosurface(x,y,z,Uplot(:,:,:,i),isovalue);
    isop3 = patch(isof3);
    set(isop3,'EdgeColor',[0,0,0],'FaceColor',[1,0,0],'FaceAlpha',0.5);
    title(sprintf('Isosurface u = %f, t = %f',isovalue,t(i)));
    pause(1);
end


%% Plot a series of isosurfaces
kplot = 1;
Uplot = Usamp{kplot};
isovals = logspace(1,log10(max(Uplot(:))/2),9);
figure; 
for i = 1:3
    for j=1:3
        subplot(3,3,(i-1)*3+j);
        isop=patch(isosurface(x,y,z,Uplot(:,:,:,end),isovals((i-1)*3+j)));
        set(isop,'EdgeColor',[0,0,0],'FaceColor',[1,0,0],'FaceAlpha',0.6,'EdgeAlpha',0.9);
        set(gca,'Clipping','off');
        xlim([0,1]);ylim([0,1]);zlim([0,1]);
        title(sprintf('u = %2.2e level set',isovals((i-1)*3+j)))
    end
end

%% Compute the final tumor sizes (units: # of cells)
final_size = zeros(Nsamp,1); 

for i=1:Nsamp
    final_size(i) = h^3*sum(Ufinal{i}(:));  % Basic 2D quadrature. 
end

if(plotting)
    figure; histogram(final_size,'normalization','pdf');
    title('Distribution of final tumor sizes (total \# of cells)')
end

