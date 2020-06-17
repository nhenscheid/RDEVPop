% File: RDE_Virtual_Pop_2D_homog_bimodal.m
% Author: Nick Henscheid 
% Date: 4-2019, 5-2020
%
% Purpose: Solves the 2D reaction diffusion (Fisher-Kolmogorov) problem
%          u_t = div(D*grad(u)) + rho*u*(1-u/k) using the method of lines 
%          and random homogeneous parameters (D,rho and k).  The
%          distribution for the parameters is assumed to be bimodal, in
%          particular a Gaussian mixture. 
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

% Number of samples to generate
Nsamp = 25;   % Number of samples to generate
Usamp = cell(Nsamp,1); % This is where samples will be stored

% Set up spatial grid
N = 128;       % Number of grid points in ea. direction 
xx = linspace(0,1,N);  
[x,y] = meshgrid(xx);  % Uniform grid on [0,1]^2  (Arbitrary spatial units L)

% Set up time points
Nt = 8;        % Number of time points to keep 
tmax = 1000;   % Arbitrary time units T
t = linspace(0,tmax,Nt);

%% Set up initial condition
s  = 0.015;    % Std. dev of initial condition
I0 = 1;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y) I0*exp(-((x-0.5).^2+(y-0.5).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(x,y);
if(plotting)
    surf(x,y,u0,'EdgeColor','none'); title(sprintf('Initial Condition (I0 = %f)',I0));
end

%%  Diffusion coefficient
D     = 1e-7;   % Diffusion coeff. (Units L^2/T) (assume constant) 

%% Lumpy background growth (bimodal: "low growth" and "aggressive growth")
%  NOTE: you can set the probability of "aggressive" with the parameter
%  prob_agg below.
gpu = gpuDeviceCount>0;
% Construct lumpy background object.  Can set parameters here or later w/
% e.g. L.Kbar, L.cov, etc.
L_rho = LumpyBgnd('N',N,'Kbar',300,'b0',1e-2,'cov',0.001,'gpu',gpu);
% Aggressive variant (note: only b0 is changed!)
L_rho_agg = LumpyBgnd('N',N,'Kbar',300,'b0',2e-2,'cov',0.001,'gpu',gpu); 
rho = @(x,y) reshape(L_rho.eval([x(:),y(:)]),size(x)); 
rho_agg = @(x,y) reshape(L_rho_agg.eval([x(:),y(:)]),size(x)); 
rho_grid = rho(x,y); % Evaluate on the grid
rho_grid_agg = rho_agg(x,y);
prob_agg  = 0.5;   % Probability of being aggressive
% Plot a sample from the growth function
if(plotting)
    fig = figure; pos = get(fig,'Position');
    set(fig,'Position',[pos(1),pos(2),1200,600]);
    subplot(1,2,1);
    imagesc(xx,xx,rho_grid);
    axis image;set(gca,'YDir','normal');colorbar;
    title('Sample growth rate $\rho(x,y)$ (units 1/T)','FontSize',22); 
    subplot(1,2,2);
    imagesc(xx,xx,rho_grid_agg);
    axis image;set(gca,'YDir','normal');colorbar;
    title('Sample aggressive growth rate $\rho(x,y)$ (units 1/T)','FontSize',22); 
end

%% Lumpy carrying cap function (growth is contained within a disk ala "petri dish")
L_kappa = LumpyBgnd('N',N,'Kbar',20,'b0',1,'cov',0.05,'gpu',gpu); 
K0 = 10;
kappa = @(x,y) 0.001+K0*((x-0.5).^2 + (y-0.5).^2 < 0.4^2).*reshape(L_kappa.eval([x(:),y(:)]),size(x));
Kgrid = kappa(x,y); % Evaluate on the grid
figure; 
surf(x,y,Kgrid);
title('Sample carrying capacity function $\kappa(x,y)$','FontSize',22);
xlabel('x'); ylabel('y'); zlabel('$\kappa(x,y)$ (Units cells/$L^2$)','FontSize',22);
%% Set up the ODE right-hand-side
h = 1/(N-1);
odefun = @(T,Y)ODEFun2D(Y,D,rho_grid,Kgrid,h,N);
%% Generate samples 
U = zeros(N,N,Nt);
rho_samp = zeros(N,N,Nsamp);   % Sampled growth fields
kappa_samp = zeros(N,N,Nsamp); % Sampled carrying caps
for isamp=1:Nsamp
    if(plotting)
        progressbar(isamp/Nsamp);
    end
    aggressive = rand()<prob_agg;  % Flip a coin  
    if(aggressive)
        L_rho_agg.randomize;   % Create a new growth function;
        rho_samp(:,:,isamp) = rho_agg(x,y); % Evaluate on the grid
    else
        L_rho.randomize;
        rho_samp(:,:,isamp) = rho(x,y);
    end
    L_kappa.randomize; % Create a new carrying cap
    kappa_samp(:,:,isamp) = kappa(x,y);
    odefun = @(T,Y)ODEFun2D(Y,D,rho_samp(:,:,isamp),kappa_samp(:,:,isamp),h,N);
    tic
    [T,Y] = ode45(odefun,t,u0(:));
    fprintf("Time elapsed = %f\n",toc)

    for i=1:Nt
        U(:,:,i) = reshape(Y(i,:),[N,N]);
    end
    Usamp{isamp} = U;
end

%% Plot the time series for a single sample
kplot = 7;
if(plotting)    
    DisplayImages([2,4],Usamp{kplot});
    figure;surf(x,y,U(:,:,end),'EdgeColor','none')
end

%% Make an array consisting of only the final cell densities
Ufinal    = zeros(N,N,Nsamp);
for i=1:Nsamp
    Ufinal(:,:,i) = Usamp{i}(:,:,end);
end
Ubarfinal = mean(Ufinal,3);
if(plotting)
    figure; imagesc(xx,xx,Ubarfinal);axis image; set(gca,'YDir','normal')
    title('Estimated mean final tumor cell density $\bar{u}(x,y)$');
end
DisplayImages([2,4],Ufinal);

%% Compute the final tumor sizes 
final_size = zeros(Nsamp,1); 

for i=1:Nsamp
    final_size(i) = h^2*sum(sum(Ufinal(:,:,i)));
end

if(plotting)
    figure; histogram(final_size,'normalization','pdf');
    title('Distribution of final tumor sizes $\int_V u(x,T) dx$')
end

