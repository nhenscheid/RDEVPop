% Demonstration of Markov Chain Monte Carlo reconstruction of a lumpy-type
% random field from noisy image data.
% The imaging system is assume to be modeled by a Gaussian PSF with Poisson
% 
close all;
plotting = 1; 

%% Spatial grid
Ngrid = 512;       % Number of grid points in ea. direction
h = 1/(Ngrid+1);
xx = linspace(0,1,Ngrid+2);  
[x,y] = meshgrid(xx);  % Uniform grid on [0,1]^2  (Spatial units are cm)

%% Time grid
Nt = 20;        % Number of time points to keep 
tmax = 365;   % Time (days)
t = linspace(0,tmax,Nt);

%% Init cond 
s  = 0.01;    % Std. dev of initial condition
I0 = 5;        % Integral of intial cell density ("number of initial cells")
initcond = @(x,y,xc,yc) I0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
u0 = initcond(x,y,0.5,0.5);
if(plotting)
    thetafig = figure; set(thetafig,'Position',[thetafig.Position(1),thetafig.Position(2),1200,800]);
    subplot(2,2,1);
    imagesc(xx,xx,u0);axis image; set(gca,'YDir','normal'); colorbar;
    %set(gca,'Position',[0.05,0.55,0.3,0.4])
    title(sprintf('Initial Condition (I0 = %f)',I0),'FontSize',14);
end

%% Lumpy isotropic diffusion coefficient D(x,y)
L_D = LumpyBgnd('N',Ngrid+2,'Kbar',20,'b',1e-7,'cov',0.04);

if(plotting)
    figure(thetafig);subplot(2,2,2);
    imagesc(xx,xx,L_D.Eval);axis image; set(gca,'YDir','normal'); colorbar;
    %set(gca,'Position',[0.55,0.55,0.3,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end

%% Lumpy growth function rho(x,y)
% Construct lumpy background object.  Can set parameters here or later w/
% e.g. L.Kbar, L.cov, etc.
L_rho = LumpyBgnd('N',Ngrid+2,'Kbar',200,'b',0.1,'cov',0.002);

% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    imagesc(xx,xx,L_rho.Eval);
    axis image;set(gca,'YDir','normal');colorbar;
    %set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Original growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
    %figure(thetafig);subplot(2,3,5); imagesc(xx,xx,L_rho_MLE.Eval);
    axis image;set(gca,'YDir','normal');colorbar;
    %set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('MLE growth rate $\hat{\rho}(x,y)$ (units cells/day)','FontSize',14); 
end
%% Lumpy carrying cap function
L_kappa = LumpyBgnd('N',Ngrid+2,'Kbar',100,'b',5e7,'cov',0.1); 
K0 = 10;  %NOTE: This isn't the carrying cap, just a scale factor.
if(plotting)
    figure(thetafig); subplot(2,2,4);
    imagesc(xx,xx,L_kappa.Eval);axis image;set(gca,'YDir','normal');colorbar;
    kappaaxis = gca;
    %set(gca,'Position',[0.55,0.05,0.4,0.4]);hold on;
    title('Sample carrying capacity function $\kappa(x,y)$','FontSize',14);
    xlabel('$x$'); ylabel('$y$'); zlabel('$\kappa(x,y)$ (Units cells/cm$^2$)','FontSize',14);
    drawnow; 
end
%% Set up RDE Solver

% With the `original' version of rho
R0 = RDE; 
R0.rho = L_rho; 
R0.kappa = L_kappa;
R0.D = L_D;
R0.u0 = u0;
R0.grid = {xx,xx};


%% Solve to get  `true' tumor profile

n0 = R0.Solve(t); 
plot(n0,1:10);

%%
Mx = 64;
My = 64;

image_FWHM  = 0.05;  % cm (0.1=1mm)
image_sigma = image_FWHM/(2*sqrt(2*log(2)));
image_amp   = 1/(2*pi*image_sigma^2);
blur_kernel = @(x,y,x0,y0) image_amp*exp(-(1/(2*image_sigma^2))*((x-x0).^2 + (y-y0).^2));
quantum_yield = 1; % Number of photons detected per cell
tumor_scale_factor = 1; 

image_time = 5; 

idx = 190:2:316;
idy = 190:2:316; 

%n0.cell_density = n0.cell_density./tumor_scale_factor;

tic
[g,gbar,h_mat] = compute_gaussian_image_tumor(n0,blur_kernel,image_time,quantum_yield,idx,idy);
fprintf('Time to compute image = %f\n',toc); 

%
xplot = linspace(0,1,Ngrid);
subplot(2,2,1);
imagesc(xplot,xplot,n0.cell_density(:,:,image_time)); set(gca,'YDir','normal'); axis image; title('Original tumor'); 

subplot(2,2,2);
imagesc(h_mat);axis image;set(gca,'YDir','normal');title('Blur Kernel');
subplot(2,2,3);
imagesc(gbar);axis image;colorbar;set(gca,'YDir','normal');title('Mean image');
subplot(2,2,4);
imagesc(g);axis image;colorbar;set(gca,'YDir','normal');title('Noisy image');
%%
L = LumpyBgnd;
L.N = 64;
L.Kbar = 100;
%L.cov  = L0.cov;                                                          
L.SetPadFactor(0);
L.israndnumlumps = 0; 
L.gpu = 0; 
L.TurnOffWarnings; 
%xmin = -3*sqrt(L.cov(1,1));
xmin = 0; 

K_recon = 10^2;

[cx,cy] = meshgrid(linspace(0.45,0.6,sqrt(K_recon)));

L.centers = [cx(:),cy(:)];

% Use this to start at a random initial condition
% L.centers = xmin + (1-2*xmin)*rand(L0.K,2);

% Use this to start at the MLE point 


% Optimize centers, b and single sigma
%X0 = [0.2+0.5*rand(2*K_recon,1);1e4*ones(K_recon,1);1e-3];
% Optimize centers, single b and single sigma
%X0 = [0.4+0.2*rand(2*K_recon,1);1e4;1e-4];
% Optimize centers and b, fix sigma
%X0 = [0.4+0.2*rand(2*K_recon,1);1e4];
% dx = 1/(K_recon+1); 
% [x_init,y_init] = meshgrid(dx:dx:(1-dx)); 
% X0 = [x_init(:);y_init(:);rand(2,1)];  

% Optimize only b, fix sigma and fix centers on a grid
X0 = [1e6*ones(K_recon,1)];

%L.centers = [X0(1:K_recon),X0((K_recon+1):2*K_recon)]; 
%L.b       = X0((2*K_recon+1):3*K_recon);
L.b = X0(1:K_recon);
%L.cov     = 0.0005; 
L.cov = (L.centers(2,2) - L.centers(1,2))^2/8;
optplotfig = figure; 
subplot(2,3,1); imagesc(n0.cell_density(:,:,image_time)); set(gca,'YDir','normal'); axis image; colorbar; 
subplot(2,3,2); plot(L);  colorbar; 
subplot(2,3,4); imagesc(gbar);set(gca,'YDir','normal'); axis image;colorbar; title('Original gbar'); 
subplot(2,3,5); imagesc(g);set(gca,'YDir','normal'); axis image; drawnow;    title('Original g'); 

objfun = @(X) poiss_obj_fun_tumor_centers_fixed(g,X,L,blur_kernel,quantum_yield,idx,idy);
maxfeval = 50000; 
maxiter  = 500; 

options = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations',maxfeval,'MaxIterations',maxiter,'Algorithm','interior-point','OptimalityTolerance',1e-8,'OutputFcn',@(x,optimValues,state)opt_plot_fun_tumor_centers_fixed(x,optimValues,state,optplotfig,n0,L,image_time,blur_kernel,idx,idy));
%options = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations',maxfeval,'MaxIterations',maxiter,'Algorithm','interior-point','OptimalityTolerance',1e-8);

%LB = [0.2*ones(2*K_recon,1);zeros(K_recon,1);0];
%LB = [0.2*ones(2*K_recon,1);0;0];
%LB = [0.2*ones(2*K_recon,1);0];
%UB = [0.8*ones(2*K_recon,1);1e9];
%UB = [0.8*ones(2*K_recon,1);1e9;.005]; 
%UB = [0.8*ones(2*K_recon,1);1e9*ones(K_recon,1);.1]; 
LB = zeros(K_recon,1); 
UB = [1e9*ones(K_recon,1)]; 
optTime = tic;
Xstar = fmincon(objfun,X0,[],[],[],[],LB,UB,[],options);
fprintf('Time elapsed = %f\n',toc(optTime)); 

%%
%n0.cell_density = tumor_scale_factor*n0.cell_density
L.b       = Xstar;
L.N = 512; % Can increase # eval points for display
[~,gbarstar,~] = compute_gaussian_image_lumpy(L,blur_kernel); 
figure; 
subplot(2,3,1); imagesc(n0.cell_density(:,:,image_time)); set(gca,'YDir','normal'); axis image;  title('Original object'); 
subplot(2,3,2); plot(L);  title('Approx. maximum likelihood estimate') 
subplot(2,3,3); imagesc(gbarstar); axis image; colorbar; set(gca,'YDir','normal'); title('Mean image (MLE)'); 
subplot(2,3,4); imagesc(gbar); axis image; colorbar; set(gca,'YDir','normal'); title('Mean image (orig.)'); 
subplot(2,3,5); 
imagesc(g); axis image; colorbar; set(gca,'YDir','normal'); title('Noisy image'); 

%%
FileNameAndLocation='LumpyMCMC';
dtime = datestr(now,'dd-mm-HHMMPM'); 
newbackup=sprintf('./savescripts/%s_%s.m',FileNameAndLocation,dtime);
currentfile=sprintf('./%s.m',FileNameAndLocation);
copyfile(currentfile,newbackup);
savefile=sprintf('./results/%s_%s.mat',FileNameAndLocation,dtime);
save(savefile); 
%%  Run MCMC
L.N = 128;   % # Eval points for MCMC
Nsamp   = 512; 

num_lumps_to_move = 1;
%sample_mtx = zeros(L.N,L.N,Nsamp); 
L_array    = cell(Nsamp,1); 
gbar_mtx   = zeros(Mx,My,Nsamp);

%Ncenter = 1000; % Should be much larger in practice
jump_sigma   = 4*L.cov(1,1);
jump_sigma_b = 3e-2; 

%L.centers = xmin + (1-2*xmin)*rand(Ncenter,2);

% Plot the initial texture
plot(L); 

f1 = figure; 
f2 = figure;
% Run markov chain 

Knew = L.K;
K_samp = zeros(Nsamp,1);
L_old          = LumpyBgnd;
L_old.N        = L.N;
L_old.cov      = L.cov;
L_old.centers  = L.centers;
L_old.b        = L.b; 
L.israndnumlumps = 0;

ratio_poisson = zeros(Nsamp,1);
ratio_prior   = zeros(Nsamp,1);

[~,gbar_mtx(:,:,1),~] = compute_gaussian_image_lumpy(L,blur_kernel);
%sample_mtx(:,:,1) = L.Eval;
L_array{1} = L.Copy; 

i = 1;
i_samp = 1;
while(i<Nsamp)
    i_samp = i_samp+1;  % Running track of number of trials
    %L.b0  = double(rand(Ncenter,1)<off_probability); 
    
    %fprintf('Number of active lumps = %i\n',Knew);
    
    % Move the old lumps
    idx = randperm(L.K,num_lumps_to_move);
    L.centers(idx,:) = L_old.centers(idx,:) + jump_sigma*randn(num_lumps_to_move,2);
    %L.b(idx)         = L_old.b(idx) + jump_sigma_b*randn(num_lumps_to_move,1);
    [~,gbar_prime,~] = compute_gaussian_image_lumpy(L,blur_kernel);
    figure(f1); subplot(2,2,1); 
    uprime = L.Eval; 
    imagesc(uprime,[0,Lmax]);colorbar;
    set(gca,'YDir','normal');axis image;axis off;title(sprintf('Current sample object, i = %i',i_samp));
    subplot(2,2,2); imagesc(gbar_prime,[0,Lmax]);colorbar;
    set(gca,'YDir','normal');axis image; axis off;title('Current sample image');
    subplot(2,2,3); imagesc(u0);colorbar;
    set(gca,'YDir','normal');axis image; axis off; title('True object');
    subplot(2,2,4); imagesc(g); 
    set(gca,'YDir','normal');axis image; axis off; title('Original image');
    drawnow

    ratio_poisson(i) = poiss_jump_ratio(g,gbar_mtx(:,:,i),gbar_prime);
    ratio_prior(i)   = prior_jump_ratio(L_old,L);
    
    jump_prob = min(1,ratio_poisson(i)*ratio_prior(i));
    fprintf('(idx,ratio_poisson(i),ratio_prior(i),jump_prob) = (%i,%f,%f,%f)\n',idx,ratio_poisson(i),ratio_prior(i),jump_prob);
    eta = rand();
    if(eta<jump_prob)
        fprintf('keeping, i = %i, eta = %f, rho = %f\n',i,eta,jump_prob)
        % store the selected centers
        L_old.centers = L.centers;
        L_array{i+1} = L.Copy; 
        %sample_mtx(:,:,i+1) = L.Eval;
        gbar_mtx(:,:,i+1) = gbar_prime;
        i = i+1;
        figure(f2);
        plot(L);axis off;title(sprintf('Sample num. %i',i));hold on; 
        plot(L.centers(idx,1),L.centers(idx,2),'r*'); 
        drawnow;
    else
        % Restore the old centers to try again
        L.centers = L_old.centers;
        L.b       = L_old.b;
    end
    
end

%%

post_mean = mean(sample_mtx,3); 
post_std  = std(sample_mtx,0,3);

figure; 
subplot(1,3,1); plot(L); title('Original object'); 
subplot(1,3,2); imagesc(post_mean); colorbar; title('Posterior mean');
set(gca,'YDir','normal'); axis image; axis off; 
subplot(1,3,3); imagesc(post_std); colorbar; title('Posterior std'); 
set(gca,'YDir','normal'); axis image; axis off; 
