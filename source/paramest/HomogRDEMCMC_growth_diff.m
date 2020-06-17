% Demonstration of Markov Chain Monte Carlo sampling of the posterior for a
% tumor growth model w/ imaging data. 
% Nick Henscheid 12-2019

%%
close all;
plotting = 1; 
%% Set up spatial grid
Ngrid = 512;               % Number of primary grid points in ea. direction
x = linspace(0,1,Ngrid+2); % NOTE: we need Ngrid+2 points for ghost cells
[xg,yg] = meshgrid(x);     % Uniform grid on [0,1]^2  (Spatial units are cm)
%% Set up time grid
Nt = 10;                    % Number of time points to keep 
tmax = 365;               % Time to simulate (days)
T_0 = linspace(0,tmax,Nt);
%% Set up initial condition
N0 = 5;       % Integral of intial cell density ("number of initial cells")
s  = 0.01;    % Std. dev of initial condition
initcond = @(x,y,xc,yc) N0*exp(-((x-xc).^2+(y-yc).^2)/(2*s^2))/(2*pi*s^2);
n_0 = initcond(xg,yg,0.5,0.5);
if(plotting)
    thetafig = figure; 
    set(thetafig,'Position',[thetafig.Position(1),thetafig.Position(2),1200,800]);
    subplot(2,2,1);
    imagesc(x,x,n_0);axis image; set(gca,'YDir','normal'); colorbar;
    set(gca,'Position',[0.05,0.55,0.4,0.4])
    title(sprintf('Initial Condition ($N_0$ = %f)',N0),'FontSize',14);
end
%% Uniform isotropic diffusion coefficient D(x,y)
D_a = 1e-7; D_b = 5e-7; 
D_dist = @()D_a + (D_b-D_a)*rand(); 
D_0 = UnifRandField('N',Ngrid+2,'dist',D_dist); 
if(plotting)
    figure(thetafig);subplot(2,2,2);
    plot(D_0); colorbar; set(gca,'CLim',[D_a,D_b]);
    set(gca,'Position',[0.55,0.55,0.4,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end
%% Uniform growth function rho(x,y)
rho_a = 0.2; rho_b = 0.5; 
rho_dist = @()rho_a + (rho_b-rho_a)*rand(); 
rho_0 = UnifRandField('N',Ngrid+2,'dist',rho_dist);
% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    plot(rho_0); colorbar; set(gca,'CLim',[rho_a,rho_b]);
    set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Sample growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
end
%% Lumpy carrying cap function
kappa_a = 1e7; kappa_b = 5e8; 
kappa_dist = @()kappa_a + (kappa_b-kappa_a)*rand(); 
kappa_0 = UnifRandField('N',Ngrid+2,'dist',kappa_dist);
if(plotting)
    figure(thetafig); subplot(2,2,4);
    plot(kappa_0);colorbar;set(gca,'CLim',[kappa_a,kappa_b]);
    kappaaxis = gca;
    set(gca,'Position',[0.55,0.05,0.4,0.4]);hold on;
    title('Sample carrying capacity function $\kappa(x,y)$','FontSize',14);
    xlabel('$x$'); ylabel('$y$'); zlabel('$\kappa(x,y)$ (Units cells/cm$^2$)','FontSize',14);
end
%% Solve the forward problem w/ the "true' 
[~,n_true] = forward_model(D_0,rho_0,kappa_0,n_0,{x,x},T_0);

if(plotting)
    plot(n_true,1:Nt); 
end
%% Compute the image data for two time points 
t0 = 2; % Time index for first sample
t1 = 3; % Time index for second sample 

image_sigma = 0.05;
image_amp   = 1/(2*pi*image_sigma^2);
blur_kernel = @(x,y,x0,y0) image_amp*exp(-(1/(2*image_sigma^2))*((x-x0).^2 + (y-y0).^2));

idx = 200:4:400;
idy = idx; 

[g,gbar,h_mat] = compute_gaussian_image_tumor(n_true,blur_kernel,[t0,t1],idx,idy);

if(plotting)
    figure; 
    ax = tight_subplot(2,4,0.08,0.05,0.05);
    axes(ax(1)); 
    imagesc(n_true.grid{1},n_true.grid{2},n_true.cell_density(:,:,t0)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$n(x,%1.1f)$',T_0(t0))); 
    axes(ax(2)); 
    imagesc(n_true.grid{1},n_true.grid{2},gbar(:,:,1)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$\\bar{g}(%1.1f)$',T_0(t0))); 
    axes(ax(3)); 
    imagesc(n_true.grid{1},n_true.grid{2},g(:,:,1)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$g(%1.1f)$',T_0(t0))); 
    axes(ax(4)); 
    imagesc(n_true.grid{1},n_true.grid{2},g(:,:,1)-gbar(:,:,1)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$g(%1.1f)-\\bar{g}(%1.1f)$',T_0(t0),T_0(t0))); 
    axes(ax(5)); 
    imagesc(n_true.grid{1},n_true.grid{2},n_true.cell_density(:,:,t1)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$n(x,%1.1f)$',T_0(t1))); 
    axes(ax(6)); 
    imagesc(n_true.grid{1},n_true.grid{2},gbar(:,:,2)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$\\bar{g}(%1.1f)$',T_0(t1))); 
    axes(ax(7)); 
    imagesc(n_true.grid{1},n_true.grid{2},g(:,:,2)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$g(%1.1f)$',T_0(t1))); 
    axes(ax(8)); 
    imagesc(n_true.grid{1},n_true.grid{2},g(:,:,2)-gbar(:,:,2)); 
    set(gca,'YDir','normal'); 
    title(sprintf('$g(%1.1f)-\\bar{g}(%1.1f)$',T_0(t1),T_0(t1))); 
end
%%   Maximum Likelihood Parameter Estimation *** for both (D,rho) *** 
D_MLE_a = 0; 
D_MLE_b = 1; 
D_MLE_dist = @()D_MLE_a + (D_MLE_b-D_MLE_a)*rand(); 
D_MLE     = UnifRandField('N',Ngrid+2,'dist',D_MLE_dist);
rho_MLE_a = 0; 
rho_MLE_b = 1; 
rho_MLE_dist = @()rho_MLE_a + (rho_MLE_b-rho_MLE_a)*rand(); 
rho_MLE   = UnifRandField('N',Ngrid+2,'dist',rho_MLE_dist);
kappa_MLE = kappa_0; 

X0 = [1e-6*rand(),rand()];
objfun = @(X) poiss_obj_fun_RDE(g,X,kappa_MLE,n_0,blur_kernel,idx,idy,x,T_0,[t0,t1]);
options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',500,'Algorithm','interior-point',...
                                'StepTolerance',1e-4);
Xstar = fmincon(objfun,X0,[],[],[],[],0,1,[],options);

D_MLE.c   = Xstar(1); 
rho_MLE.c = Xstar(2); 

%   Compute resulting forward model at the MLE estimate 
[~,nstar] = forward_model(D_MLE,rho_MLE,kappa_MLE,n_0,{x,x},T_0); 
Nstar = nstar.TumorBurden; 
%% Plot the result 
if(plotting)
    figure; 
    subplot(2,2,1); imagesc(x,x,n_true.cell_density(:,:,t0)); title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t0),Nstar(t0)));
    subplot(2,2,2); imagesc(x,x,nstar.cell_density(:,:,t0));  title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t0),Nstar(t0)));
    subplot(2,2,3); imagesc(x,x,n_true.cell_density(:,:,t1)); title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t1),Nstar(t1)));
    subplot(2,2,4); imagesc(x,x,nstar.cell_density(:,:,t1));  title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t1),Nstar(t1)));
end

%% Estimate the sampling distribution of the MLE by re-sampling g ~ poi(gbar)

N_samp_mle = 64;
rhostar_mle = zeros(N_samp_mle,1); 
for i=1:N_samp_mle
    fprintf('Computing MLE estimate number %i\n',i); 
    tic
    X0 = rand(); 
    g = poissrnd(gbar); 
    objfun = @(X) poiss_obj_fun_RDE(g,D_MLE,X,kappa_MLE,n_0,blur_kernel,idx,idy,x,T_0,[t0,t1]);
    rhostar_mle(i) = fmincon(objfun,X0,[],[],[],[],0,1,[],options);
    fprintf('Time elapsed = %f\n',toc); 
end
fname = sprintf('test_%s.mat', datestr(now,'mm-dd-yyyy HH-MM'));
save(fname);
%%  Run MCMC on a particular sample

N_samp_post  = 2048; 

rho_post   = zeros(N_samp_post,1); % Storing the posterior samples here
D_post     = zeros(N_samp_post,1); % Storing the posterior samples here
samples    = cell(N_samp_post,1);  % Storing the resulting cell densities here
gbar_post  = cell(N_samp_post,1);  % Storing the mean image data here
accepted   = zeros(N_samp_post,1); % 0-1 vector to compute acceptance ratio
accepted(1) = 1; 

jump_sigma_rho = 1e-7;
jump_sigma_D   = 1e-11; 

f1 = figure('Position',[100,100,1000,400]); ax1 = tight_subplot(2,4,0.08,0.08,0.08);
%f2 = figure;

% Run markov chain 

rho_old    = UnifRandField('N',Ngrid+2,'a',0,'b',1);
rho_new    = UnifRandField('N',Ngrid+2,'a',0,'b',1);
rho_old.c  = rho_MLE.c; 
rho_new.c  = rho_MLE.c;
D_old      = UnifRandField('N',Ngrid+2,'a',0,'b',1);
D_new      = UnifRandField('N',Ngrid+2,'a',0,'b',1);
D_old.c    = D_MLE.c; 
D_new.c    = D_MLE.c;

ratio_poisson = zeros(N_samp_post,1);
ratio_prior   = zeros(N_samp_post,1);

[~,gbar_post{1},~] = compute_gaussian_image_tumor(nstar,blur_kernel,[t0,t1],idx,idy);
samples{1}         = nstar;

axes(ax1(1));  imagesc(nstar.cell_density(:,:,t0));
set(gca,'YDir','normal');axis off;title(sprintf('Proposal object (t = %1.2f), i = %i',T_0(t0),i));
axes(ax1(2));  imagesc(gbar_post{1}(:,:,1));
set(gca,'YDir','normal');axis off;title('Proposal image');
axes(ax1(3));  imagesc(nstar.cell_density(:,:,t1));
set(gca,'YDir','normal');axis off; title(sprintf('Proposal object (t = %1.2f)',T_0(t1)));
axes(ax1(4));  imagesc(gbar_post{1}(:,:,2));
set(gca,'YDir','normal');axis off;title('Proposal image');
axes(ax1(5));  imagesc(n_true.cell_density(:,:,t0));
set(gca,'YDir','normal');axis off;title(sprintf('Original object (t = %1.2f)',T_0(t0)));
axes(ax1(6));  imagesc(g(:,:,1));
set(gca,'YDir','normal');axis off;title('Original image');
axes(ax1(7));  imagesc(n_true.cell_density(:,:,t1));
set(gca,'YDir','normal');axis off;title(sprintf('Original object (t = %1.2f)',T_0(t1)));
axes(ax1(8));  imagesc(g(:,:,2));
set(gca,'YDir','normal');axis off;title('Original image');
drawnow

i = 1;
i_samp = 1;
for i=2:N_samp_post
    % Propose a new rho
    rho_new.c = rho_old.c + jump_sigma_rho*randn();
    D_new.c   = D_old.c   + jump_sigma_D*randn(); 
    % Compute forward model for proposed rho
    [~,n_proposal]    = forward_model(D_new,rho_new,kappa_MLE,n_0,{x,x},T_0); 
    % Compute image data for the proposed rho
    [~,gbar_prime,~]  = compute_gaussian_image_tumor(n_proposal,blur_kernel,[t0,t1],idx,idy);
    
    set(ax1(1).Children,'CData',n_proposal.cell_density(:,:,t0));
    set(gca,'YDir','normal');axis image;axis off;
    set(ax1(1).Title,'String',sprintf('$(D,\\rho)=(%1.4e,%1.4e)$, t = %1.2f), i = %i',D_new.c,rho_new.c,T_0(t0),i));
    set(ax1(2).Children,'CData',gbar_prime(:,:,1));
    set(gca,'YDir','normal');axis image; axis off;title('Proposal image');
    set(ax1(3).Children,'CData',n_proposal.cell_density(:,:,t1));
    set(gca,'YDir','normal');axis image; axis off; 
    title(sprintf('Proposal object (t = %1.2f)',T_0(t1)));
    set(ax1(4).Children,'CData',gbar_prime(:,:,2));
    set(gca,'YDir','normal');axis image; axis off;title('Proposal image');
    set(ax1(5).Children,'CData',n_true.cell_density(:,:,t0));
    set(gca,'YDir','normal');axis image; axis off;title(sprintf('Original object (t = %1.2f)',T_0(t0)));
    set(ax1(6).Children,'CData',g(:,:,1));
    set(gca,'YDir','normal');axis image; axis off;title('Original image');
    set(ax1(7).Children,'CData',n_true.cell_density(:,:,t1)); 
    set(gca,'YDir','normal');axis image; axis off;title(sprintf('Original object (t = %1.2f)',T_0(t1)));
    set(ax1(8).Children,'CData',g(:,:,2));
    set(gca,'YDir','normal');axis image; axis off;title('Original image');
    drawnow

    ratio_poisson(i) = poiss_jump_ratio(g,gbar_post{i-1},gbar_prime);
    % NEED TO UPDATE THIS...
    ratio_prior(i)   = prior_jump_ratio_rho_D(rho_new.c,rho_old.c,D_new.c,D_old.c);
    
    jump_prob = min(1,ratio_poisson(i)*ratio_prior(i));
    fprintf('(ratio_poisson(i),ratio_prior(i),jump_prob) = (%f,%f,%f)\n',ratio_poisson(i),ratio_prior(i),jump_prob);
    eta = rand();
    if(eta<jump_prob)
        fprintf('Keeping, i = %i, (D,rho) = (%1.6e,%1.6e), eta = %f, jump_prob = %f\n',i,D_new.c,rho_new.c,eta,jump_prob)
        % store the selected centers
        rho_old.c      = rho_new.c;
        rho_post(i)    = rho_new.c; 
        D_old.c        = D_new.c; 
        D_post(i)      = D_new.c; 
        % New sample added
        samples{i}     = n_proposal;
        gbar_post{i}   = gbar_prime;
        accepted(i)      = 1; 
    else
        fprintf('Rejecting, i = %i, (D,rho) = (%1.6e,%1.6e), eta = %f, jump_prob = %f\n',i,D_new.c,rho_new.c,eta,jump_prob)
        % Restore the old rho value to try again
        rho_new.c    = rho_old.c;
        rho_post(i)  = rho_old.c; 
        D_new.c      = D_old.c; 
        D_post(i)    = D_old.c; 
        % Old sample added
        samples{i}   = samples{i-1}; 
        gbar_post{i} = gbar_post{i-1}; 
    end
end

%%
fname = sprintf('MCMC_test_%s.mat', datestr(now,'mm-dd-yyyy HH-MM'));
save(fname);
%%

post_mean = mean(sample_mtx,3); 
post_std  = std(sample_mtx,0,3);

figure; 
subplot(1,3,1); plot(L); title('Original object'); 
subplot(1,3,2); imagesc(post_mean); colorbar; title('Posterior mean');
set(gca,'YDir','normal'); axis image; axis off; 
subplot(1,3,3); imagesc(post_std); colorbar; title('Posterior std'); 
set(gca,'YDir','normal'); axis image; axis off; 
