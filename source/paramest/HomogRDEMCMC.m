% Demonstration of Markov Chain Monte Carlo sampling of the posterior for a
% tumor growth model w/ imaging data. 
% Nick Henscheid 12-2019

%%
close all;
plotting = 1; 
%% Set up spatial grid
Ngrid = 256;               % Number of primary grid points in ea. direction
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
D_0 = UnifRandField('N',Ngrid+2,'a',1e-7,'b',5e-7); 
if(plotting)
    figure(thetafig);subplot(2,2,2);
    plot(D_0); colorbar; set(gca,'CLim',[D_0.a,D_0.b]);
    set(gca,'Position',[0.55,0.55,0.4,0.4])
    title('Sample diffusion coefficient $D(x,y)$ (units cm$^2/$day)','FontSize',14);
end
%% Uniform growth function rho(x,y)
rho_0 = UnifRandField('N',Ngrid+2,'a',0.2,'b',0.5);
% Plot a sample from the growth function
if(plotting)
    figure(thetafig); subplot(2,2,3); 
    plot(rho_0); colorbar; set(gca,'CLim',[rho_0.a,rho_0.b]);
    set(gca,'Position',[0.05,0.05,0.4,0.4])
    title('Sample growth rate $\rho(x,y)$ (units cells/day)','FontSize',14); 
end
%% Lumpy carrying cap function
kappa_0 = UnifRandField('N',Ngrid+2,'a',1e7,'b',5e8);
if(plotting)
    figure(thetafig); subplot(2,2,4);
    plot(kappa_0);colorbar;set(gca,'CLim',[kappa_0.a,kappa_0.b]);
    kappaaxis = gca;
    set(gca,'Position',[0.55,0.05,0.4,0.4]);hold on;
    title('Sample carrying capacity function $\kappa(x,y)$','FontSize',14);
    xlabel('$x$'); ylabel('$y$'); zlabel('$\kappa(x,y)$ (Units cells/cm$^2$)','FontSize',14);
end
%% Solve the forward problem 
[~,n_true] = forward_model(D_0,rho_0,kappa_0,n_0,{x,x},T_0);
if(plotting)
    plot(n_true,1:Nt); 
end
%% Compute the image data for two time points 
t0 = 2; % Time index for first sample
t1 = 8; % Time index for second sample 

image_sigma = 0.01;
image_amp   = 1/(2*pi*image_sigma^2);
blur_kernel = @(x,y,x0,y0) image_amp*exp(-(1/(2*image_sigma^2))*((x-x0).^2 + (y-y0).^2));

idx = 64:4:250;
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
%%   Maximum Likelihood Parameter Estimation for growth parameter only 

D_MLE     = D_0; 
rho_MLE   = UnifRandField('N',Ngrid+2,'a',0,'b',1);
kappa_MLE = kappa_0; 

X0 = rand();
objfun = @(X) poiss_obj_fun_RDE(g,D_MLE,X,kappa_MLE,n_0,blur_kernel,idx,idy,x,T_0,[t0,t1]);
options = optimoptions(@fmincon,'Display','iter','MaxFunctionEvaluations',500,'Algorithm','interior-point');
Xstar = fmincon(objfun,X0,[],[],[],[],0,1,[],options);
rho_MLE.c = Xstar; 
%%
[~,nstar] = forward_model(D_MLE,rho_MLE,kappa_MLE,n_0,{x,x},T_0); 
%%
Nstar = nstar.TumorBurden; 
%%
subplot(2,2,1); imagesc(x,x,n_true.cell_density(:,:,t0)); title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t0),Nstar(t0)));
subplot(2,2,2); imagesc(x,x,nstar.cell_density(:,:,t0));  title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t0),Nstar(t0)));
subplot(2,2,3); imagesc(x,x,n_true.cell_density(:,:,t1)); title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t1),Nstar(t1)));
subplot(2,2,4); imagesc(x,x,nstar.cell_density(:,:,t1));  title(sprintf('True cell density (t=%1.1f), N = %1.3e',T_0(t1),Nstar(t1)));

%%  Run MCMC

Nsamp   = 1024; 

num_lumps_to_move = 1;
sample_mtx = zeros(L.N,L.N,Nsamp); 
gbar_mtx   = zeros(Mx,My,Nsamp);

%Ncenter = 1000; % Should be much larger in practice
jump_sigma   = 0.5*L.cov(1,1);
jump_sigma_b = 1e-3; 

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
L_old.b0       = L.b0; 
L.israndnumlumps = 0;

ratio_poisson = zeros(Nsamp,1);
ratio_prior   = zeros(Nsamp,1);

[~,gbar_mtx(:,:,1),~] = compute_gaussian_image(L,blur_kernel);
sample_mtx(:,:,1) = L.eval;

i = 1;
i_samp = 1;
while(i<Nsamp)
    i_samp = i_samp+1;  % Running track of number of trials
    %L.b0  = double(rand(Ncenter,1)<off_probability); 
    
    %fprintf('Number of active lumps = %i\n',Knew);
    
    % Move the old lumps
    idx = randperm(L.K,num_lumps_to_move);
    L.centers(idx,:) = L_old.centers(idx,:) + jump_sigma*randn(num_lumps_to_move,2);
    L.b0(idx)        = L_old.b0(idx) + jump_sigma_b*randn(num_lumps_to_move,1);
    [~,gbar_prime,~] = compute_gaussian_image(L,blur_kernel);
    figure(f1); subplot(2,2,1); 
    uprime = L.eval; 
    imagesc(uprime,[0,Lmax]);colorbar;
    set(gca,'YDir','normal');axis image;axis off;title(sprintf('Current sample object, i = %i',i_samp));
    subplot(2,2,2); imagesc(gbar_prime,[0,Lmax]);colorbar;
    set(gca,'YDir','normal');axis image; axis off;title('Current sample image');
    subplot(2,2,3); imagesc(n_0);colorbar;
    set(gca,'YDir','normal');axis image; axis off; title('True object');
    subplot(2,2,4); imagesc(g); 
    set(gca,'YDir','normal');axis image; axis off; title('Original image');
    drawnow

    ratio_poisson(i) = poiss_jump_ratio(g,gbar_mtx(:,:,i),gbar_prime);
    ratio_prior(i)   = prior_jump_ratio(L_old,L);
    
    jump_prob = min(1,ratio_poisson(i)*ratio_prior(i));
    fprintf('(ratio_poisson(i),ratio_prior(i),jump_prob) = (%f,%f,%f)\n',ratio_poisson(i),ratio_prior(i),jump_prob);
    eta = rand();
    if(eta<jump_prob)
        fprintf('keeping, i = %i, eta = %f, rho = %f\n',i,eta,jump_prob)
        % store the selected centers
        L_old.centers = L.centers;
        sample_mtx(:,:,i+1) = L.eval;
        gbar_mtx(:,:,i+1) = gbar_prime;
        i = i+1;
        figure(f2);
        plot(L);axis off;title(sprintf('Sample num. %i',i));drawnow;
    else
        % Restore the old centers to try again
        L.centers = L_old.centers;
        L.b0      = L_old.b0;
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
