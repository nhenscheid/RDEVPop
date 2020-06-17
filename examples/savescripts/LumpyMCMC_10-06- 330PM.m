% Demonstration of Markov Chain Monte Carlo reconstruction of a lumpy-type
% random field from noisy image data.
% The imaging system is assume to be modeled by a Gaussian PSF with Poisson
% 
close all;

FileNameAndLocation='LumpyMCMC';
dtime = datestr(now,'dd-mm-HHMMPM'); 
newbackup=sprintf('./savescripts/%s_%s.m',FileNameAndLocation,dtime);
currentfile=sprintf('./%s.m',FileNameAndLocation);
copyfile(currentfile,newbackup);

%%
L0 = LumpyBgnd; 
L0.Kbar = 50;
L0.b = 2; 
%L.cov = diag([9e-3,3e-4]);
L0.cov = 2.5e-3;
L0.SetPadFactor(0);  % Want to keep lumps fully inside
L0.Randomize;
L0.gpu = 0; 

Mx = 64;
My = 64;

image_FWHM  = 0.1;  % cm (0.1=1mm)
image_sigma = image_FWHM/(2*sqrt(2*log(2)));
image_amp   = 1/(2*pi*image_sigma^2);
blur_kernel = @(x,y,x0,y0) image_amp*exp(-(1/(2*image_sigma^2))*((x-x0).^2 + (y-y0).^2));

tic
[g,gbar,h_mat] = compute_gaussian_image_lumpy(L0,blur_kernel);
fprintf('Time to compute image = %f\n',toc); 

%u0 = L0.eval.data_array;0
u0 = L0.Eval;

Lmax = max(u0(:)); 
%
subplot(2,2,1);
plot(L0);colorbar; title('Object');
subplot(2,2,2);
imagesc(h_mat);axis image;set(gca,'YDir','normal');title('Blur Kernel');
subplot(2,2,3);
imagesc(gbar);axis image;colorbar;set(gca,'YDir','normal');title('Mean image');
subplot(2,2,4);
imagesc(g);axis image;colorbar;set(gca,'YDir','normal');title('Noisy image');
%%
L = LumpyBgnd;
L.N = 64;
L.Kbar = L0.Kbar;
L.cov  = L0.cov;                                                          
L.SetPadFactor(0);
L.israndnumlumps = 0; 
L.gpu = 1; 
%xmin = -3*sqrt(L.cov(1,1));
xmin = 0; 

K_recon = 50;
L.centers = rand(K_recon,2);

% Use this to start at a random initial condition
% L.centers = xmin + (1-2*xmin)*rand(L0.K,2);

% Use this to start at the MLE point 

X0 = rand(3*K_recon,1);
objfun = @(X) poiss_obj_fun_ampl(g,X,L,blur_kernel);
maxiter = 10000; 
options = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations',maxiter,'Algorithm','interior-point');
LB = zeros(size(X0));
UB = ones(size(X0)); 
optTime = tic;
Xstar = fmincon(objfun,X0,[],[],[],[],LB,UB,[],options);
fprintf('Time elapsed = %f\n',toc(optTime)); 
%
L.centers = [Xstar(1:K_recon),Xstar((K_recon+1):2*K_recon)]; 
L.b       = Xstar((2*K_recon+1):end);
L.N = 512; % Can increase # eval points for display
figure; 
subplot(2,2,1); plot(L0); title('Original object'); 
subplot(2,2,2); plot(L);  title('Approx. maximum likelihood estimate') 
subplot(2,2,3); 
imagesc(gbar); axis image; colorbar; set(gca,'YDir','normal'); title('Mean image'); 
subplot(2,2,4); 
imagesc(g); axis image; colorbar; set(gca,'YDir','normal'); title('Noisy image'); 

savefile=sprintf('./results/%s_%s.mat',FileNameAndLocation,dtime);
save(savefile); 
%%  Run MCMC
L.N = 64;   % # Eval points for MCMC
Nsamp   = 512; 

num_lumps_to_move = 1;
sample_mtx = zeros(L.N,L.N,Nsamp); 
gbar_mtx   = zeros(Mx,My,Nsamp);

%Ncenter = 1000; % Should be much larger in practice
jump_sigma   = L.cov(1,1);
jump_sigma_b = 5e-3; 

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
sample_mtx(:,:,1) = L.Eval;

i = 1;
i_samp = 1;
while(i<Nsamp)
    i_samp = i_samp+1;  % Running track of number of trials
    %L.b0  = double(rand(Ncenter,1)<off_probability); 
    
    %fprintf('Number of active lumps = %i\n',Knew);
    
    % Move the old lumps
    idx = randperm(L.K,num_lumps_to_move);
    L.centers(idx,:) = L_old.centers(idx,:) + jump_sigma*randn(num_lumps_to_move,2);
    L.b(idx)         = L_old.b(idx) + jump_sigma_b*randn(num_lumps_to_move,1);
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
        sample_mtx(:,:,i+1) = L.Eval;
        gbar_mtx(:,:,i+1) = gbar_prime;
        i = i+1;
        figure(f2);
        plot(L);axis off;title(sprintf('Sample num. %i',i));drawnow;
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
