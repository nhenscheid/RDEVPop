function l = poiss_obj_fun_RDE_rho(g,D,rhoval,kappa,n_0,h,idx,idy,x,T_solve,T_eval)
    % Computes the negative log-likelihood for a poisson 
    % rhoval is the proposed value of the growth parameter 
    % T_solve is a vector of times to solve the PDE 
    % T_eval  is an index vector of times at which the images are taken

    % First, convert rhoval to a UnifRandField: 
    rho_k   = UnifRandField('N',D.N,'a',0,'b',1);
    rho_k.c = rhoval;   
    fprintf('Current rho = %1.6e\n',rhoval); 
    % Solve the forward model 
    [~,n_k] = forward_model(D,rho_k,kappa,n_0,{x,x},T_solve);
    % Compute mean image data 
    [~,gbar,~] = compute_gaussian_image_tumor(n_k,h,T_eval,idx,idy);

    l = sum(gbar(:) - g(:).*log(gbar(:)));
end

