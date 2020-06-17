function l = poiss_obj_fun_RDE(g,X,kappa,n_0,h,idx,idy,x,T_solve,T_eval)
    % Computes the negative log-likelihood for a poisson 
    % rhoval is the proposed value of the growth parameter 
    % T_solve is a vector of times to solve the PDE 
    % T_eval  is an index vector of times at which the images are taken

    % First, convert rhoval to a UnifRandField: 
    Dval    = X(1);
    rhoval  = X(2); 
    rho_k   = UnifRandField('N',kappa.N);
    D_k     = UnifRandField('N',kappa.N);
    rho_k.c = rhoval; 
    D_k.c   = Dval; 
    fprintf('Current (D,rho) = (%1.6e,%1.6e)\n',Dval,rhoval); 
    % Solve the forward model 
    [~,n_k] = forward_model(D_k,rho_k,kappa,n_0,{x,x},T_solve);
    % Compute mean image data 
    [~,gbar,~] = compute_gaussian_image_tumor(n_k,h,T_eval,idx,idy);

    l = sum(gbar(:) - g(:).*log(gbar(:)));
end

