function l = poiss_obj_fun_tumor(g,X,L,h,Q)
    % Computes the negative log-likelihood for a poisson 
    fprintf('Computing Poisson log-likelihood function...');
    tic
    % First, convert X to a center array: 
    L.centers = [X(1:L.K),X((L.K+1):2*L.K)];
    L.b       = X((2*L.K+1):3*L.K);
    L.cov     = X(end); 
    
    [~,gbar,~] = compute_gaussian_image_lumpy(L,h);
    gbar = Q*gbar; 
    
    gbar(gbar==0) = eps; % prevent NaN for 0*log(0)
    
    l = double(sum(gbar(:) - g(:).*log(gbar(:))));
    fprintf('l = %f, time = %f\n',l,toc); 
end

