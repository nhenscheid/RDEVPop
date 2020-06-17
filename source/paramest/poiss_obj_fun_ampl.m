function l = poiss_obj_fun_ampl(g,X,L,h)
    % Computes the negative log-likelihood for a poisson 
    fprintf('Computing Poisson log-likelihood function...');
    tic
    % First, convert X to a center array: 
    L.centers = [X(1:L.K),X((L.K+1):2*L.K)];
    L.b      = X((2*L.K+1):end);
    
    [~,gbar,~] = compute_gaussian_image_lumpy(L,h);

    l = double(sum(gbar(:) - g(:).*log(gbar(:))));
    fprintf('l = %f, time = %f\n',l,toc); 
end

