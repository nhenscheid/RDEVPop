function l = poiss_obj_fun(g,X,L,h)
    % Computes the negative log-likelihood for a poisson 

    % First, convert X to a center array: 
    L.centers = [X(1:L.K),X((L.K+1):2*L.K)];

    [~,gbar,~] = compute_gaussian_image(L,h);

    l = sum(gbar(:) - g(:).*log(gbar(:)));
end

