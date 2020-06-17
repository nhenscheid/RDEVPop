function r_prior = prior_jump_ratio_rho_D(rho,rhoprime,D,Dprime)

% Uniform prior
    if(rho>0 && rho<1 && D>0 && D<1e-4)
        r_prior = 1;
    else
        r_prior = 0; 
    end
end