function r_prior = prior_jump_ratio_rho(rho,rhoprime)

% Uniform prior
    if(rho>0 && rho<1)
        r_prior = 1;
    else
        r_prior = 0; 
    end
end