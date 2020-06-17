function rp = poiss_jump_ratio(g,gbar,gbarprime)

% Computes the transition probability ratio 
% rp = pr(g|gbar(thetaprime))/pr(g|gbar(theta))
% this is part of a Metropolis-Hastings jump transition ratio calc
% the other part is the prior ratio pr(thetaprime)/pr(theta)
% it's assumed that gbar and gbarprime are computed elsewhere 

L      = log_poisspdf_numerator(g(:),gbar(:));
Lprime = log_poisspdf_numerator(g(:),gbarprime(:));

% Note the usage of sum-of-log to avoid underflow
rp = exp(sum(Lprime)-sum(L));

end



function L = log_poisspdf_numerator(g,gbar)
% Computes poisson PDF without factorial denominator
L = -gbar + g.*log(gbar);
end