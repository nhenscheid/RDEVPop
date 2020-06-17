function r_prior = prior_jump_ratio(L,Lprime)

    Kbar = L.Kbar;
    Kn   = L.K;
    Kprime = Lprime.K;

    r_prior = Kbar^(Kprime-Kn);

end