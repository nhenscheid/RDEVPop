function [N_beta,n_beta] = forward_model(D,rho,kappa,n0,grid,T) 

    RDESim = RDE; 
    RDESim.grid = grid; 
    RDESim.D = D; 
    RDESim.rho = rho; 
    RDESim.kappa = kappa;
    RDESim.u0 = n0;

    fprintf('Solving forward model...'); 
    n_beta = RDESim.solve(T); 
    fprintf('done!\n'); 
    
    N_beta = n_beta.TumorBurden;
    N_beta = N_beta(end);