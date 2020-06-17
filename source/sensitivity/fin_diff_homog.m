function [delta,N_beta_0,N_beta] = fin_diff_homog(RDESim,beta,T) 

    % Computes the finite difference (M(beta0+beta)-M(beta0))/norm(beta)
    % Assume that beta = (D,rho,kappa) is a homogeneous parameter
    % RDESim 

    RDESim2 = RDE; 
    RDESim2.grid = RDESim.grid; 
    RDESim2.D = RDESim.D;
    RDESim2.rho = RDESim.rho; 
    RDESim2.kappa = RDESim.kappa;
    RDESim2.u0 = RDESim.u0;

    disp('Solving at beta0');
    n_beta_0 = RDESim2.solve(T); 
    N_beta_0 = n_beta_0.TumorBurden;
    disp('Done');
    RDESim2.D.c     = RDESim2.D.c + beta(1);
    RDESim2.rho.c   = RDESim2.rho.c + beta(2); 
    RDESim2.kappa.c = RDESim2.kappa.c + beta(3); 
    
    disp('Solving at beta0+beta');
    n_beta = RDESim2.solve(T); 
    N_beta = n_beta.TumorBurden; 
    disp('Done');
    
    delta = (N_beta(end)-N_beta_0(end))/norm(beta); 