function [delta_n,delta_N] = local_fin_diff_lumpy(RDESim,T_N,n_times) 

    % Computes the finite difference (M(beta0+beta)-M(beta0))/norm(beta)
    % Assume that beta = (D,rho,kappa) is a lumpy background parameter
    % RDESim 
    
    RDESim2 = RDE; 
    RDESim2.grid = RDESim.grid; 
    RDESim2.D = RDESim.D;
    RDESim2.rho = RDESim.rho; 
    RDESim2.kappa = RDESim.kappa;
    RDESim2.u0 = RDESim.u0;

    disp('Solving at beta0');
    n_beta_0 = RDESim2.solve(T_N); 
    n_beta_density_0 = n_beta_0.cell_density(:,:,n_times); 
    N_beta_0 = n_beta_0.TumorBurden;
    disp('Done');
    disp('Computing sensitivity matrix'); 
    
    P = 3*RDESim2.D.K + 3*RDESim2.rho.K + 3*RDESim2.kappa.K; 
    
    delta_N = zeros(length(T_N),P);
    delta_n = zeros(prod(RDESim2.n_grid)*length(n_times),P);
    
    dx          = 1e-6; 
    dy          = dx; 
    db_D        = 1e-8;     
    D_tilde     = 1e-7;    % Unit factor for diffusion (cm^2/day)
    db_rho      = 1e-6;
    rho_tilde   = 1;       % Cells/day
    db_kappa    = 1e-4;
    kappa_tilde = 1e8;       % Cells
    % D x pos
    p = 0; 
    disp('Computing sensitivity for D x positions...'); 
    for i=1:RDESim2.D.K
        disp(i);
        p = p+1; 
        % Move this lump center...
        RDESim2.D.centers(i,1) = RDESim2.D.centers(i,1) + dx;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dx;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dx; 
        % ...then move it back
        RDESim2.D.centers(i,1) = RDESim2.D.centers(i,1) - dx;
    end
    disp('...done!');
    
    % D y pos
    disp('Computing sensitivity for D y positions...'); 
    for i=1:RDESim2.D.K
        disp(i)
        p = p+1; 
        % Move this lump center...
        RDESim2.D.centers(i,2) = RDESim2.D.centers(i,2) + dy;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dy;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dy; 
        % ...then move it back
        RDESim2.D.centers(i,2) = RDESim2.D.centers(i,2) - dy;
    end
    disp('...done!'); 
    
    % D amp
    disp('Computing sensitivity for D amplitudes...'); 
    for i=1:RDESim2.D.K
        disp(i)
        p = p+1; 
        % Adjust lump amplitude...
        RDESim2.D.b0(i) = RDESim2.D.b0(i) + db_D;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = D_tilde*(n_beta_density(:)-n_beta_density_0(:))/db_D;
        delta_N(:,p) = D_tilde*(N_beta(:)-N_beta_0(:))/db_D; 
        % ...adjust it back
        RDESim2.D.b0(i) = RDESim2.D.b0(i) - db_D;
    end
    disp('...done!'); 
    
    % rho x pos
    disp('Computing sensitivity for rho x positions...'); 
    for i=1:RDESim2.rho.K
        disp(i);
        p = p+1; 
        RDESim2.rho.centers(i,1) = RDESim2.rho.centers(i,1) + dx;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dx;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dx; 
        RDESim2.rho.centers(i,1) = RDESim2.rho.centers(i,1) - dx;
    end
    disp('...done!');
    
    % rho y pos
    disp('Computing sensitivity for rho y positions...'); 
    for i=1:RDESim2.rho.K
        disp(i);
        p = p+1; 
        RDESim2.rho.centers(i,2) = RDESim2.rho.centers(i,2) + dy;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dy;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dy; 
        RDESim2.rho.centers(i,1) = RDESim2.rho.centers(i,1) - dy;
    end
    disp('...done!');
    
    % rho amp
    disp('Computing sensitivity for rho amplitudes...'); 
    for i=1:RDESim2.rho.K
        disp(i)
        p = p+1; 
        % Adjust lump amplitude...
        RDESim2.rho.b0(i) = RDESim2.rho.b0(i) + db_rho;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = rho_tilde*(n_beta_density(:)-n_beta_density_0(:))/db_rho;
        delta_N(:,p) = rho_tilde*(N_beta(:)-N_beta_0(:))/db_rho; 
        % ...adjust it back
        RDESim2.rho.b0(i) = RDESim2.rho.b0(i) - db_rho;
    end
    disp('...done!'); 
    
    % kappa x pos
    disp('Computing sensitivity for kappa x positions...'); 
    for i=1:RDESim2.kappa.K
        disp(i);
        p = p+1; 
        RDESim2.kappa.centers(i,1) = RDESim2.kappa.centers(i,1) + dx;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dx;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dx; 
        RDESim2.kappa.centers(i,1) = RDESim2.kappa.centers(i,1) - dx;
    end
    disp('...done!');
    
    % kappa y pos
    disp('Computing sensitivity for kappa y positions...'); 
    for i=1:RDESim2.kappa.K
        disp(i);
        p = p+1; 
        RDESim2.kappa.centers(i,2) = RDESim2.kappa.centers(i,2) + dy;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = (n_beta_density(:)-n_beta_density_0(:))/dy;
        delta_N(:,p) = (N_beta(:)-N_beta_0(:))/dy; 
        RDESim2.kappa.centers(i,1) = RDESim2.kappa.centers(i,1) - dy;
    end
    disp('...done!');
    
    % kappa amp
    disp('Computing sensitivity for kappa amplitudes...'); 
    for i=1:RDESim2.kappa.K
        disp(i)
        p = p+1; 
        % Adjust lump amplitude...
        RDESim2.kappa.b0(i) = RDESim2.kappa.b0(i) + db_kappa;
        n_beta = RDESim2.solve(T_N); 
        n_beta_density = n_beta.cell_density(:,:,n_times);
        N_beta = n_beta.TumorBurden; 
        delta_n(:,p) = kappa_tilde*(n_beta_density(:)-n_beta_density_0(:))/db_kappa;
        delta_N(:,p) = kappa_tilde*(N_beta(:)-N_beta_0(:))/db_kappa; 
        % ...adjust it back
        RDESim2.kappa.b0(i) = RDESim2.kappa.b0(i) - db_kappa;
    end
    disp('...done!'); 
   
   
    