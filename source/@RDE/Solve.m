function n = Solve(obj,T)
    % Solves the PDE on the grid for the times specified
    % Returns an object of type RDESolutionPath

    % Set up the mesh grid
    x = obj.grid{1}; y = obj.grid{2};
    [xx,yy] = meshgrid(x,y);
    nx = obj.n_grid(1); ny = obj.n_grid(2); 
    hx = obj.h_grid(1); hy = obj.h_grid(2);
    n_time = length(T);

    % Compute fields on the grid
    rho_grid = obj.rho.Eval([xx(:),yy(:)],[nx,ny]);
    kappa_grid = obj.kappa.Eval([xx(:),yy(:)],[nx,ny]);
    % NOTE: Dx and Dy must be evaluated on a shifted grid.  See
    % documentation. 
    Dx = obj.D.Eval([xx(:)+hx/2,yy(:)],[nx,ny]);  Dx = Dx(1:end-1,2:end-1);
    Dy = obj.D.Eval([xx(:),yy(:)+hy/2],[nx,ny]);  Dy = Dy(2:end-1,1:end-1);
    Dcell{1} = Dx;
    Dcell{2} = Dy;

    % Set up solution

    n = RDESolutionPath;
    n.grid = obj.grid;
    n.times = T;
    n.cell_density = zeros(nx,ny,length(T));

    % Set up the ODE right-hand-side
    options = odeset; 
    odeset(options,'AbsTol',1e-10,'RelTol',1e-6); 

    if(obj.dim==2)
        odefun = @(T,Y)obj.ODEFun2D(Y,Dcell,rho_grid,kappa_grid,hx,hy,nx,ny);
    elseif(obj.dim==3)
        odefun = @(T,Y)obj.ODEFun3D(Y,Dcell,rho_grid,kappa_grid,hx,hy,nx,ny);
    else
        error('Something is wrong; RDE.dim should only be 2 or 3!'); 
    end

    % Solve the ODE system
    [~,Y] = ode23(odefun,T,obj.u0(:),options);  
    if(n_time==2)
            n.cell_density(:,:,1) = reshape(Y(1,:),[nx,ny]);
            n.cell_density(:,:,2) = reshape(Y(end,:),[nx,ny]);
    else
        for i=1:n_time
            n.cell_density(:,:,i) = reshape(Y(i,:),[nx,ny]); 
        end
    end
end