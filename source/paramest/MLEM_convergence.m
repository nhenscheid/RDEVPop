% Computes MLEM solution to g = H(theta) + n
% Poisson noise model 
% Linear H 
% Compares to a known parameter theta_true to assess convergence

function [theta,err,erel] = MLEM_convergence(g,H,theta0,niter,reltol,plotting,L,n_true)
err = zeros(niter,1); 
erel = zeros(niter,1); 
Ltest = L.Copy; 
Ltest.N = 512; 
L.gpu = 1; 
theta = theta0; 

S = sum(H)'; 
    i=0;
    while i<niter
        Hthetak = H*theta; 
        adjk    = H'*(g(:)./Hthetak);
        thetaprime = theta.*adjk./S; 
        e = norm(theta-thetaprime);
        erel(i+1) = e/norm(theta); 
        Ltest.b = thetaprime; 
        nprime = Ltest.Eval; 
        err(i+1) = norm(nprime(:)-n_true(:))/norm(n_true(:)); 
        theta = thetaprime; 
        L.b = theta; 
        if(plotting==2)
            subplot(2,2,2);imagesc(L.Eval); set(gca,'YDir','normal'); axis image; title(sprintf('MLEM Reconstruction, i = %i/%i, (e,erel) = (%f,%f)',i,niter,e,erel)); colorbar;  
            subplot(2,2,4);imagesc(reshape(Hthetak,[Mx,My])); set(gca,'YDir','normal'); axis image; title('MLE image');colorbar;  drawnow; 
        elseif(plotting==1)
            fprintf('i = %i/%i, (erel,err) = (%e,%e)\n',i,niter,erel(i+1),err(i+1)); 
        end
        if(erel(i+1)<reltol)
            i = niter; fprintf('Relative tolerance reached, quitting now!\n'); 
        else 
            i = i+1; 
        end
    end
end