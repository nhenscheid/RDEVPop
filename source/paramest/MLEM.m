% Computes MLEM solution to g = H(theta) + n
% Poisson noise model 
% Linear H 

function theta = MLEM(g,H,theta0,niter,reltol,plotting)
theta = theta0; 
S = sum(H)'; 
    i=0;
    while i<niter
        Hthetak = H*theta; 
        adjk    = H'*(g(:)./Hthetak);
        thetaprime = theta.*adjk./S; 
        e = norm(theta-thetaprime);
        erel = e/norm(theta); 
        theta = thetaprime; 
        L.b = theta; 
        if(plotting==2)
            subplot(2,2,2);imagesc(L.Eval); set(gca,'YDir','normal'); axis image; title(sprintf('MLEM Reconstruction, i = %i/%i, (e,erel) = (%f,%f)',i,niter,e,erel)); colorbar;  
            subplot(2,2,4);imagesc(reshape(Hthetak,[Mx,My])); set(gca,'YDir','normal'); axis image; title('MLE image');colorbar;  drawnow; 
        elseif(plotting==1)
            fprintf('i = %i/%i, (e,erel) = (%e,%e)\n',i,niter,e,erel); 
        end
        if(erel<reltol)
            i = niter; fprintf('Relative tolerance reached, quitting now!\n'); 
        else 
            i = i+1; 
        end
    end
end