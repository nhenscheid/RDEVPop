function u = Eval(obj,X,XSize)
    % This function will evaluate the texture for the sample points
    % defined in the array X.
    % This should be updated once the GridData object is ready to
    % go!
    L = obj.L;   % Bounding box
    N = obj.N;
    if(numel(N)==1)
        N = N*ones(1,obj.dim);
    end
    if(nargin<2) % No array provided
        if(obj.dim==2)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            [xx,yy] = meshgrid(xtemp,ytemp);
            X = [xx(:),yy(:)];
            XSize = N;
        elseif(obj.dim==3)
            xtemp = linspace(L(1,1),L(1,2),N(1));
            ytemp = linspace(L(2,1),L(2,2),N(2));
            ztemp = linspace(L(3,1),L(3,2),N(3));
            [xx,yy,zz] = meshgrid(xtemp,ytemp,ztemp);
            X = [xx(:),yy(:),zz(:)];
            XSize = N;
        end
    end
    X = X';
    nEval = size(X,2);
    c_vec = obj.centers';
    if(numel(obj.b)==1)
        b_vec = obj.b*ones(1,obj.K);
    elseif(numel(obj.b)==obj.K)
        b_vec = reshape(obj.b,[1,obj.K]);
    else
        error('Incorrectly formatted b! Must be 1-by-K or K-by-1');
    end
    
    % Process cov into the vectorized form 
    
    if(all(size(obj.cov)==[obj.dim,obj.dim]))
        % Uniform covariance
        temp = inv(obj.cov);
        mask = tril(true(size(obj.cov)));
        temp = temp(mask); 
        cov_vec = temp;
        unif = 1;
        nParam = 1+obj.dim*(obj.dim+1)/2;
    elseif(size(obj.cov,2)==obj.dim*obj.K)
        % Non-uniform covariance
        nParam = 1+obj.dim*(obj.dim+1)/2;
        cov_vec = zeros(nParam-1,obj.K);
        mask = tril(true([obj.dim,obj.dim]));
        unif = 0;
        for i=1:obj.K
            temp = inv(obj.cov(:,(obj.dim*(i-1)+1):(obj.dim*i)));
            temp = temp(mask);
            cov_vec(:,i) = temp;
        end
    else
        error('Incorrectly formatted cov! Must be dim-by-dim or dim-by-K*dim.  Size(cov) = (%i,%i)\n dim = %i\n K = %i',size(obj.cov,1),size(obj.cov,2),obj.dim,obj.K);
    end

    if(obj.gpu)
        if(all(size(cov_vec)==[obj.dim*(obj.dim+1)/2,1]))
            cov_vec = repmat(cov_vec,[1,obj.K]);
        end
        theta_vec = [b_vec;cov_vec];
        u = lumpy_mex_gpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),uint32(nParam),single(c_vec(:)),single(X(:)),single(obj.B0),single(theta_vec(:)));
    else
        u = lumpy_mex_cpu(uint32(obj.dim),uint32(obj.K),uint32(nEval),uint32(unif),c_vec(:),X(:),obj.B0,b_vec(:),cov_vec(:));
    end
    if(any(isnan(u)))
        warning('nan detected!') 
        obj.centers
        obj.K
        b_vec
    end
    if(any(isinf(u)))
        warning('Inf detected!') 
        obj.centers
        obj.K
        b_vec
    end
    if(nargin<2)  % No array provided, reshaping to locally generated grid size
        u = reshape(u,XSize);
    end
    if(nargin==3) % User supplied a grid size, reshape to it
        u = reshape(u,XSize);
    end
end
