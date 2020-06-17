function u = Eval(obj,XSize)
    % 
    N = obj.N;
    if(numel(N)==1)
        N = N*ones(1,obj.dim);
    end
    
    if(nargin<1) % No array provided
        N = XSize;
    end
    
    u = obj.c*ones(N);
end


 