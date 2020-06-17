classdef Grid < handle
    properties
        N = [128,128];   % 1-by-d vector of grid sizes
        L = [0,1;0,1];   % d-by-2 matrix of grid limits
    end
    
    properties (Dependent)
        dim;
        Dx;    % Mesh size 
    end
    
    methods
        function obj = Grid
            
            
        end

        function val = get.dim(obj)
            val = size(obj.L,1);
        end
        function val = get.Dx(obj)
            val = (obj.L(:,2) - obj.L(:,1))./(obj.N(:)+1);
        end
        function X = meshgrid(obj)
            if(obj.dim == 2)
                dx = obj.Dx(1); dy = obj.Dx(2); nx = obj.N(1); ny = obj.N(2);
                [xx,yy] = meshgrid(dx:dx:nx*dx,dy:dy:ny*dy);
                X = cat(3,xx,yy);
            end
            
        end
    end
end