% -------------------------------------------------------------------------
% ************************* RDEVPop Package *******************************
% File:     RDE.m
% Author:   Nick Henscheid
% Date:     10-2019, 2-2020, 5-2020
% Info:     This is the class definition file for the Matlab RDE solver           
% Inputs:         
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------
classdef RDE < handle
    properties 
        rho;   % Growth factor     (Assume is LumpyBgnd class)
        kappa; % Carrying capacity (Assume is LumpyBgnd class)
        u0;    % Initial condition 
        D;     % Diffusion coefficient (Assume is LumpyBgnd class)
        bdy_cond % Boundary condition type 
        grid;  % Grid {x,y} or {x,y,z} (assume rectalinear grid)
    end
    
    properties (Dependent)
        dim;       % Ambient dimension (2 or 3 - inherited from grid)
        n_grid;    % Number of grid points in each direction e.g. [128,128]
        h_grid;    % Mesh width in each direction e.g. [0.01,0.01].  
    end
   
    methods
        function obj = RDE
            % Nothing is set by default - user must define all class
            % properties manually. See example file RDE_Virtual_Pop_2D. 
        end
        
        % Externally defined functions
        n = Solve(obj,T);   % Solve the RDE for the time vector T
        [n,rho,kappa,D] = Sample(obj,N,T,which_to_rand); % Create a virtual population 
                                           % of N subjects via
                                           % randomization
        function Randomize(obj,which_to_rand)
            if(nargin<2)
                which_to_rand = struct('rho',1,'kappa',1,'D',1); 
            end
            % Randomize the coefficients.   Note: this assumes all of these
            % have a method called "randomize".  
            if(which_to_rand.rho)
                obj.rho.Randomize; 
            end
            if(which_to_rand.kappa)
                obj.kappa.Randomize; 
            end
            if(which_to_rand.D)
                obj.D.Randomize;  
            end
        end
        
        % Get methods for dependent properties
        function y = get.dim(obj)
            y = length(obj.grid);  % Ambient dimension set by number of grid vectors
        end
        function y = get.n_grid(obj)
            if(obj.dim==2)
                y = [length(obj.grid{1}),length(obj.grid{2})];
            elseif(obj.dim==3)
                y = [length(obj.grid{1}),length(obj.grid{2}),length(obj.grid{3})];
            end
        end
        function y = get.h_grid(obj)
            % NOTE: this assumes that the grid is uniform. 
            if(obj.dim==2)
                y = [obj.grid{1}(2)-obj.grid{1}(1),...
                     obj.grid{2}(2)-obj.grid{2}(1)];
            elseif(obj.dim==3)
                y = [obj.grid{1}(2)-obj.grid{1}(1),...
                     obj.grid{2}(2)-obj.grid{2}(1),...
                     obj.grid{3}(2)-obj.grid{3}(1)];
            end
        end
    end
    
    methods (Static) 
       Y_o = ODEFun2D(Y_i,D,Aa,Kk,hx,hy,nx,ny); 
       Y_0 = ODEFun3D(Y_i,D,Aa,Kk,hx,hy,hz,nx,ny,nz); 
    end
end