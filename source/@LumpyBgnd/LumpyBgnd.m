% -------------------------------------------------------------------------
% ************************* RDEVPop Package *******************************
% File:     LumpyBgnd.m
% Author:   Nick Henscheid
% Date:     9-2016 
% Info:     u = LumpyBgnd(varargin). The lumpy background object.
%           The name-value pairs in varargin determine the random field
%           
% Inputs: 
%           'b0' (a scalar)   The lump amplitude l(x) = b0*l0(x)
%           'B0' (a scalar)   The DC offset i.e. f(x) = B0 + sum()
%           'cov' (a scalar or dim-by-1 vector)
%           'Kbar' (a positive scalar) 
%           'gpu',(0 or 1) If gpu = 1 and a gpu is available, the CUDA 
%                          version is used.
%           'N' (an integer or dim-by-1 vector of integers)
%               
% Contact: nph@email.arizona.edu
% This software is in the public domain, furnished "as is", without 
% technical support, and with no warranty, express or implied, as to its 
% usefulness for any purpose.
% -------------------------------------------------------------------------

classdef LumpyBgnd < handle
    properties (SetObservable = true)
        b   % Lump "Amplitude" (assuming constant for classical LBG)
        B0   % DC offset
        cov  % Lump covariance matrix 
        Kbar % Mean number of lumps
        centers % Lump centers
        gpu  % Evaluate using gpu or not
        N    % Default number of eval. pts in ea. direction 
        israndnumlumps = 1; % Randomize lumps?
    end
    
    properties (Dependent)
        L;% Bounding box for evaluation (lump centers can extend slightly beyond to avoid edge effect issues)
        dim; 
        K; 
    end

    properties (SetAccess=private)
        Support    % Support set for evaluation.  Can only set when LumpyBgnd is first initialized!
        padfactor = 3;      % Padding factor (so that boundary effects don't occur)
        showwarnings = 1;
    end
    
    % Standard methods
    methods
        function obj = LumpyBgnd(varargin)
            p = obj.ParseInputs(varargin{:});
            obj.Kbar    = p.Results.Kbar;
            obj.b       = p.Results.b;
            obj.B0      = p.Results.B0;
            obj.Support = p.Results.Support;
            if(numel(p.Results.cov)==1)
                obj.cov = p.Results.cov*eye(obj.Support.dim);
            else
                obj.cov = p.Results.cov;
            end
            obj.centers = p.Results.centers;
            obj.centers;
            
            if(gpuDeviceCount()>0)
                obj.gpu     = p.Results.gpu;
            else
                obj.gpu     = 0;
            end
            obj.N       = p.Results.N;
            if(numel(obj.centers)==0)
                % Need to generate random lump centers
                obj.Randomize();
            end
            addlistener(obj,'Kbar','PostSet',@LumpyBgnd.HandlePropertyEvents);
            addlistener(obj,'cov','PostSet',@LumpyBgnd.HandlePropertyEvents);
        end
        
        % Get methods for dependent properties
        function val = get.L(obj)
            val = obj.Support.L;
        end
        function val = get.dim(obj)
            val = obj.Support.dim;
        end
        function val = get.K(obj)
            val = size(obj.centers,1);
        end
        function set.K(obj,value)
            if(isnumeric(value)&&value>=0&&(mod(value,1)==0))
                obj.K = value;
            else
                error('Invalid K value');
            end
        end
        % Set methods 
        function set.gpu(obj,value)
            if(isnumeric(value)||islogical(value))
                if(value==1||value==true)
                    if(exist('gpuDeviceCount'))
                        if(gpuDeviceCount>0)
                            obj.gpu = value;
                        else
                            warning('No CUDA-capable GPU Found! Setting gpu to 0.')
                            obj.gpu = 0; 
                        end
                    else
                        warning('Parallel computing toolkit not found! Setting gpu to 0.'); 
                        obj.gpu = 0; 
                    end
                elseif(value==0||value==false)
                    obj.gpu = value;
                else
                    error('Invalid value for gpu!');
                end
            else
                error('Invalid value for gpu!');
            end
        end
        
        % Misc utilities
        function TurnOffWarnings(obj)
            obj.showwarnings = 0;
        end
        
        function TurnOnWarnings(obj)
            obj.showwarnings = 1;
        end
        
        function SetPadFactor(obj,x)
            obj.padfactor = x;
        end
        
        % Externally defined functions
        p = ParseInputs(varargin);
        u = Eval(obj,X,XSize);
        obj = Randomize(obj)
        varargout = plot(obj,alpha);
        z = minus(x,y);  % Method to subtract two lumpy backgrounds
        z = plus(x,y);   % Method to add two lumpy backgrounds
        
        function U = Sample(obj,Ns)
            % Generates Ns samples of the lumpy background, returning it in
            % a cell array U (U{i} is a sample for each i)
            U = cell(Ns,1);
            for i=1:Ns
                obj.Randomize;
                U{i} = obj.Eval;
            end
        end
        
        function v = Copy(obj)
           % Makes a copy of the current object 
            v = LumpyBgnd('S',obj.Support); v.TurnOffWarnings;
            v.centers = obj.centers;
            v.b   = obj.b;
            v.cov = obj.cov;
            v.gpu = obj.gpu;
        end  
    end
    
    % Static methods
    methods (Static)
        HandlePropertyEvents(src,evnt);   % Defined externally
    end
    
end