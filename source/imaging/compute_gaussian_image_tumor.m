function [g,gbar,h_mat] = compute_gaussian_image_tumor(n,h,t,Q,idx,idy)

% Assume that n is an RDESolutionPath
% h is a function handle (the blur kernel/PSF)
% t is a time index vector 
% idx and idy are the measurement indices in ea. direction (x,y)

    n_arr    = Q*n.cell_density(:,:,t);  
    mt = length(t); 
    %idx = 1:8:512;
    %idy = 1:8:512; 
    mx = length(idx); 
    my = length(idy); 

    xp = n.grid{1};
    dx = xp(2); 
    [xx,yy] = meshgrid(xp);
    x0 = 0.5;
    y0 = 0.5;
    h_mat  = h(xx,yy,x0,y0);
    %h_idx = 64:128;
    %h_idy = 64:128;
    h_idx = 193:320; 
    h_idy = 193:320; 
    
    h_mat  = h_mat(h_idx,h_idy);  % Makes the convolution faster.

    gbar = zeros(mx,my,mt);
    fprintf('Simulating image(s)...'); 
    for it=1:mt
        gbartemp = dx^2*conv2(n_arr(:,:,it),h_mat,'same');
        gbar(:,:,it) = gbartemp(idx,idy); 
    end
    g = poissrnd(gbar);   
    fprintf('done!\n'); 
end