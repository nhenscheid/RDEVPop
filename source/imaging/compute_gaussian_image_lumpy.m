function [g,gbar,h_mat] = compute_gaussian_image_lumpy(f,h,idx,idy,Q)

% Assume that f is a lumpy background
% h is a function handle 

    N_fine   = 512;
    N_eval   = f.N; 
    f.N = N_fine;
    % These are hard-coded for now
    %idx = 1:8:512;
    %idy = 1:8:512; 
    %idx = 64:x_step:442;
    %idy = 64:y_step:442;
    xp = linspace(0,1,N_fine);
    dx = xp(2); 
    [xx,yy] = meshgrid(xp);
    x0 = 0.5;
    y0 = 0.5;
    h_mat  = h(xx,yy,x0,y0);
    h_idx = 193:320;
    h_idy = 193:320;
    h_mat  = h_mat(h_idx,h_idy);  % Makes the convolution faster.

    u0 = Q*f.Eval; 
    gbar = dx^2*conv2(u0,h_mat,'same');
    gbar = gbar(idx,idy);

    g = poissrnd(gbar);
    f.N = N_eval;  % Change back 
end