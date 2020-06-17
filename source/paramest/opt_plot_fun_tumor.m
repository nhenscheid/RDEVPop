function stop = opt_plot_fun_tumor(x,optimValues,state,fig,n,time,h)

nx = length(x); 
nk = (nx-2)/2; 
Lplot = LumpyBgnd; 

Lplot.centers = [x(1:nk),x((nk+1):2*nk)]; 
Lplot.b       = x(2*nk+1);
Lplot.cov     = x(end); 
[~,gbarplot,~] = compute_gaussian_image_lumpy(Lplot,h);

figure(fig); 
subplot(2,3,1); imagesc(n.cell_density(:,:,time)); set(gca,'YDir','normal'); axis image; colorbar; title('Original n(x)');
subplot(2,3,2); plot(Lplot); colorbar;  title('Current MLE'); 
subplot(2,3,3); imagesc(gbarplot); set(gca,'YDir','normal'); axis image; colorbar;  title('Current gbar'); 
drawnow; 

stop = false; 