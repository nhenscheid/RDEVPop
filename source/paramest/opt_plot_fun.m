function stop = opt_plot_fun(x,optimValues,state,fig,L0)

nx = length(x); 
nk = (nx-2)/2; 
Lplot = LumpyBgnd; 

Lplot.centers = [x(1:nk),x((nk+1):2*nk)]; 
Lplot.b       = x(2*nk+1);
Lplot.cov     = x(end); 

figure(fig); 
subplot(2,2,1); plot(L0); 
subplot(2,2,2); plot(Lplot); drawnow; 

stop = false; 