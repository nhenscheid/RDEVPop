% File: LumpyBgnd_2D_stresstest.m
% Author: Nick Henscheid 
% Date: 6-2020
% Purpose: Testing LumpyBgnd performance and trying to sort out crash bug

% Initialize a default lumpy background object & display its properties
close all;

fid = fopen('log.txt','w'); 
%%

L = LumpyBgnd; 

Ns = [256,512,1024,2048,4096,8192]; 

NNs = length(Ns); 
U = cell(NNs,1); 

for i=1:NNs
    fprintf('Computing %i/%i with N = %i...',i,NNs,Ns(i)); 
    L.N = Ns(i); 
    tic
    U{i} = L.Eval; 
    fprintf('...time elapsed = %f\n',toc); 
end
    
%%
L.N = 1024; 
Nsamp = 1024; 

U2 = L.Sample(Nsamp); 



%% 
image_FWHM  = 0.1;  % cm (0.1=1mm)
image_sigma = image_FWHM/(2*sqrt(2*log(2)));
image_amp   = 1/(2*pi*image_sigma^2);
blur_kernel = @(x,y,x0,y0) image_amp*exp(-(1/(2*image_sigma^2))*((x-x0).^2 + (y-y0).^2));
fprintf(fid,'Computing gaussian images from compute_gaussian_image_lumpy'); 
for i=1:Nsamp
    L.Randomize; 
    [~,~,~] = compute_gaussian_image_lumpy(L,blur_kernel);
    fprintf('%i/%i\n',i,Nsamp); 
    fprintf(fid,'%i/%i\n',i,Nsamp); 
end

fclose(fid); 