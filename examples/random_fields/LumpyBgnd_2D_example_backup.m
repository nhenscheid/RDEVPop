% File: LumpyBgnd_2D_example.m
% Author: Nick Henscheid 
% Date: 4-2019, 5-2020
% Purpose: Demonstrates the lumpy background object functionality (d = 2)

% Initialize a default lumpy background object & display its properties
close all;
L = LumpyBgnd('S',RectSupport([0,1;0,1]));  
L.TurnOffWarnings;  % Don't show warning messages (comment out if you want to see them).
L.Kbar = 200;
L.N = 256;

% Display a realization 
plot(L);title('Realization of LumpyBgnd');
%%
% Randomize the field and re-plot
L.Randomize;
fig2 = figure;
subplot(2,4,1);
plot(L);title('L.cov = 0.005','FontSize',12);

% Change the lump width, but keep it isotropic and keep the same centers
% Note the warning displayed in the command window! The original lump
% centers were selected based on the original lump width 

L.cov = 0.001;
subplot(2,4,2);
plot(L);title('L.cov = 1e-3','FontSize',12);

L.cov = 0.0001;
subplot(2,4,3);
plot(L);title('L.cov = 1e-4','FontSize',12);

L.cov = 0.1;
subplot(2,4,4);
plot(L);title('L.cov = 0.1','FontSize',12);
%%
% Change the lump to be anisotropic
L.cov = diag([0.001,0.01]);
subplot(2,4,5);
plot(L);title('L.cov = diag([1e-3,1e-2])','FontSize',12);

%%
L.cov = diag([0.05,0.0005]);
subplot(2,4,6);
plot(L);title('L.cov = diag([5e-2,5e-4])','FontSize',12);
%%
% Rotate lumps by a fixed rotation matrix 
theta = pi/4;  % Rotation angle
var1  = 0.01;   % var (rotated x-hat)
var2  = 0.001; % var (rotated y-hat)

R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
K = R*diag([var1,var2])*R';

L.cov = K;

subplot(2,4,7);
plot(L);title(sprintf('L.cov = [%1.1e,%1.1e;%1.1e,%1.1e]',L.cov(1,1),L.cov(1,2),L.cov(2,1),L.cov(2,2)),'FontSize',12);


%% Make the lumps non-uniform in size (randomly generated isotropic covariance)
lump_mean = log(0.07);
lump_std  = log(0.12);
L.cov = exp(lump_mean+lump_std*randn())*eye(2);  % Lognormally distributed
L.b   = 1/sqrt((2*pi)^2*det(L.cov));  %PDF/Gaussian mixture model scaling
for i=1:L.K-1
    covtemp = exp(lump_mean+lump_std*randn())*eye(2);
    L.cov = [L.cov,covtemp];
    L.b  = [L.b,1/sqrt((2*pi)^2*det(covtemp))];
end

subplot(2,4,8);
plot(L);title('Indepdent cov for ea. lump','FontSize',12);


