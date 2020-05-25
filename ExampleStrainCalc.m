% This code demonstrates how to calculation using the function
% nanobubble_straintensor_solve.m. This script calls a coarse AFM
% topography, then applys 2D spline interpolation and a FFT based Gaussian
% filter as preprocessing before passing to the function that will
% calculate the strain. The strain map here will roughly reproduce the map
% shown in Figs. 3 and 4a of the main text. 

imL=csvread('Fig3AFM.txt'); % Load AFM topography

dx0=50/3; % Original pixel spacing of AFM in nm. 
dxf=1; % final desired spacing 

r1=dx0/dxf;


% Create x and y vectors for 
[yd,xd]=size(imL);

x=dx0*(0:xd-1)';
y=dx0*(0:yd-1)';

x=x-median(x);
y=y-median(y);

[X0,Y0]=meshgrid(x,y);

xu=linspace(min(x),max(x),ceil(r1*xd));
yu=linspace(min(y),max(y),ceil(r1*yd));

[Xu,Yu]=meshgrid(xu,yu);

% Spline interpolation 

imL_u=interp2(X0,Y0,imL,Xu,Yu,'spline');

% define the Gaussian Filter Kernal 

sigma = 10; % Note filter peforms poorly for sigma values <1. 

gkern = exp(-0.5*(Xu.^2+Yu.^2)/sigma.^2)*1/2/pi/sigma.^2;

% figure, imagesc(x,y,gkern), axis image 

% Appy the filter via the Fourier Convolution Theorem
imu_f=real(ifftshift(ifft2(fft2(imL_u).*fft2(gkern))));

% Plot orignal and filtered AFM images to check performance 
figure(3), subplot(121), imagesc(imL_u), axis image, colormap hot; colorbar
subplot(122), imagesc(imu_f), axis image, colormap hot; colorbar

% Call function nanobubble_straintensor_solve.m to calculate the strain
% tensor 

[eps_xx,eps_yy,eps_xy,chi]=nanobubble_straintensor_solve(imu_f,xu,0.22,80);

