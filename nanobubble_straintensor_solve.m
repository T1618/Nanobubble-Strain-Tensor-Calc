function [eps_xx,eps_yy,eps_xy,chi]=nanobubble_straintensor_solve(fAFM,x,pr,N)
% solve of the strain of a nanobubble assuming a thin plate governed by the
% Foppl-von Karman equations. Using a Chebyshev spectral collation method
% in cartesian coordiants
% Inputs: 
%       fAFM = matrix of AFM height values. Matrix must be square in this
%           code.
%       x = vector of the x,y coordinates for each element of fAFM. Take
%           care that the units of x are same as the hieght in fAFM.
%       pr = Possoin ratio of the material. 
%       N = number of collocation points. More points will give greater
%       accuracy, but note computation load scales as N^2. Collocation
%       points are unevenly spaced and are concentrated toward the edges.
%       We recommend N > 60, which should solve on order a few seconds. If
%       computational load is a concern, we recommend breaking the sample
%       space into finite elements. 
%       
%  Outputs: 
%       eps_xx, eps_yy, eps_xy = components of the strain tensor at each 
%           pixel. 
%       chi = the Airy stress function normalized to the materials Young's
%       modulus. 

xlr=(max(x)-min(x)); % range of values of lateral dimension

%% Gauss-Lobbatto Points and Chebyshev derivative matrix 

[D,gl]=cheb(N); % function written by Nathanial Jewell based on Lloyd Trefen's Spectral Methods in Matlab, reproduced below.
D2=D^2;


% Define 2D Laplacian Operator Using kron.m take compute the tensor product
% of the second derivative matrix with the identiy matrix. This creates the
% Laplcian which is an (N+1)^2 x (N+1)^2 matrix for the whole space. 
Id=eye(N+1);
L=kron(Id,D2)+kron(D2,Id); 

%% Source function (gaussian curvature)

[Xgl,Ygl]=meshgrid(gl,gl);

xx=Xgl(:); yy=Ygl(:); % colon operator here stacks all columns into a 1D vector.


% Numerically differentiate the AFM topography 

[gx,gy]=gradient(fAFM/xlr*2,2/size(fAFM,1)); 
% [Dividing by xlr/2 normalzes the lateral and vertical lengths into
% dimensionless units.]

[gxx,gxy]=gradient(gx,2/size(fAFM,1));
[~,gyy]=gradient(gy,2/size(fAFM,1));

xi=linspace(-1,1,size(gxx,1));

[Xi,Yi]=meshgrid(xi,xi);

% Use 2D interpolation to calculate the curvature at the collocation points
gxxi=interp2(Xi,Yi,gxx,Xgl,Ygl); 
gyyi=interp2(Xi,Yi,gyy,Xgl,Ygl);
gxyi=interp2(Xi,Yi,gxy,Xgl,Ygl);

f=gxxi.*gyyi-(gxyi).^2;
f=-f(:);

%% Boundary Conditions 

b=find(abs(xx)==1 | abs(yy)==1);

L(b,b)=eye(length(b));

f(b)=0;

%% Numerical Solution and Interpolation 

%%% Solve using  '\' operator.  
tic, 
v=L\f;
v(b)=0; % boundary conidition of equation (9) in Supplimentary Information
chi_t=L\v;
toc  

%% Interpolate and differentiate for the strain.
%%% Interpolation points between the collocation points of chi (Airy
%%% stress function). This is needed for the subsequent differntiation to find
%%% the strain tensor elements. 

Ni=size(fAFM,1); 
% Ni determines the sampling of the outputed strain matrices. Here we choose 
% the sampling to be the same as the input AFM, but the value of Ni
% is arbitrary. 
xf=linspace(-1,1,Ni);

[Xf,Yf]=meshgrid(xf,xf);

% reshape the Airy stress function from [1 x (N+1)^2] to [(N+1) x (N+1)]
chi=reshape(chi_t,[N+1,N+1]); 
chi=interp2(Xgl,Ygl,chi,Xf,Yf,'spline');

dx=2/Ni;

[gix,giy]=gradient(chi,dx);
[gixx,giyx]=gradient(gix,dx);
[~,giyy]=gradient(giy,dx);

% Strain components in terms of the Airy Function curvatures. See Landau Ch. 14. 
eps_xx=giyy-pr*gixx; 
eps_yy=gixx-pr*giyy;
eps_xy=2*(1+pr)*giyx;

 %% Plot Results 
figure(1), imagesc(xf*xlr,xf*xlr,chi), axis image 
title('Airy stress function')
figure(2), p1=subplot(121); imagesc((x-median(x)),(x-median(x)),fAFM), axis image, colorbar('Southoutside')
title('nanobubble topography')
p2=subplot(122); imagesc(xf*xlr,xf*xlr,100*(eps_xx+eps_yy)), axis image, colorbar('Southoutside'), % multiply by 100 to convert to %
title('strain tensor trace')
colormap(p1,'hot')
colormap(p2,'jet')

function [D,x] = cheb(N)
% CHEB  Trefethen's Chebyshev spectral collocation scheme over [-1,1].
%
% [D,x] = cheb(N)
%
%  x - Chebyshev collocation points in descending order
%      [column-vector of length N+1]
%
% D - differentiation matrix
%      [square matrix, size N+1]
%


if N==0, D=0; x=1; return, end

x = cos(pi*(0:N)/N)';
c = [2 ones(1,N-1) 2] .* ((-1).^(0:N));  % c = [2 -1 1 -1 ... 1 -1 2]
X = repmat(x,1,N+1);
D_numer = c' * (1./c);
D_denom = (X - X') + eye(N+1);
D = D_numer./D_denom;
D = D - diag(sum(D'));
