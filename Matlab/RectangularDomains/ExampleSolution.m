addpath('Subroutines')

%%% Solution to the Dirichlet Problem for MA on a square

x0 = -1; x1 = 1; 
y0 = -1; y1 = 1; 
N = 2^3; h = (x1-x0)/(N+1);
depth = 1; 
% Parameters needed to generate grid

[Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);
% Build the mesh
% Points - (Np x 2) arrary of node coordinates 
% Interior - (1 x Ni) interior node indicator 
% Boundary - (Nb x 1) boundary node indicator 
% NMatSDD - (Ni x 3*Nt) indeces of stencil neighbors  
% CMatSDD - (Ni x 3*Nt) stencil coef 
% theta - (Nt x 1) angles of stencil directions (in radians) 

order = 1;
epsilon = h^2;
% Parameters needed to solve problem

DirBC = @(x,y) (exp((x.^2+y.^2)/2));
contF = @(x,y) ((1+x.^2+y.^2).*exp(x.^2+y.^2));
% (10) Example with Radially Smooth Solution

% pos = @(x) max(x,0);
% DirBC = @(x,y) (.5*pos(((x-.5).^2+(y-.5).^2).^.5-.2).^2);
% contF = @(x,y) (pos(1-.2./sqrt((x-.5).^2+(y-.5).^2)));
% (11) Example with C1 Solution


% DirBC = @(x,y) (-sqrt(2-(x.^2+y.^2)));
% contF = @(x,y) (2*(2-(x.^2+y.^2)).^-2);
% (12) Example with blowup on Bdry

F = contF(Points(Interior,1),Points(Interior,2));
% source evaluated at the interior grid points

uBdry = DirBC(Points(Boundary,1),Points(Boundary,2));
% boundary values 

weight = quadWeights(theta,order);
% weights from angle discretization

[uSoln, perf] = quadSolver(NMatSDD,CMatSDD,F,uBdry,epsilon,weight,h);
% Solve with newton's method and no given initial guess. 

figure(2)
plot3(Points(:,1),Points(:,2),uSoln,'.')