function[u,xgrid,ygrid,Boundary] = sgndist_Rect(x0,x1,y0,y1,h,tolSteps)
% Builds a larger computational domain and computes the signed distance
% function of the rectangle on the computational domain. 
%
% Note: A signed distance function for a given finite domain and a point x
% returns the shortest distance to the boundary of the domain. If the point
% x is outside the domain, then it returns the shortest distance with a
% negative sign.

% INPUTS
% x0,x1,y0,y1 denote the dimensions of a rectangle, expects x0<x1,y0<y1
% h is the spatial resolution
% tolSteps is how many grid points we need outside of the rectangle for the
% computational domain

% OUTPUTS
% u is the signed distance values on the computational domain
% x,ygrid are row vectors containting the x,y (resp.) disctizations
% Boundary is a point cloud index array of the boundary points, each row
% corresponds to a point of the boundary

if (x0>=x1)||(y0>=y1)
    error('Specified rectangle has negative lengths')
end

xgrid = [(x0-tolSteps*h):h:x0-h x0 x0+h:h:x1-h x1 x1+h:h:(x1+tolSteps*h)];
ygrid = [(y0-tolSteps*h):h:y0-h y0 y0+h:h:y1-h y1 y1+h:h:(y1+tolSteps*h)];
% x,ygrid is the x,y (resp.) discretization including tol

xN = length(xgrid);
yN = length(ygrid);
% x,yN is the length of each discretization

x0_ind = tolSteps + 1; x1_ind = xN-tolSteps;
y0_ind = tolSteps + 1; y1_ind = yN-tolSteps;
% Creates the i,j index for the x and y boundaries

Boundary = [(x0_ind:x1_ind)' repmat(y0_ind,x1_ind-x0_ind+1,1);
            repmat(x1_ind,y1_ind-y0_ind,1) (y0_ind+1:y1_ind)';
            (x1_ind-1:-1:x0_ind)' repmat(y1_ind,x1_ind-x0_ind,1);
            repmat(x0_ind,y1_ind-y0_ind-1,1) (y1_ind-1:-1:y0_ind+1)'];
% Boundary lists the i,j indices of the boundary points.
% Built from bottom -> right -> top -> left

xBody = [1:x0_ind-1 x0_ind+1:x1_ind-1 x1_ind+1:xN]; 
yBody = [1:y0_ind-1 y0_ind+1:y1_ind-1 y1_ind+1:yN];

xVert = [1:x0_ind-1 x1_ind+1:xN];
yVert = [1:y0_ind-1 y1_ind+1:yN];

Main = [repmat(xBody',length(yBody),1) sort(repmat(yBody',length(xBody),1));
            repmat(xVert',2,1) sort(repmat([y0_ind y1_ind]',length(xVert),1));
            sort(repmat([x0_ind x1_ind]',length(yVert),1)) repmat(yVert',2,1)];
% x,yBody and x,yVert are specific i,j indices for the
% non-boundary points for x or y. 
% Main lists the i,j indices of the non-boundary points


gridMax = 10*(max(xgrid)-min(xgrid)+max(ygrid)-min(ygrid));
F = repmat(h^2,xN,yN);
% gridMin is an upper bound on the distance to the boundary
% F is the "source" on the computational boundary. 
% F needs to be h^2*1 for signed distance

u = zeros(xN,yN);            
u(sub2ind([xN yN],Main(:,1),Main(:,2))) = gridMax;
% u is the initializion of our xN by yN solution.
% u at the non-boundary points is set to gridMax (arb large). 
% u is 0 on boundary



II = 1:xN; JJ = 1:yN;
u = sweeps(II,JJ,u,F,xN,yN);
% Computes distance function with fast sweeping method

findInt = (xgrid(Main(:,1))>x0)&(xgrid(Main(:,1))<x1)...
            &(ygrid(Main(:,2))>y0)&(ygrid(Main(:,2))<y1);
        xInt = Main(findInt,1);
        yInt = Main(findInt,2);
u(sub2ind([xN yN],xInt,yInt)) = -u(sub2ind([xN yN],xInt,yInt));
% Currently need to use rectangular structure to turn distance function
% into signed distance

u = -u;
% Distance in interior is positive, exterior is negative 

function[x] = pQuad(a,b,c)
if abs(a-b) >= c
    x = min(a,b) + sqrt(c);
else
    x = (a+b+sqrt(2*c-(a-b)^2))/2;
end
end

function[u] = sweeps(II,JJ,u,F,xN,yN)


for i = II
    for j = JJ
        u(i,j) = min(pQuad(uhxmin(u,i,j,xN),uhymin(u,i,j,yN),F(i,j)),u(i,j));
    end
end
for i = flip(II)
    for j = JJ
        u(i,j) = min(pQuad(uhxmin(u,i,j,xN),uhymin(u,i,j,yN),F(i,j)),u(i,j));
    end
end

for i = flip(II)
    for j = flip(JJ)
        u(i,j) = min(pQuad(uhxmin(u,i,j,xN),uhymin(u,i,j,yN),F(i,j)),u(i,j));
    end
end

for i = II
    for j = flip(JJ)
        u(i,j) = min(pQuad(uhxmin(u,i,j,xN),uhymin(u,i,j,yN),F(i,j)),u(i,j));
    end
end     

function[x] = uhxmin(u,i,j,xN)
    if i == 1
        x = u(2,j);
    elseif i == xN
        x = u(xN-1,j);
    else
        x = min(u(i-1,j),u(i+1,j));
    end
end
function[x] = uhymin(u,i,j,yN)
    if j == 1
        x = u(i,2);
    elseif j == yN
        x = u(i,yN-1);
    else
        x = (min(u(i,j-1),u(i,j+1)));
    end
end
end
% See https://doi.org/10.1090/S0025-5718-04-01678-3 for fast sweeping
% algorithm

% Need more general signed distance solver that doesn't rely on structure
end