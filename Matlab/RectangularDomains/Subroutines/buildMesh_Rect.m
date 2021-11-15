function [Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth)
% Builds the discretized domain and the Neighbor,Coefficient Matrices

% INPUTS
% x0,x1,y0,y1 denote the dimensions of a rectangle
% h is the spatial resolution
% depth determines the size of the stencil

% OUTPUTS
% Points is a point cloud of all the points in the domain, each row is a pt
% Interior, Boundary are the Linear indices of the Interior, Boundary
% NMatSDD is the neighbor matrix (of linear indices) where each row 
% corresponds to a linear index of an interior point and columns are in 
% groups of three to represent the neighbors of a
% SDD; in each group the first column is the center (which should be the
% row index) the second column is the index of the forward direction and
% the third column is the index of the backward direction. The column
% groups are arranged in increasing angle of direction from 0 to pi. 
% CMatSDD is the matrix of FD coeffecients and has the same structure as
% NMatSDD, and they can be used by cmat(~)*f(nmat(~)).
% theta is the angles of the search direction in increasing order


bdryTol = h^2;
% Min distance allowed for interior points from the boundary

[direction, dCount, thetaTemp] = cartStencil(depth);
% cartDirection generates all the directions corresponding to a given
% depth. dCount is the number of distinct directions (counting both
% forward and backward directions), directions is a vCountTemp x 2 array where
% the first column is the x direction and the second column is the y
% direction. They are organized in increasing angles from 0 to 2pi. 
% thetaTemp is the angle in radians of the corresponding direction

[D,xgrid,ygrid,Bdry] = sgndist_Rect(x0,x1,y0,y1,h,depth+1);
% sgndistRect takes in the parameters of a rectangle, a spatial
% resolution, and an extension length (depth+1). It outputs an extended
% rectangular grid (xgrid and ygrid), the signed distance field (u) and the
% indices of the boundary. 
% Bdry is a point cloud of subindices of points on the boundary

xIndExtn = repmat((1:length(xgrid))',length(ygrid),1);
yIndExtn = sort(repmat((1:length(ygrid))',length(xgrid),1));
% The subindex vectors for the grid, going from left to right then bottom
% to top

% kIndExtn = xIndExtn+(yIndExtn-1)*length(xgrid);
% the linear index vector for the whole grid

PointsExtn = [repmat(xgrid',length(ygrid),1) sort(repmat(ygrid',length(xgrid),1))];
% Creating a pointcloud from the grids, the index of each point correspond
% to the linear index in kIndExtn

kBdryIndOld= Bdry(:,1)+(Bdry(:,2)-1)*length(xgrid);
% A vector of linear indices for the boundary points

dist = D(sub2ind(size(D),xIndExtn,yIndExtn));
% Mapping the distance field array to a vector, again the index here
% corresponds to the linear index in kIndExtn

IntLog = dist > bdryTol;
% Points less than the boundary tolerance will form the interior

kIntIndOld = find(IntLog);
% A vector of linear indices of the interior points, using the initial
% (old) indices from kIndExtn

Interior = 1:length(kIntIndOld);
% We have the interior points at the start of our main point cloud, so we
% reindex them to start at 1 (new indices)

ExactBdry = ((length(kIntIndOld)+1):(length(kIntIndOld)+length(Bdry)))';
% New linear index of points exactly on the boundary (don't want to make
% duplicates of these points)

NMatTemp = repmat(kIntIndOld,1,dCount);
NMatTemp = NMatTemp + repmat((direction(:,1)+direction(:,2)*length(xgrid))',length(kIntIndOld),1);
% Create a Neighbor Matrix with the old linear indices, relation with
% neighboring points in the linear index if found using the direction
% vector. 
% Only interior points have rows 

old2new = zeros(1,length(PointsExtn));
old2new(kIntIndOld) = Interior;
old2new(kBdryIndOld) = ExactBdry;
% old2new is a vector that allows us to plug in old linear indices and get
% the correspond new ones. For now, we want points not in the interior or
% exactly on the boundary to map to 0 so we can easily find them

NMat = old2new(NMatTemp);
% create a new Neighbor matrix which has 0's for points that need to be
% projected

Points = PointsExtn([kIntIndOld; kBdryIndOld],:);
% Create new point cloud where interior and exact boundary points are now 
% at the proper (new) linear index

bdry_old = setdiff(find(dist<=bdryTol),kBdryIndOld);
% old linear index of points that will projected onto the boundary

bdryIndOld = find(ismember(NMatTemp,bdry_old)==1);
% NMat linear indices of points that will be on the boundary

BoundaryProj = (length(Points)+1:length(Points)+length(bdryIndOld))';
% the new linear indices of the new boundary points

NMat(bdryIndOld) = BoundaryProj;
% Assigning the new index to the boundary points

fixBdryInd = NMatTemp(bdryIndOld);
% old linear indices of boundary points

fixBdryPts = PointsExtn(fixBdryInd,:);
% old coordinate points for new boundary points

[~,fixIntIndTemp,~] = intersect(NMat,BoundaryProj);
[fixIntInd,~] = ind2sub(size(NMat),fixIntIndTemp);
% new linear indices of interior points whose neighbor will be a new boundary
% point

fixIntPts = Points(fixIntInd,:);
% coordinates of interior points whose neighbor will be a new boundary
% point

F = scatteredInterpolant(PointsExtn(:,1),PointsExtn(:,2),dist(:),'linear');
% Interpolated signed distance function that can take in points within the
% comuptational domain

damp = .8;
% damping coef for projection to not overshoot from interpolation error

nMax = 100;
nStep = 0;
% max number of projections
tol = 10^-15; 
% Can probably make this like h^2 but need to check 

bp_new = fixBdryPts;
% initiate projection
IndEr = abs(F(bp_new(:,1),bp_new(:,2)))> tol;
ErCount = sum(IndEr);
while (ErCount > 0) && (nStep < nMax)
    nStep = nStep + 1;
    bp_old = bp_new;
    % assign previous projected value
    
    bp_new(IndEr,:) = bp_old(IndEr,:)+damp*F(bp_old(IndEr,1),bp_old(IndEr,2)).*(bp_old(IndEr,:)-fixIntPts(IndEr,:))./vecnorm((bp_old(IndEr,:)-fixIntPts(IndEr,:)),2,2);
    % damped projection towards the boundary along the corresponding 
    % search direction 
    IndEr = abs(F(bp_new(:,1),bp_new(:,2)))> tol;
    ErCount = sum(IndEr);
end


Points(BoundaryProj,:) = bp_new;
% Assign new coordinate points to new boundary points

Boundary = [ExactBdry;BoundaryProj];

% Since the integral has period pi (rather than 2pi), I'm going to only
% consider the SDDs in directions from [0,pi]

vCount = dCount/2+1;
% vCount is the number quadrature points

NMatSDD = zeros(length(Interior),(vCount)*3);
NMatSDD(:,1:3:(vCount)*3) = repmat(Interior',1,length(1:3:(vCount)*3));
NMatSDD(:,2:3:(vCount)*3) = NMat(:,1:vCount);
NMatSDD(:,3:3:(vCount)*3) = NMat(:,[vCount:dCount 1]);
% The neighbor index matrix for SDDs. Each direction is grouped into three
% neighbors. In each group, the first column is the center, the second
% column is the forward direction neighbor, and the third column is the
% backward direction neighbor. 

X1 = Points(NMatSDD(:,1:3:(vCount)*3),1)-Points(NMatSDD(:,2:3:(vCount)*3),1);
Y1 = Points(NMatSDD(:,1:3:(vCount)*3),2)-Points(NMatSDD(:,2:3:(vCount)*3),2);
H1 = sqrt(X1.^2+Y1.^2);
% X1,Y1 is coordinate differences between the center and the forward
% directions, H1 is the actual distance

X2 = Points(NMatSDD(:,1:3:(vCount)*3),1)-Points(NMatSDD(:,3:3:(vCount)*3),1);
Y2 = Points(NMatSDD(:,1:3:(vCount)*3),2)-Points(NMatSDD(:,3:3:(vCount)*3),2);
H2 = sqrt(X2.^2+Y2.^2);
% X2,Y2 is coordinate differences between the center and the backward
% directions, H2 is the actual distance

A = 2./(H1.^2+H1.*H2);
B = 2./(H2.^2+H1.*H2);
C = -(A+B);
% A, B, and C are the standard finite difference coeffecients 

CMatSDD = zeros(size(NMatSDD));
CMatSDD(:,1:3:(vCount)*3) = reshape(C,size(CMatSDD(:,1:3:(vCount)*3)));
CMatSDD(:,2:3:(vCount)*3) = reshape(A,size(CMatSDD(:,2:3:(vCount)*3)));
CMatSDD(:,3:3:(vCount)*3) = reshape(B,size(CMatSDD(:,3:3:(vCount)*3)));
% CMatSDD is the coeffecient matrix that has same structure as NMatSDD. So
% if we want to use the coef we have CMatSDD(~)*f(NMatSDD(~))
% The reshape just makes assigning the values easier

theta = thetaTemp(1:vCount);
% direction angles from 0 to pi

% Builds a mesh with uniform spatial resolution in both directions, can
% make a more general mesh generator at some point
end

