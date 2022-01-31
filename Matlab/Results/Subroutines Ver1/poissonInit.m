function [uPM] = poissonInit(NMatSDD,CMatSDD,F,uBdry,v1,v2)
% One step of the Poisson method for an initialization of Newton's Method

% INPUTS
% NMatSDD, CMatSDD are matrices of neighbors and coefficients 
% F is the source function in the interior
% uBdry is the Dirichlet data on the bdry
% v1 and v2 are the directions corresponding to x and y axis

rho = sqrt(2*F);
% rho is the "data" term in the poisson equation

intLength = length(F); 
Interior = 1:intLength;
pLength = length(uBdry)+intLength; 
Boundary = intLength+1:pLength;
% Assumption on data is that the boundary is after the interior points

Dvv1 = sparse(repmat(1:intLength,1,3),[NMatSDD(:,v1*3-2) NMatSDD(:,v1*3-1) NMatSDD(:,v1*3)],[CMatSDD(:,v1*3-2) CMatSDD(:,v1*3-1) CMatSDD(:,v1*3)],intLength,pLength);
Dvv2 = sparse(repmat(1:intLength,1,3),[NMatSDD(:,v2*3-2) NMatSDD(:,v2*3-1) NMatSDD(:,v2*3) ],[CMatSDD(:,v2*3-2) CMatSDD(:,v2*3-1) CMatSDD(:,v2*3)],intLength,pLength);
DvvTemp = Dvv1+Dvv2;
% Dvv1 is the difference matrix for x axis and Dvv2 is for the y axis.
% Combine to get the laplace matrix

uExt = zeros(pLength,1);
uExt(Boundary) = uBdry;
RHS = rho- DvvTemp*uExt;
% uExt is the values only on the boundary, rho-DvvTemp*uExt shifts the known
% boundary data to the function data

Dvv = DvvTemp(Interior,Interior);
uPM = Dvv\sparse(RHS);
% we solve the linear system for our initalization 
% Only gives soln in the Interior

% in 2D the method comes from det(H) = 1/2(trace(H)^2-trace(H^2))
% Currently requires that 1 search direction aligns with x-axis and 1
% search direction aligns with y-axis, but this can be changed.
end

