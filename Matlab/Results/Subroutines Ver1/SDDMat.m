function [Mat] = SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon)
% Compute matrix of all SDDs

% INPUTS
% NMatSDD, CMatSDD are matrices of neighbors and coefficients 
% u is the current estimate of the grid function
% vCount is number of quadrature points
% epsilon is the regularization factor

% OUTPUT
% Mat if the matrix of SDDs, each column corresponds to a direct and each
% row correspons to an interior grid point

Q = CMatSDD.*u(NMatSDD);
% Q is the coefs times u at it's neighbors

Q1 = reshape(Q',3,[]);
% Q1 has columns that are sets of neighbors * coefs

Q2 = sum(Q1);
% Q2 sums the cols of Q1 giving the directional derivatives at each point

Q3 = max(Q2,epsilon);
% Q3 replaces any SDD less than epsilon with epsilon

Mat = reshape(Q3,vCount,[])';
% Mat gives columns of regularized SDDs

% There may be a more efficient way to do this computation
end