function [Mat] = SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon)
% Compute matrix of all SDDs - second direction derivative
% For each of the Ni interior points, there are Ns forward angular
% directions. The goal is to compute the SDD of each interior point in each
% angular direction. The result is a Ni by Ns matrix of SDD values. 

% Note: the 4th arg epsilon regularizes the SSD. 
% If we let epsilon = -Inf, then we disable the regularization. 

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
% Q is Ni by 3*Ns 

Q1 = reshape(Q',3,[]);
% Q1 has columns that are sets of neighbors * coefs
% Q1 is 3 by Ni*Ns 

Q2 = sum(Q1); % col sum 
% Q2 sums the cols of Q1 giving the directional derivatives at each point
% Q2 is 1 by Ni*Ns

Q3 = max(Q2,epsilon);
% Q3 replaces any SDD less than epsilon with epsilon
% dimension does not change 

Mat = reshape(Q3,vCount,[])';
% Mat gives columns of regularized SDDs
% Mat is Ni by Ns 

% There may be a more efficient way to do this computation
end