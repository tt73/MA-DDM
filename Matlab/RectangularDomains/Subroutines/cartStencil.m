function [directions, dCount, theta] = cartStencil(depth)
% Generate all of the search directions of a given depth in 2D

% INPUTS
% *depth (positive int) is the L1 radius of the stencil

% OUTPUTS
% *directions (dCount x 2) is a marix of grid aligned vectors, each row is a direction
% *dCount (int) is the number of search directions corresponding to a given depth
% *theta (dCount x 1) is a column vector of the angles corresponding to directions

if (mod(depth,1) ~= 0) || (depth < 1)
    error('Depth must be a whole number')
end

x = [depth:-1:-depth, (1-depth):(depth-1)]';
y = [0:depth, (depth-1):-1:-depth, (-depth+1):-1]';
% Generate all pairs of coordinates which solve |x|+|y| = depth (L1 Circle)
% Goes counter clockwise from (depth,0)

directions = [x y];
% directions which have unique angles

dCount = length(directions);
% number of individual directions

theta = atan2(y,x);
% angles in radians

theta(theta <0) = 2*pi + theta(theta<0);
[theta, sortInd] = sort(theta);
% angles of search directions from 0 to 2pi in ascending order 

directions = directions(sortInd,:);
% put directions in order of their angles

% depth = 1 gives 5 point stencil
% depth = 2 gives 9 point stencil
% depth = 3 gives 13 point stencil
% ...
end

