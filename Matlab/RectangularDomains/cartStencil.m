function [directions, dCount, theta] = cartStencil(depth)
% Generate all of the search directions of a given depth in 2D

% INPUTS
% depth is the sup norm of the stencil, must be a whole number

% OUTPUTS
% directions is a marix of grid aligned vectors, each row is a direction
% dCount is the number of search directions corresponding to a given depth
% theta is a column vector of the angles corresponding to directions

if (mod(depth,1) ~= 0) || (depth < 1)
    error('Depth must be a whole number')
end

x_temp1 = repmat((-depth:depth)',2*depth+1,1);
y_temp1 = sort(repmat((-depth:depth)',2*depth+1,1));
% Generate all pairs of coordinates where each |coordinate| <= depth

R_temp1 = abs(x_temp1)+abs(y_temp1);
% Calculate L1 norm for all directions

ind2 = (R_temp1 > 0);
% Look for points that aren't the origin

x_temp2 = x_temp1(ind2);
y_temp2 = y_temp1(ind2);
% Set of coordinates that aren't the origin 

R_temp2 = R_temp1(ind2);
% L1 norm of remaining points

theta_temp2 = atan2(y_temp2,x_temp2);
% Angles of all remaining points

thetaTol = .1/(4*(4*depth+1)^2);
theta_temp3 = uniquetol(theta_temp2,thetaTol);
% Set of unique angles, want to get rid of duplicate angles such as (1,1)
% and (2,2). We will drop the one with the larger L1 norm

indTrue = zeros(size(theta_temp2));
% Initiate index of directions we keep

for t = theta_temp3'
    % Loop through unique angles 
    
    t_ind = find(abs(theta_temp2-t)<=thetaTol);
    % find index of all directions of a given angle
    
    R_t = R_temp2(t_ind);
    % check L1 norm of those directions
    
    [~,R_ind] = min(R_t);
    % index of minimum norm in the subset of directions
    
    theta_ind = t_ind(R_ind);
    % find index of smallest norm out of all directions
    
    indTrue(theta_ind) = 1;
    % set index to true
end
indTrue = indTrue == 1;
% convert array to actual logical array

directions = [x_temp2(indTrue) y_temp2(indTrue)];
% directions which have unique angles

dCount = length(directions);
% number of individual directions

theta = theta_temp2(indTrue);
theta(theta <0) = 2*pi + theta(theta<0);
[theta, sortInd] = sort(theta);
% angles of search directions from 0 to 2pi in ascending order 

directions = directions(sortInd,:);
% put directions in order of their angles

% depth = 1 gives 9 point stencil
% depth = 2 gives 17 point stencil
% depth = 3 gives 33 point stencil
% ...
end

