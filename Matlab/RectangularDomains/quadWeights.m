function [weight] = quadWeights(theta,order)
% Build the weights for a quadrature on the unit circle

% INPUTS
% Theta should be angles in radians of the discretization of the unit
% circle
% Order is the order of the quadrature that we want

% OUTPUTS
% weight is the corresponding quadrature weight, will be column vector

if nargin == 1
    
    order = 1;
    % Default to order 1
    
elseif nargin == 2
    
    if (order ~= 1) && (order ~=2)
        error('Currently only supports order 1 or order 2')
        % only have scheme for trapezoid and simpson
        
    end
    
end

if isrow(theta)
    theta = theta';
end
% Expects theta to be a column

if order == 1

    dtheta = diff(theta);
    weight = 1/2*([0 dtheta'] + [dtheta' 0])';
    % Order 1 is trapezoid rule
    
elseif order == 2

    dtheta = diff(theta);
    M = length(dtheta);
    weight = zeros(M+1,1);
    
    for i = 0:(M/2-1)
        weight([2*i+1, 2*i+2, 2*i+3])=weight([2*i+1, 2*i+2, 2*i+3])+(dtheta(2*i+1)+dtheta(2*i+2))/6*[(2-dtheta(2*i+2)/dtheta(2*i+1));(dtheta(2*i+1)+dtheta(2*i+2))^2/(dtheta(2*i+1)*dtheta(2*i+2));(2-dtheta(2*i+1)/dtheta(2*i+2))];
        % Can make a little neater...
    end
    % Order 2 is Simpsons rule
    
end

if min(weight) < 0
    
    warning('Negative weights lead to non monotone schemes!')
    % We require all the weights to be non-negative in order for the scheme to
    % be monotone.

end

% You can guarantee the weights are positive if the ratio of consecutive
% angle differences is always less than 2. 

% May add higher order at some point
% Would be curious to know if there is a more general form for these 
end

