function [uNewt,t,stepcount] = quadSolver3(NMatSDD,CMatSDD,Dvvs,F,uBdry,epsilon,weight,h,uInit)
% Solves Dirichlet problem for Monge-Ampere using the quadrature based
% finite difference scheme via Newton's Method

% INPUTS
% NMatSDD, CMatSDD are matrices of neighbors and coefficients
% F is the source evaluated at the interior
% uBdry is the Dirichlet data at the boundary
% epsilon is the regularization parameter
% weight is the column vector of quadrature weights
% uInit is an optional initial guess for Newton's Method

% OUTPUTS
% uNewt is the result of running Newton's Method
% t is the computational time for Newton's Method
% stepcount is the number if steps it took for convergence

vCount = length(weight);
% number of points in the quadrature

if (~exist('uInit','var'))
    uInit = poissonInit(NMatSDD,CMatSDD,F,uBdry,1,(vCount+1)/2);
    % Does a single step of the poisson iteration to initialize the solver
end

aproxMAOp = @(u)(pi^2*((SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon).^(-1))*weight).^(-2)+min(min(SDDMat(NMatSDD,CMatSDD,u,vCount,-Inf),epsilon),[],2));
% regularized quadrature based finite difference operator (See writeup)

uNewt = [uInit; uBdry];
% Initialize Newton's Method
bZ = zeros(length(uBdry),1);

resid = norm(aproxMAOp(uNewt)-F,inf);
stepcount = 0;
tic
    deltaU = [newtUpdate2(NMatSDD,CMatSDD,Dvvs,weight,uNewt,F,epsilon);bZ];
    % Only updated in the Interior
    
    alpha = 1;
    uNewtTemp = uNewt -alpha*deltaU;
    residTemp = norm(aproxMAOp(uNewtTemp)-F,inf);
    % alpha is our dampening coeffecient, reset it to 1 at each iteration,
    % newtUpdate_modQuad takes in the current u and the values of f at the interior
    % and gives an update for newton's method. uNewtTemp and residTemp are
    % the values of u and the residual with this update
    
    while residTemp > resid
        
        alpha = alpha/2;
        uNewtTemp = uNewt -alpha*deltaU;
        residTemp = norm(aproxMAOp(uNewtTemp)-F,inf);
        % We bisect alpha until the residual decreases
        
        if alpha < 10^-16
            break
        end
        
    end

    
%     resid = residTemp;
    uNewt = uNewtTemp;
    stepcount = stepcount + 1;
t = toc;
% Track the time and number of steps for NM to converge

end
