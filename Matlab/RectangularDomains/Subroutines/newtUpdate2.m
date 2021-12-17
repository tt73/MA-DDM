function [deltaU] = newtUpdate2(NMatSDD,CMatSDD,weight,uNewt,F,epsilon)
% Compute an update step of Newton's method for the quadrature scheme

% INPUTS
% NMatSDD, CMatSDD are matrices of neighbors and coefficients 
% weight is vector of quadrature weights
% uNewt is the current estimate of the grid function
% F is the source function in the interior
% epsilon is the regularization factor

% OUTPUT
% newtUpdate the a colum

vCount = length(weight);
% Number of points in quadrature

intLength = length(F); Interior = 1:intLength;
pLength = length(uNewt); 

spmat = spdiags((SDDMat(NMatSDD,CMatSDD,uNewt,vCount,epsilon).^-1*weight),0,intLength,intLength)^3;
% correction = 2*pi^2*spdiags((SDDMat(NMatSDD,CMatSDD,uNewt,vCount,epsilon).^-1*weight),0,intLength,intLength)^-3;
% correction is a term that shows up in the jacobian. Each row is to be
% scaled by a specific value, we write as a diagon matrix for faster
% multiplication

% jacobian = sparse(intLength,intLength);
% jacobian = sparse(i,j,vals)... make i,j,vals 
iJac = [];
jJac = [];
vJac = [];

reg = sparse(intLength,intLength);
% reg will be the jacobian of the regularization term

Mk = inf*ones(intLength,1);
% Mk is to keep track of the smallest direction for each point in u, which
% corresponds to the smallest eigenvalue at each grid point


for i = 1:vCount
    % Loop over each search direction, found that 'vectorizing' was too
    % memory intensive, but the loops do start to create a bottle neck as
    % the depth increases
    
    Dvv = sparse(repmat(Interior,1,3),[NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3) ],[CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)],intLength,pLength);
    % Dvv is the SDD difference matrix to the corresponding to current
    % search direction
    
    UvvTemp= full(Dvv*uNewt); Uvv = max(UvvTemp,epsilon);
    % Uvv is the set of regularized SDDs for the given direction
    
    D = weight(i)*spdiags(Uvv.^-2,0,intLength,intLength);
    dTemp = D*Dvv; dTemp = dTemp(Interior,Interior);
    % D,dTemp are just terms that follow out of the jacobian formula
    % Only care about jacobian at the interior points
    
    % add rows dTemp which satisfy this condition
    dTemp(UvvTemp <= epsilon,:) = 0;
    
    [I, J, V] = find(dTemp);
	iJac = cat(1,iJac,I); 
    jJac = cat(1,jJac,J);
    vJac = cat(1,vJac,V);
    % Only updating at UvvTemp > epsilon corresponds to scaling by the
    % heaviside function we get in the derivative of max(epsilon,uvv)
    
    DvvTemp = Dvv(Interior,Interior);
    [Mk, IndKTemp] = min([Mk UvvTemp],[],2);
    Ind = (UvvTemp <= epsilon)&(IndKTemp == 2);
    % There are 2 mins in regularization, need 2 separate heaviside updates
    
    reg(Ind,:) = DvvTemp(Ind,:);
    % reg updates can be redundant, but don't create any bottlenecks
    
    % Can ignore MatLab warnings about reg,jacobian sparse indexing.
end
% Can potentially make loops faster by reducing terms that need to be
% updated, but not priority at the moment.

jacobian = sparse(iJac,jJac,vJac);
jacobian = 2*pi^2*(spmat\jacobian) + reg;
% jacobian = correction*jacobian+reg;
% Combine all terms of jacobian

resid = sparse(pi^2*((SDDMat(NMatSDD,CMatSDD,uNewt,vCount,epsilon).^(-1))*weight).^(-2)+min(min(SDDMat(NMatSDD,CMatSDD,uNewt,vCount,-Inf),epsilon),[],2)-F);
% Could be neater...

deltaU = jacobian\resid;
% Write newton update as linear solve instead of matrix inverse


% Currently getting warnings about ill conditioned matrices, but the
% corresponding matrix is only used in a matrix product, not any linear
% solves.

% Scheme assumes that the quadrature is the same at each grid point and
% that the stencil is the same at each grid point
end