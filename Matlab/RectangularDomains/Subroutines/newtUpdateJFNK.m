function [deltaU,r] = newtUpdateJFNK(NMatSDD,CMatSDD,Dvvs,weight,uNewt,F,epsilon,delta)
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

% correction = 2*pi^2*spdiags((SDDMat(NMatSDD,CMatSDD,uNewt,vCount,epsilon).^-1*weight),0,intLength,intLength)^-3;
% correction is a term that shows up in the jacobian. Each row is to be
% scaled by a specific value, we write as a diagon matrix for faster
% multiplication
aproxMAOp = @(u)(pi^2*((SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon).^(-1))*weight).^(-2)+min(min(SDDMat(NMatSDD,CMatSDD,u,vCount,-Inf),epsilon),[],2));


N = length(F);
tol = 1e-16;
APj  = zeros(N,N);
res  = zeros(N,N);
mu   = zeros(N,N);
P    = zeros(N,N);
beta = zeros(N,N);
fu = aproxMAOp(uNewt);
b = fu-F;
P(:,1) = b;
res(:,1) = b;
for j = 1:N
   
   %--------------Compute alpha_j
   APj(:,j) = (aproxMAOp([uNewt(1:N) + delta*P(:,j)/norm(P(:,j));uNewt(N+1:end)])-fu)/delta*norm(P(:,j));
   alpha = dot(res(:,j),APj(:,j))/dot(APj(:,j),APj(:,j));
   
   if isnan(alpha)
       conv_iter = j; 
       break
   elseif (abs(alpha)<tol)
       conv_iter = j;
       break
   end
   
   %------ Compute mu_slow a solution
   mu(:,j+1) =  mu(:,j) + alpha*P(:,j);
   
   %------ Compute residual
   res(:,j+1) =  res(:,j) - alpha*APj(:,j);
       
   %----- Compute Beta_ij
   temp = (aproxMAOp([uNewt(1:N) + delta*APj(:,j)/norm(APj(:,j));uNewt(N+1:end)])-fu)/delta*norm(APj(:,j));
   for i = 1:j
      beta(i,j) = -dot(temp,APj(:,i))/dot(APj(:,i),APj(:,i));
   end
   
   % ----- Compute next basis 
   P(:,j+1) = APj(:,j);
   for i = 1:j
      P(:,j+1) = P(:,j+1) + beta(i,j)*P(:,i);
   end
end
res = res(:,1:conv_iter);
mu = mu(:,1:conv_iter);
% Currently getting warnings about ill conditioned matrices, but the
% corresponding matrix is only used in a matrix product, not any linear
% solves.

% Scheme assumes that the quadrature is the same at each grid point and
% that the stencil is the same at each grid point
deltaU = mu(:,conv_iter);
r = vecnorm(res);
end