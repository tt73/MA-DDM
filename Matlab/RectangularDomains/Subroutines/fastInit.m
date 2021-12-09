function [uPrec,stepCount,resid,dr] = fastInit(NMatSDD,CMatSDD,F,uBdry,epsilon,weight,h,uInit)
% 


% number of points in the quadrature
vCount = length(weight);


if nargin == 7
    
    uInit = poissonInit(NMatSDD,CMatSDD,F,uBdry,1,(vCount+1)/2);
    % Does a single step of the poisson iteration to initialize the solver
    
end

uPrec = [uInit; uBdry];

rNew = Inf;
dr = -Inf;
linF = sqrt(F);
intLength = length(uInit);
pLength = length(uBdry)+intLength;
Boundary = (intLength+1):pLength;

aproxMAOp = @(u)(pi^2*((SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon).^(-1))*weight).^(-2));
stepCount = 1;
stepMax = 50;
while (rNew > h) && (dr <0) && (stepCount < stepMax)

    rOld = rNew;
Mat = SDDMat(NMatSDD,CMatSDD,uPrec,vCount,epsilon);
[L2,Ind2] = max(Mat,[],2);
[L1,Ind1] = min(Mat,[],2);
z=(sqrt(L1.*L2) -L1)./(L2-L1);
z(isnan(z)) = 1/2;


D1 = sparse(intLength,pLength);
D2 = D1;

for i = 1:vCount
    Dvv = sparse(repmat(1:intLength,1,3),[NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3) ],[CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)],intLength,pLength);
    Temp1 = find(Ind1 == i);
    D1(Temp1,:) = Dvv(Temp1,:);
    Temp2 = find(Ind2 == i);
    D2(Temp2,:) = Dvv(Temp2,:);
end
D = spdiags((1-z),0,intLength,intLength)*D1+spdiags(z,0,intLength,intLength)*D2;

rho = linF-D(:,Boundary)*uBdry;
Dsolve = D(1:intLength,1:intLength);
uPrec = [Dsolve\rho; uBdry];
rNew = norm(linF-sqrt(aproxMAOp(uPrec)),Inf);
dr = rNew-rOld;
stepCount = stepCount+1;
end
resid = rNew;
end

