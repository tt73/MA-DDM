function [uInterp] = InterpInit(x0,x1,y0,y1,DirBC,contF,order,N)


if (~exist('N','var'))
    N = 2^4+1;
end
h = (x1-x0)/(N+1);
depth = ceil(h^(-1/3));
[pCoarse,iCoarse,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);
epsilon = (h*depth)^2;
weight = quadWeights(theta,order);
Dvvs = cell(length(theta),1);

for i = 1:length(theta)
   Dvvs{i} = sparse( repmat(iCoarse,1,3), [NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3)], [CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)], length(iCoarse), length(pCoarse));
end
F = contF(pCoarse(iCoarse,1),pCoarse(iCoarse,2));
uBdry = DirBC(pCoarse(Boundary,1),pCoarse(Boundary,2));
uCoarse = quadSolver2(NMatSDD,CMatSDD,Dvvs,F,uBdry,epsilon,weight,h);
uInterp = scatteredInterpolant(pCoarse(:,1),pCoarse(:,2),full(uCoarse));
end

