%%
% This is a script that runs schwarz alternating method with 2
% overlapping subdomains. The domain is split along the x-axis.

addpath('Subroutines')
clear

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^6+1;
h = (x1-x0)/(N+1);

% requirement: overlap + depth - 1 <= (N-1)/2
depth = ceil(h^(-1/3));
overlap = ceil((N-1)/5);
% overlap = 3;
if (overlap + depth - 1 > (N-1)/2)
   error("overlap + depth exceeds mesh size")
end
% for around 10% overlap, choose overlap ~ (N-1)/20 


% choose F
choice = 1;
switch(choice)
   case 1
      DirBC = @(x,y) (exp((x.^2+y.^2)/2));
      contF = @(x,y) ((1+x.^2+y.^2).*exp(x.^2+y.^2));
   case 2
      pos = @(x) max(x,0);
      DirBC = @(x,y) (.5*pos(((x-.5).^2+(y-.5).^2).^.5-.2).^2);
      contF = @(x,y) (pos(1-.2./sqrt((x-.5).^2+(y-.5).^2)));
   case 3
      DirBC = @(x,y) (-sqrt(2-(x.^2+y.^2)));
      contF = @(x,y) (2*(2-(x.^2+y.^2)).^-2);
end

% Call mesh builder
[Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);

% solver settings
order = 2;
epsilon = (h*depth)^2;
weight = quadWeights(theta,order);
delta = h*overlap;

[subInd,subOL] = fourSquares(x0,x1,y0,y1,delta,Points,Interior);


vCount = length(theta);
aproxMAOp = @(u)(pi^2*((SDDMat(NMatSDD,CMatSDD,u,vCount,epsilon).^(-1))*weight).^(-2)+min(min(SDDMat(NMatSDD,CMatSDD,u,vCount,-Inf),epsilon),[],2));


% DDM settings
max_iter = 200;
conv_iter = max_iter;
% tol = 1e-6;
tol = h;
relax = 1; % must be between 0 and 2, default = 1

% Precompute Dvv matrices
Dvvs = cell(length(theta),1);
for i = 1:length(theta)
   Dvvs{i} = sparse( repmat(Interior,1,3), [NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3)], [CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)], length(Interior), length(Points));
end 

% exact solutions
exact = DirBC(Points(Interior,1),Points(Interior,2));

% solve without doing DDM
F = contF(Points(Interior,1),Points(Interior,2));
uBdry = DirBC(Points(Boundary,1),Points(Boundary,2));

dirtime = tic;
[uSoln, ~,dirCount] = quadSolver2(NMatSDD,CMatSDD,Dvvs,F,uBdry,epsilon,weight,h);
dirCount
toc(dirtime)
% smallest possible error with DDM
min_err = norm(exact-uSoln(Interior),inf)

uDDM = [zeros(length(Interior),1);uBdry];

Doms = subDoms(NMatSDD,CMatSDD,Dvvs,subInd,subOL,F,uDDM);


%% DDM Iteration
% This section can take a while to run.

ress = zeros(max_iter,1);
err_direct = zeros(max_iter,1);
err_exact = zeros(max_iter,1);
newt_iters = zeros(max_iter,1);

ddmtime = tic;
output = cell(1,4);

for k = 1:max_iter
    [out1,newt1] = dom2qs(Doms(1),epsilon,weight,h); % call solver
    [out2,newt2] = dom2qs(Doms(2),epsilon,weight,h); % call solver
    [out3,newt3] = dom2qs(Doms(3),epsilon,weight,h);
    [out4,newt4] = dom2qs(Doms(4),epsilon,weight,h);
    newt_iters(k) = newt1+newt2+newt3+newt4;
    [Doms,uDDM] = updater(out1,out2,out3,out4,Doms,uDDM);
% 

res = norm(aproxMAOp(uDDM)-F,Inf);

    %    % compute error from exact solution
       err_exact(k) = norm(uDDM(Interior)-exact,inf);
    %
    %    % compute discrepancy from direct numerical solver
       err_direct(k) = norm(uDDM(Interior)-uSoln(Interior),inf);
    %
       ress(k) = res;
    if (res < tol)
        conv_iter = k;
        break
    end
end
toc(ddmtime)


%% Residue and Error

figure, hold on, title(sprintf('Error over Iteration with \n N = %d, depth = %d, overlap = %d',N,depth,overlap))
plot(1:conv_iter,log10(ress(1:conv_iter)),'g--','linewidth',1.5)
plot(1:conv_iter,log10(err_direct(1:conv_iter)),'r:','linewidth',1.5)
plot(1:conv_iter,log10(err_exact(1:conv_iter)),'b:','linewidth',1.5)
yline(log10(min_err),'k--')
legend('log10(residue)','log10(err direct)','log10(err exact)','min err')
xlabel('Iterations')


function[output,newt] = dom2qs(Dom,epsilon,weight,h)

[output,~,newt] = quadSolver2(Dom.NMat_l,Dom.CMat_l,Dom.Dvvs_l,Dom.F_l,Dom.uBdry_l,epsilon,weight,h,1,Dom.uInt_l);

end

function[Doms,uDDM] = updater(out1,out2,out3,out4,Doms,uDDM)
    uDDM(Doms(1).nOL_g) = out1(Doms(1).nOL_l);
    uDDM(Doms(2).nOL_g) = out2(Doms(2).nOL_l);
    uDDM(Doms(3).nOL_g) = out3(Doms(3).nOL_l);
    uDDM(Doms(4).nOL_g) = out4(Doms(4).nOL_l);
    for kk = 1:4
        Doms(kk).uInt_l = uDDM(Doms(kk).I_g);
        Doms(kk).uBdry_l = uDDM(Doms(kk).B_g);
    end
end
