%%
% This is a script that runs schwarz alternating method with 2
% overlapping subdomains. The domain is split along the x-axis.

addpath('Subroutines')
clear

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^4+1;
h = (x1-x0)/(N+1);

% requirement: overlap + depth - 1 <= (N-1)/2
depth = 1;
overlap = 1;
if (overlap + depth - 1 > (N-1)/2)
    error("overlap + depth exceeds mesh size")
end
% for around 10% overlap, choose overlap ~ (N-1)/20
fprintf("Overlap:   %4.2f%%\n",overlap/((N-1)/2)*100)


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
order = 1;
epsilon = h^2;
weight = quadWeights(theta,order);

% Call domain decomposition
[Dom1, Dom2] = splitDomain(Points,Interior,Boundary,h,theta,NMatSDD,depth,overlap,0);

% DDM settings
max_iter = 100;
conv_iter = max_iter;
tol = 1e-6;
relax = 1.0; % must be between 0 and 2, default = 1

% Initialize ddm
u1 = zeros(size(Dom1.l2g));
u2 = zeros(size(Dom2.l2g));

F1 = contF(Points(Dom1.Interior,1),Points(Dom1.Interior,2));
uBdry1 = DirBC(Points(Dom1.Boundary,1),Points(Dom1.Boundary,2));
uBdry1 = [uBdry1; zeros(size(Dom1.Interface))];

F2 = contF(Points(Dom2.Interior,1),Points(Dom2.Interior,2));
uBdry2 = DirBC(Points(Dom2.Boundary,1),Points(Dom2.Boundary,2));
uBdry2 = [uBdry2; zeros(size(Dom2.Interface))];

% Precompute Dvv matrices
Dvvs = cell(length(theta),1);
Dvvs1 = cell(length(theta),1);
Dvvs2 = cell(length(theta),1);
for i = 1:length(theta)
    Dvvs{i} = sparse( repmat(Interior,1,3), [NMatSDD(:,i*3-2) NMatSDD(:,i*3-1) NMatSDD(:,i*3)], [CMatSDD(:,i*3-2) CMatSDD(:,i*3-1) CMatSDD(:,i*3)], length(Interior), length(Points));
    Dvvs1{i} = sparse( repmat(Dom1.Interior,1,3), [Dom1.NMatLoc(:,i*3-2) Dom1.NMatLoc(:,i*3-1) Dom1.NMatLoc(:,i*3)], [CMatSDD(Dom1.Interior,i*3-2) CMatSDD(Dom1.Interior,i*3-1) CMatSDD(Dom1.Interior,i*3)], Dom1.Ni, Dom1.Ni+Dom1.Nb+Dom1.Ns);
    Dvvs2{i} = sparse( repmat(Dom1.Interior,1,3), [Dom2.NMatLoc(:,i*3-2) Dom2.NMatLoc(:,i*3-1) Dom2.NMatLoc(:,i*3)], [CMatSDD(Dom2.Interior,i*3-2) CMatSDD(Dom2.Interior,i*3-1) CMatSDD(Dom2.Interior,i*3)], Dom2.Ni, Dom2.Ni+Dom2.Nb+Dom2.Ns);
end


%% Error in the direct solution

% exact solutions
exact = DirBC(Points(Interior,1),Points(Interior,2));
exact1 = DirBC(Points(Dom1.Interior,1),Points(Dom1.Interior,2));
exact2 = DirBC(Points(Dom2.Interior,1),Points(Dom2.Interior,2));

% solve without doing DDM
F = contF(Points(Interior,1),Points(Interior,2));
uBdry = DirBC(Points(Boundary,1),Points(Boundary,2));
[uSoln, ~] = quadSolver2(NMatSDD,CMatSDD,Dvvs,F,uBdry,epsilon,weight,h);

% smallest possible error with DDM
min_err = norm(exact-uSoln(Interior),inf);

% direct solution
direct1 = uSoln(Dom1.Interior);
direct2 = uSoln(Dom2.Interior);

% init with exact solution
% u1 = DirBC(Points(Dom1.l2g,1),Points(Dom1.l2g,2));
% u2 = DirBC(Points(Dom2.l2g,1),Points(Dom2.l2g,2));





%% Orthodir
% The solutions on the interface being sent make up the vector that is
% being solved.

ns1 = length(Dom1.send);
ns2 = length(Dom2.send);
odir_dim = ns1 + ns2;

mu = zeros(odir_dim,1);
APj  = zeros(odir_dim,max_iter);
P = zeros(odir_dim,max_iter);
beta = zeros(odir_dim,max_iter);

% first computation
[out1,~,~] = quadSolver2(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),Dvvs1,F1,uBdry1,epsilon,weight,h,u1(1:length(Dom1.Interior))); % call solver
[out2,~,~] = quadSolver2(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),Dvvs2,F2,uBdry2,epsilon,weight,h,u2(1:length(Dom2.Interior))); % call solver

% initialize `res` and first `P`
b = zeros(odir_dim,1);
b(1:ns1) = out1(Dom1.send);
b(ns1+1:ns1+ns2) = out2(Dom2.send);
res = b;
P(:,1) = b;

ress = zeros(max_iter,1);
err_direct = zeros(max_iter,1);
err_exact = zeros(max_iter,1);
newt_iters = zeros(max_iter,1);

ttt = tic;
for k = 1:max_iter
    
    
    %--------------Compute alpha_j
    APj(:,k) = project(P(:,k),Dom1,Dom2,CMatSDD,Dvvs1,Dvvs2,F1,F2,epsilon,weight,h);
    alpha = dot(res,APj(:,k))/dot(APj(:,k),APj(:,k));
    
    if isnan(alpha)
        conv_iter = k;
        break
    elseif (abs(alpha)<tol)
        conv_iter = k;
        break
    end
    
    %------ Compute mu_slow a solution
    mu =  mu + alpha*P(:,k);
    
    %------ Compute residual
    res =  res - alpha*APj(:,k);
    
    %----- Compute Beta_ij
    temp = project(APj(:,k),Dom1,Dom2,CMatSDD,Dvvs1,Dvvs2,F1,F2,epsilon,weight,h);
    
    for i = 1:k
        beta(i,k) = -dot(temp,APj(:,i))/dot(APj(:,i),APj(:,i));
    end
    
    % ----- Compute next basis
    P(:,k+1) = APj(:,k);
    for i = 1:k
        P(:,k+1) = P(:,k+1) + beta(i,k)*P(:,i);
    end
    
    ress(k) = norm(res);
end
ttime = toc(ttt);

%% one last solve


uBdry1(end-Dom1.Ns+1:end) = mu(Dom1.Ns+1:end);
uBdry2(end-Dom2.Ns+1:end) = mu(1:Dom1.Ns);

[u1,~,~] = quadSolver2(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),Dvvs1,F1,uBdry1,epsilon,weight,h,u1(1:length(Dom1.Interior))); % call solver
[u2,~,~] = quadSolver2(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),Dvvs2,F2,uBdry2,epsilon,weight,h,u2(1:length(Dom2.Interior))); % call solver



    % compute the global DDM solution
    u_interface = (u1(Dom1.ai) + u2(Dom2.ai))/2;
    u_global = [u1(Dom1.oi); u_interface; u2(Dom2.oi)];
    


%% info

if (conv_iter < max_iter)
    disp("Conv:    "+k+" iterations")
else
    disp("Conv:    hit max iterations")
end

%% Visualize solution
figure, hold on, view(3)

% plot boundary in black
plot3(Points(Dom1.Boundary,1),Points(Dom1.Boundary,2),u1(Dom1.Ni+1:Dom1.Ni+Dom1.Nb),'.','Color','k')
plot3(Points(Dom2.Boundary,1),Points(Dom2.Boundary,2),u2(Dom2.Ni+1:Dom2.Ni+Dom2.Nb),'.','Color','k')

% plot the original interior (oi) solutions for each domain in different colors
plot3(Points(Dom1.Interior(Dom1.oi),1),Points(Dom1.Interior(Dom1.oi),2),u1(Dom1.oi),'.','Color','b')
plot3(Points(Dom2.Interior(Dom2.oi),1),Points(Dom2.Interior(Dom2.oi),2),u2(Dom2.oi),'.','Color','g')

% plot the interface - average of both domains
u_interface = (u1(Dom1.ai) + u2(Dom2.ai))/2;
plot3(Points(Dom1.Interior(Dom1.ai),1),Points(Dom1.Interior(Dom1.ai),2),u_interface,'.','Color','m')
xlabel('x'), ylabel('y'), zlabel('u')
title('DDM Solution')


%% Residue and Error

figure, hold on, title(sprintf('Error over Iteration with \n N = %d, depth = %d, overlap = %d',N,depth,overlap))
plot(1:conv_iter,log10(ress(1:conv_iter)),'g--','linewidth',1.5)
plot(1:conv_iter,log10(err_direct(1:conv_iter)),'r:','linewidth',1.5)
plot(1:conv_iter,log10(err_exact(1:conv_iter)),'b:','linewidth',1.5)
yline(log10(min_err),'k--')
legend('log10(residue)','log10(err direct)','log10(err exact)','min err')
xlabel('Iterations')



%%
function Ap = project(p,Dom1,Dom2,CMatSDD,Dvvs1,Dvvs2,F1,F2,epsilon,weight,h)

% Ap
Ap = p;

% prepare uBdry
uBdry1 = zeros(Dom1.Nb+Dom1.Ns,1);
uBdry2 = zeros(Dom2.Nb+Dom2.Ns,1);

uBdry2(end-Dom2.Ns+1:end) = p(Dom1.Ns+1:end);
uBdry1(end-Dom1.Ns+1:end) = p(1:Dom1.Ns);

% solve
[out1,~,~] = quadSolver2(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),Dvvs1,F1,uBdry1,epsilon,weight,h); % call solver
[out2,~,~] = quadSolver2(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),Dvvs2,F2,uBdry2,epsilon,weight,h);
Ap(1:Dom1.Ns) = Ap(1:Dom1.Ns) + out1(Dom1.send);
Ap(Dom1.Ns+1:end) = Ap(Dom1.Ns+1:end) + out2(Dom2.send);

end


%%
% This function splits the domain along the x-axis Dom1 at the bottom.
% Then it extends the each domain over the x-axis delta layers of nodes.
% Then it labels which nodes are interfaces.
% Finally, it labels which of its own nodes should be sent over.
% Note that delta has to be > 0.
function[Dom1, Dom2] = splitDomain(Points,Interior,Boundary,h,theta,NMat,depth,delta,check)

if (delta == 0)
    error("Input delta must be a positive integer")
end

% interior pts along the line y=0 are the original interface
Interface = find(Points(Interior,2)==0);
N = length(Interface); % number of interior pts per row

% identify node type y-coordinates and multiples of h
Interior1 = Interior(Points(Interior,2) <  delta*h-eps); % epsilons are needed for rounding error
Interior2 = Interior(Points(Interior,2) > -delta*h+eps);
Boundary1 = Boundary(Points(Boundary,2) <=  (delta+depth-1)*h+eps);
Boundary2 = Boundary(Points(Boundary,2) >= -(delta+depth-1)*h-eps);

Interface1 = Interface+N*delta;
Interface2 = Interface-N*delta;

for i = 1:depth-1
    Interface1 = [Interface1; Interface+N*(delta+i)];
    Interface2 = [Interface2; Interface-N*(delta+i)];
end

Ni1 = length(Interior1);
Nb1 = length(Boundary1);
Ns1 = length(Interface1);
Ni2 = length(Interior2);
Nb2 = length(Boundary2);
Ns2 = length(Interface2);


% store everything in Dom1 and Dom2
Dom1.Interior = Interior1;
Dom1.Boundary = Boundary1;
Dom1.Interface = Interface1;
Dom1.Ni = Ni1;
Dom1.Nb = Nb1;
Dom1.Ns = Ns1;
Dom2.Interior = Interior2;
Dom2.Boundary = Boundary2;
Dom2.Interface = Interface2;
Dom2.Ni = Ni2;
Dom2.Nb = Nb2;
Dom2.Ns = Ns2;

% this just plots the domains
if(check)
    c1 = 'b'; c2 = 'g';
    figure
    subplot(121)
    scatter(Points(:,1),Points(:,2),1,'k')
    hold on
    scatter(Points(Interior1,1),Points(Interior1,2),[],c1)
    scatter(Points(Boundary1,1),Points(Boundary1,2),10,c1,"filled")
    scatter(Points(Interface1,1),Points(Interface1,2),10,"black","filled")
    for i = 1:length(Interior1)
        text(Points(Interior1(i),1),Points(Interior1(i),2),num2str(i))
    end
    for i = 1:length(Interface1)
        text(Points(Interface1(i),1),Points(Interface1(i),2),num2str(i),'color','r')
    end
    title('\Omega_1^{*}')
    
    subplot(122)
    scatter(Points(:,1),Points(:,2),1,'k')
    hold on
    scatter(Points(Interior2,1),Points(Interior2,2),[],c2)
    scatter(Points(Boundary2,1),Points(Boundary2,2),10,c2,'filled')
    scatter(Points(Interface2,1),Points(Interface2,2),10,"black","filled")
    for i = 1:length(Interior2)
        text(Points(Interior2(i),1),Points(Interior2(i),2),num2str(i))
    end
    for i = 1:length(Interface2)
        text(Points(Interface2(i),1),Points(Interface2(i),2),num2str(i),'color','r')
    end
    title('\Omega_2^{*}')
    sgtitle("Local ordering of Nodes")
    set(gcf,'position',[1001  559  2*560  420])
    
    disp("Check the plot. Press any key to continue.")
    pause
    disp("Resuming...")
end

% This is just for plotting:
% original interior (oi) is the index of the interior nodes if delta was 0
oi1 = 1:N*(N-1)/2;
oi2 = oi1 + N*delta;
% artificial interface (ai) is the index of the x-axis nodes
ai1 = (Ni1-N+1:Ni1) - (delta-1)*N;
ai2 = (1:N) + (delta-1)*N;


% deterimine which nodes are sent
send1 = zeros(Ns1,1);
send2 = zeros(Ns2,1);
for i = 1:depth
    send1((1:N)+(i-1)*N) = (Ni1-N+1:Ni1)-N*(2*delta-2+i);
    send2((1:N)+(i-1)*N) = (1:N)+N*(2*delta-2+i);
end

% the local-to-global index mapping
l2g1 = [Interior1'; Boundary1; Interface1];
l2g2 = [Interior2'; Boundary2; Interface2];

% the local NMat matrix
Ns = length(theta);
NMatLoc1 = NMat(Interior1,:);
for i = 1:Ni1
    for j = 1:3*Ns
        ind = NMat(i,j);
        loc = find(ind == l2g1);
        NMatLoc1(i,j) = loc;
    end
end
NMatLoc2 = NMat(Interior2,:);
for i = 1:Ni2
    for j = 1:3*Ns
        ind = NMatLoc2(i,j);
        loc = find(ind == l2g2);
        NMatLoc2(i,j) = loc;
    end
end

% store the rest of the info
Dom1.send = send1;
Dom1.oi = oi1;
Dom1.ai = ai1;
Dom1.l2g = l2g1;
Dom1.NMatLoc = NMatLoc1;
Dom2.send = send2;
Dom2.oi = oi2;
Dom2.ai = ai2;
Dom2.l2g = l2g2;
Dom2.NMatLoc = NMatLoc2;
end
