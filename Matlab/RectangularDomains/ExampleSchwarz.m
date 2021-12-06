%%
% This is a script that runs schwarz alternating method with 2
% overlapping subdomains. The domain is split along the x-axis.

addpath('Subroutines')
clear

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^7+1;
h = (x1-x0)/(N+1);

% requirement: overlap + depth - 1 <= (N-1)/2
depth = 3;
overlap = 10;
if (overlap + depth - 1 > (N-1)/2)
   error("overlap + depth exceeds mesh size")
end

% choose F
choice = 3;
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

%% Error in the direct solution

% exact solutions
exact = DirBC(Points(Interior,1),Points(Interior,2));
exact1 = DirBC(Points(Dom1.Interior,1),Points(Dom1.Interior,2));
exact2 = DirBC(Points(Dom2.Interior,1),Points(Dom2.Interior,2));

% solve without doing DDM
F = contF(Points(Interior,1),Points(Interior,2));
uBdry = DirBC(Points(Boundary,1),Points(Boundary,2));
[uSoln, ~] = quadSolver(NMatSDD,CMatSDD,F,uBdry,epsilon,weight,h);

% smallest possible error with DDM
min_err = norm(exact-uSoln(Interior),inf);

% direct solution
direct1 = uSoln(Dom1.Interior);
direct2 = uSoln(Dom2.Interior);

% init with exact solution
% u1 = DirBC(Points(Dom1.l2g,1),Points(Dom1.l2g,2));
% u2 = DirBC(Points(Dom2.l2g,1),Points(Dom2.l2g,2));

%% DDM Iteration
% This section can take a while to run.

ress = zeros(max_iter,1);
err_direct = zeros(max_iter,1);
err_exact = zeros(max_iter,1);
newt_iters = zeros(max_iter,1);

for k = 1:max_iter
   
   % prepare uBdry
   uBdry2(end-Dom2.Ns+1:end) = u1(Dom1.send);
   uBdry1(end-Dom1.Ns+1:end) = u2(Dom2.send);
   
   % solve
   if (k>1)
      [out1,~,newt1] = quadSolver(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),F1,uBdry1,epsilon,weight,h,u1(1:length(Dom1.Interior))); % call solver
      [out2,~,newt2] = quadSolver(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),F2,uBdry2,epsilon,weight,h,u2(1:length(Dom2.Interior))); % call solver
      newt_iters(k) = newt1+newt2;
   else
      [out1,~,newt1] = quadSolver(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),F1,uBdry1,epsilon,weight,h);
      [out2,~,newt2] = quadSolver(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),F2,uBdry2,epsilon,weight,h);
      newt_iters(k) = newt1+newt2;
   end
   
   % residue
   res = norm(uBdry2(end-Dom2.Ns+1:end) - out1(Dom1.send)) + ...
      norm(uBdry1(end-Dom1.Ns+1:end) - out2(Dom2.send));
   if(k == 1)
      res0 = res;
   end
   res = res/res0;
   if (res < tol)
      conv_iter = k;
      break
   end
   
   % relaxation, relax = 1 is the regular schwarz method
   u1 = relax*out1 + (1-relax)*u1;
   u2 = relax*out2 + (1-relax)*u2;
   
   % compute the global DDM solution
   u_interface = (u1(Dom1.ai) + u2(Dom2.ai))/2;
   u_global = [u1(Dom1.oi); u_interface; u2(Dom2.oi)];
   
   % compute error from exact solution
   err_exact(k) = norm(u_global-exact,inf);
   
   % compute discrepancy from direct numerical solver
   err_direct(k) = norm(u_global-uSoln(Interior),inf);
   
   ress(k) = res;
end

%% info

fprintf("Overlap: %4.2f%%\n",overlap/((N-1)/2)*100)

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
