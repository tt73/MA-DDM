% This is a script that runs schwarz alternating method with 2
% overlapping subdomains.

addpath('Subroutines')

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^5+1;
h = (x1-x0)/(N+1);

depth = 2;        
overlap = 3;  

% PDE
DirBC = @(x,y) (exp((x.^2+y.^2)/2));
contF = @(x,y) ((1+x.^2+y.^2).*exp(x.^2+y.^2));

% Call mesh builder
[Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);

% Call domain decomposition 
[Dom1, Dom2] = splitDomain(Points,Interior,Boundary,h,theta,NMatSDD,depth,overlap,0);

% solver settings
order = 1;
epsilon = h^2;
weight = quadWeights(theta,order);

% DDM settings
max_iter = 1000;
tol = 1e-8; 

% Initialize ddm
u1 = zeros(size(Dom1.l2g));
u2 = zeros(size(Dom2.l2g));

F1 = contF(Points(Dom1.Interior,1),Points(Dom1.Interior,2)); 
uBdry1 = DirBC(Points(Dom1.Boundary,1),Points(Dom1.Boundary,2));  
uBdry1 = [uBdry1; zeros(size(Dom1.Interface))];           

F2 = contF(Points(Dom2.Interior,1),Points(Dom2.Interior,2));  
uBdry2 = DirBC(Points(Dom2.Boundary,1),Points(Dom2.Boundary,2));
uBdry2 = [uBdry2; zeros(size(Dom2.Interface))];          

% DDM Iteration
for k = 1:max_iter
   
    % Step 1: Compute residue 
    res =       norm(uBdry2(end-Dom2.Ns+1:end) - u1(Dom1.send)); 
    res = res + norm(uBdry1(end-Dom1.Ns+1:end) - u2(Dom2.send));
    if (k>1 && res < tol) 
       disp("Converged in "+k+" iterations")
       break
    end
    
    % Step 2: prepare uBdry 
    uBdry2(end-Dom2.Ns+1:end) = u1(Dom1.send); 
    uBdry1(end-Dom1.Ns+1:end) = u2(Dom2.send);
    
    % Step 3: solve 
    [u1, ~] = quadSolver(Dom1.NMatLoc,CMatSDD(Dom1.Interior,:),F1,uBdry1,epsilon,weight,h); % call solver
    [u2, ~] = quadSolver(Dom2.NMatLoc,CMatSDD(Dom2.Interior,:),F2,uBdry2,epsilon,weight,h); % call solver
end


%% Visualize solution 
figure, hold on, view(3)

% plot boundary in black 
plot3(Points(Dom1.Boundary,1),Points(Dom1.Boundary,2),u1(Dom1.Ni+1:Dom1.Ni+Dom1.Nb),'.','Color','k')
plot3(Points(Dom2.Boundary,1),Points(Dom2.Boundary,2),u2(Dom2.Ni+1:Dom2.Ni+Dom2.Nb),'.','Color','k')

% plot interior solutions in different colors 
plot3(Points(Dom1.Interior,1),Points(Dom1.Interior,2),u1(1:Dom1.Ni),'.','Color','b')
plot3(Points(Dom2.Interior,1),Points(Dom2.Interior,2),u2(1:Dom2.Ni),'.','Color','g')

xlabel('x'), ylabel('y'), zlabel('u')
title('DDM Solution')


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

Interface = find(Points(Interior,2)==0);
N = length(Interface); % number of interior pts per row

Interior1 = Interior(Points(Interior,2) <  delta*h);
Interior2 = Interior(Points(Interior,2) > -delta*h);
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

% deterimine which nodes are sent
send1 = zeros(Ns1,1);
send2 = zeros(Ns2,1);
for i = 1:depth
    send1((1:N)+(i-1)*N) = (Ni1-N+1:Ni1)-N*(delta+i);
    send2((1:N)+(i-1)*N) = (1:N)+N*(delta+i);
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

% store everything in Dom1 and Dom2 
Dom1.Interior = Interior1;
Dom1.Boundary = Boundary1;
Dom1.Interface = Interface1;
Dom1.Ni = Ni1;
Dom1.Nb = Nb1;
Dom1.Ns = Ns1;
Dom1.send = send1;
Dom1.l2g = l2g1;
Dom1.NMatLoc = NMatLoc1; 

Dom2.Interior = Interior2;
Dom2.Boundary = Boundary2;
Dom2.Interface = Interface2;
Dom2.Ni = Ni2;
Dom2.Nb = Nb2;
Dom2.Ns = Ns2;
Dom2.send = send2; 
Dom2.l2g = l2g2;
Dom2.NMatLoc = NMatLoc2; 

% this just plots the domain 
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
end

end
