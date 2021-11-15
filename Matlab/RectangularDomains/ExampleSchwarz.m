% This is a script that runs schwarz alternating method with 2
% overlapping subdomains.

addpath('Subroutines')

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 7;
h = (x1-x0)/(N+1);
depth = 2;        
overlap = 2;  

[Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);

[Dom1, Dom2] = splitDomain(Points,Interior,Boundary,h,depth,overlap,1);




% This function splits the domain along the x-axis Dom1 at the bottom.
% Then it extends the each domain over the x-axis delta layers of nodes.
% Then it labels which nodes are interfaces.
% Finally, it labels which of its own nodes should be sent over.
% Note that delta 
function[Dom1, Dom2] = splitDomain(Points,Interior,Boundary,h,depth,delta,check)

if (delta == 0) 
    error("Input delta must be a positive integer")
end

Interface = find(Points(Interior,2)==0);
N = length(Interface); % number of interior pts per row
thickness = depth-1;   % thickness of interface layer

Interior1 = Interior(Points(Interior,2) <  delta*h);
Interior2 = Interior(Points(Interior,2) > -delta*h);
Boundary1 = Boundary(Points(Boundary,2) <=  (delta+thickness)*h+eps);
Boundary2 = Boundary(Points(Boundary,2) >= -(delta+thickness)*h-eps);

Interface1 = Interface+N*delta;
Interface2 = Interface-N*delta;

for i = 1:thickness
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
for i = 1:thickness
    
end


Dom1.Interior = Interior1;
Dom1.Boundary = Boundary1;
Dom1.Interface = Interface1;
Dom1.Ni = Ni1;
Dom1.Nb = Nb1;
Dom1.Ns = Ns1;
Dom1.send = send1;

Dom2.Interior = Interior2;
Dom2.Boundary = Boundary2;
Dom2.Interface = Interface2;
Dom2.Ni = Ni2;
Dom2.Nb = Nb2;
Dom2.Ns = Ns2;
Dom2.send = send2; 


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