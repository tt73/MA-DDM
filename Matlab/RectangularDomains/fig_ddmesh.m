
addpath('Subroutines') % the functions are all stored in the subfolder 'Subroutines'
clear
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;

padding = 0.2;
lims = [x0-padding, x1+padding, y0-padding, y1+padding]; % useful for plotting

% Parameters needed to generate grid
N = 9;
h = (x1-x0)/(N+1);
delta = 1;
width = 2;

[Points,Interior,Boundary,NMatSDD1,CMatSDD,theta1] = buildMesh_Rect(x0,x1,y0,y1,h,width);


Interface = find(Points(Interior,2)==0); % points in the interior that lie on the artificial interface
Interior1 = Interior(Points(Interior,2)<0);
Interior2 = Interior(Points(Interior,2)>0);
Boundary1 = Boundary(Points(Boundary,2)<=0);
Boundary2 = Boundary(Points(Boundary,2)>=0);


Boundary1 = Boundary(Points(Boundary,2) <= h*(delta+1)); % extend upward
Boundary2 = Boundary(Points(Boundary,2) >= -h*(delta+1)); % extend downward
Interface1 = [Interface+N*delta; Interface+N*2];
Interface2 = [Interface-N*delta; Interface-N*2];

c1 = [0 0 1];
c2 = [0 1 0];

for i = 1:delta
    Interior1 = [Interior1, Interface'+N*(i-1)];
    Interior2 = [Interior2, Interface'-N*(i-1)];
end



tiledlayout(1,2) % trying something new

nexttile
scatter(Points(:,1),Points(:,2),4,'k')
hold on
scatter(Points(Interior1,1),Points(Interior1,2),15,'r',"filled")
scatter(Points(Boundary1,1),Points(Boundary1,2),20,'b',"filled")
scatter(Points(Interface1,1),Points(Interface1,2),20,"black","filled")
axis(lims)
axis square
title('\Omega_1')

nexttile
scatter(Points(:,1),Points(:,2),4,'k')
hold on
s = gobjects(1,3);
s(1) = scatter(Points(Interior2,1),Points(Interior2,2),15,'r','filled')
s(2) = scatter(Points(Boundary2,1),Points(Boundary2,2),20,'b','filled')
s(3) = scatter(Points(Interface2,1),Points(Interface2,2),20,"black","filled")
hold off
axis(lims)
axis square
title('\Omega_2')

l = legend(s,{'\Omega_i^h','\Gamma_i^h','B_i^h'});
l.Orientation = 'horizontal';
l.Layout.Tile = 'south';

p = get(gcf,'position');
p(3) = p(3)*2; % double the width
set(gcf,'position',p)
