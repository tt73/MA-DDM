addpath('Subroutines') % the functions are all stored in the subfolder 'Subroutines'
clear
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;

padding = 0.2;
lims = [x0-padding, x1+padding, y0-padding, y1+padding]; % useful for plotting

% Parameters needed to generate grid
N = 2^3+1;
h = (x1-x0)/(N+1);

% 3
[Points,Interior,Boundary,NMatSDD1,CMatSDD,theta1] = buildMesh_Rect(x0,x1,y0,y1,h,3);

figure
scatter(Points(:,1),Points(:,2),4,'k')
hold on
axis(lims)
axis off

ind = 2*N + 10; row = NMatSDD1(ind,:);
for i = 1:length(theta1)
   pts = Points( row((1:3)+(i-1)*3),:);
   plot(pts(:,1),pts(:,2),'k','LineWidth',1.5)
   scatter(pts(:,1),pts(:,2),20,'k','filled')
end

axis square
