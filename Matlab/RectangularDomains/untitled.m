addpath('Subroutines')

%%% Solution to the Dirichlet Problem for MA on a square

x0 = -1; x1 = 1; 
y0 = -1; y1 = 1; 
N = 2^3; h = (x1-x0)/N;
depth = 2; 
% Parameters needed to generate grid

[Points,Interior,Boundary,NMatSDD,CMatSDD,theta] = buildMesh_Rect(x0,x1,y0,y1,h,depth);
% Build the mesh
% Points - (Np x 2) arrary of node coordinates 
% Interior - (1 x Ni) interior node indicator 
% Boundary - (Nb x 1) boundary node indicator 
% NMatSDD - (Ni x ?) 
% CMatSDD - (Ni x ?)
% theta - (? x 1) 
Np = length(Points);
Ni = Interior(end);
Nb = length(Boundary);


figure(1)

xi = Points(Boundary,1);
yi = Points(Boundary,2);

vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for ind = 1:Nb
   hold on 
      scatter(Points(Boundary,1),Points(Boundary,2),[],'b','filled')
      text(xi(ind),yi(ind),num2str(ind))
   hold off
   bb = 0.2;
   axis([x0-bb x1+bb y0-bb y1+bb])
   title('Oredering of Boundary Nodes')
   drawnow
   ff = getframe(gcf);
   for j = 1:10
      writeVideo(vidfile,ff);
   end
end
close(vidfile)
