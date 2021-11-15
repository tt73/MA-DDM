addpath('Subroutines')

%%% Solution to the Dirichlet Problem for MA on a square

x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^3; h = (x1-x0)/N;
depth = 3;
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


% figure
% xi = Points(Boundary,1);
% yi = Points(Boundary,2);
% vidfile = VideoWriter('testmovie1.mp4','MPEG-4');
% open(vidfile);
% for ind = 1:Nb
%     hold on
%     scatter(Points(Boundary,1),Points(Boundary,2),[],'b','filled')
%     text(xi(ind),yi(ind),num2str(ind))
%     hold off
%     bb = 0.2;
%     axis([x0-bb x1+bb y0-bb y1+bb])
%     title('Oredering of Boundary Nodes')
%     drawnow
%     ff = getframe(gcf);
%     for j = 1:10
%         writeVideo(vidfile,ff);
%     end
% end
% close(vidfile)



ind = N + 1; 
row = NMatSDD(ind,:); 

figure
vidfile = VideoWriter('testmovie2.mp4','MPEG-4');
open(vidfile);
for i = 1:length(theta)
    scatter(Points(:,1),Points(:,2),1,'k')
    hold on
    c = rand(1,3); % random color
    pts = Points( row((1:3)+(i-1)*3),:);
    scatter(pts(:,1),pts(:,2),20,c,'filled')
    hold off
    title("Angle: "+i)
    drawnow
    ff = getframe(gcf);
    for j = 1:10
        writeVideo(vidfile,ff);
    end
end

pts = Points(row,:);
close(vidfile)