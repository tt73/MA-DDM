function[subInd,subOL] = rectSubs(x0,x1,y0,y1,Points,delta,Nx,Ny)
% Set of corners of subdomains
xCorns = linspace(x0,x1,Nx+1);
yCorns = linspace(y0,y1,Ny+1);

% subInd is the global indices of the points in each disjoint subdomain
% subOL is the global indices of the exterior overlap of each subdomain
subInd = cell(Ny,Nx);
subOL = subInd;

for i = 1:Ny
    % top and bottom bounds
    yB = yCorns(i); yT = yCorns(i+1);
    % yInd is a Logical array of points which fall into the [yB,yT] strip
    % yOL is a logical array of points which fall into [yB-delta,yT+delta]
    if i  == 1
        yInd = Points(:,2) <= yT;
        yOL = Points(:,2) <= yT+delta;
    else
        yInd = (Points(:,2) <= yT) & (yB < Points(:,2));
        yOL = (Points(:,2) <= yT+delta)&(yB-delta < Points(:,2));
    end
    for j = 1:Nx
        % left and right bounds
        xL = xCorns(j); xR = xCorns(j+1);
        % same as y but horizontal strip
        if j  == 1
            xInd = Points(:,1) <= xR;
            xOL = Points(:,1) <= xR+delta;
        else
            xInd = (Points(:,1) <= xR) & (xL < Points(:,1));
            xOL = (Points(:,1) <= xR+delta)&(xL-delta < Points(:,1));
        end
        
        % Intersection of the strips
        subInd{i,j} = find(yInd & xInd);
        % Overlap minus the disjoint subdomain
        subOL{i,j} = setdiff(find(yOL&xOL),subInd{i,j});
    end
end
end