function [subInd,subOL] = fourSquares(x0,x1,y0,y1,delta,Points,Interior)
% Midpoints of domain
xm = (x1+x0)/2; ym = (y1+y0)/2;

% Ind is the subdomain with no overlap
% OL is the indices only in the overlap
% Together the form the interior of the subdomain

subInd = cell(1,4);
subOL = cell(1,4);

% Bottom left 
subInd{1} = find((Points(Interior,1)<xm)&(Points(Interior,2)<ym));
subOL{1} = setdiff(find((Points(Interior,1)<xm+delta)&(Points(Interior,2)<ym+delta)),subInd{1});

% Top left
subInd{2} = find((Points(Interior,1)<xm)&(Points(Interior,2)>=ym));
subOL{2} = setdiff(find((Points(Interior,1)<xm+delta)&(Points(Interior,2)>=ym-delta)),subInd{2});

% Top Right
subInd{3} = find((Points(Interior,1)>=xm)&(Points(Interior,2)>=ym));
subOL{3} = setdiff(find((Points(Interior,1)>=xm-delta)&(Points(Interior,2)>=ym-delta)),subInd{3});

% Bottom Right
subInd{4} = find((Points(Interior,1)>=xm)&(Points(Interior,2)<ym));
subOL{4} = setdiff(find((Points(Interior,1)>=xm-delta)&(Points(Interior,2)<ym+delta)),subInd{4});
end