d = 5;
xd = -d:(d+1);
yd = -d:(d+1);
[xx,yy] = meshgrid(xd,yd);
jj = 1:(2*d-1);
M = (d-abs(d-jj))./(d-jj);
alpha = zeros(2*d+2,2*d+2,2*d-1);
for j = jj
    alpha(:,:,j) = xx-yy/M(j);
end

dist = abs(xx-alpha)+abs(yy);
ind1 = find(dist <= d);
ptemp = uniquetol(alpha(ind1));
p = (ptemp((0 <= ptemp)&(ptemp <1)))';

%%
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^4+1;
h = (x1-x0)/(N+1);
xold = x0:h:x1;
yold = y0:h:y1;
m = length(p);
xnew = sort(repmat(xold,1,m))+repmat(p,1,length(xold))*h;
xnew = xnew(1:end-m+1);
ynew = sort(repmat(yold,1,m))+repmat(p,1,length(yold))*h;
ynew = ynew(1:end-m+1);

[X,Y] = meshgrid(xnew,ynew);

[D,~,t] = cartStencil(d);


