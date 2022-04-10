%% Find boundary neighbors
% Should work independently of h, can scale with h afterwards

% Paramaters
N = 2^4+1;
Nx = N; Ny = N;
d = 8;

% Need d < N otherwise we get weight bugs
if d > N
    error('Stencil width is larger than grid')
end

% Grid Points
Xn = 0:Nx+1;
Yn = 0:Ny+1;

% Stencil
[dirs,dcount] = cartStencil(d);
%
Sf = dirs(1:dcount/2,:);
Sb = dirs(dcount/2+1:end,:);

% Near boundary points
nb = [repmat(Xn(2:end-1)',2*d,1) repmat(([1:d (Nx-d+1):Nx])',length(Xn)-2,1); ...
        repmat(([1:d (Nx-d+1):Nx])',length(Yn)-2*d,1) repmat(Yn(d+1:end-d)',2*d,1)];
% where we store the table data, inefficiently concat.
tabf = [];
tabb = [];

% Loop over all directions
for k = 1:dcount/2

    % Slope can be 0 or inf, matlab makes Nans and infs work 
    m = Sf(k,2)/Sf(k,1);

    % find nb points which have this neighbor in the exterior
    indf = nb(:,1)+Sf(k,1) < 0 | nb(:,1)+Sf(k,1)> Nx+1 | nb(:,2)+Sf(k,2)> Ny+1 | nb(:,2)+Sf(k,2) < 0;

    % compute all possible neighbors for each point
    xtf = [zeros(sum(indf),1) (Nx+1)*ones(sum(indf),1) nb(indf,1)-nb(indf,2)/m nb(indf,1)+(Ny+1-nb(indf,2))/m];
    ytf = [nb(indf,2)-m*nb(indf,1), m*(Nx+1)+nb(indf,2)-m*nb(indf,1),zeros(sum(indf),1) (Ny+1)*ones(sum(indf),1)];

    % Correct neighbor has to point in the same direction
    [rf,cf] = find(sign(xtf-nb(indf,1))==sign(Sf(k,1))&sign(ytf-nb(indf,2))==sign(Sf(k,2)));
    xtf = xtf(unique(rf),unique(cf));
    ytf = ytf(unique(rf),unique(cf));

    % Minimize over neighbors which point in the right direction
    df = sqrt((xtf-nb(indf,1)).^2+(ytf-nb(indf,2)).^2);
    [hf,mtf] = min(df,[],2);
    mf = sub2ind(size(df),(1:length(df))',mtf);

    % fill out the table w/ [i,j,k,xp,yp,hp]
    tabf = [tabf;nb(indf,1) nb(indf,2) k*ones(sum(indf),1) xtf(mf) ytf(mf) hf];

    % Same process for backwards (could make function)
    indb = nb(:,1)+Sb(k,1) < 0 | nb(:,1)+Sb(k,1)> Nx+1 | nb(:,2)+Sb(k,2)> Ny+1 | nb(:,2)+Sb(k,2) < 0;
    xtb = [zeros(sum(indb),1) (Nx+1)*ones(sum(indb),1) nb(indb,1)-nb(indb,2)/m nb(indb,1)+(Ny+1-nb(indb,2))/m];
    ytb = [nb(indb,2)-m*nb(indb,1), m*(Nx+1)+nb(indb,2)-m*nb(indb,1),zeros(sum(indb),1) (Ny+1)*ones(sum(indb),1)];
    [rb,cb] = find(sign(xtb-nb(indb,1))==sign(Sb(k,1))&sign(ytb-nb(indb,2))==sign(Sb(k,2)));
    xtb = xtb(unique(rb),unique(cb));
    ytb = ytb(unique(rb),unique(cb));
    db = sqrt((xtb-nb(indb,1)).^2+(ytb-nb(indb,2)).^2);
    [hb,mtb] = min(db,[],2);
    mb = sub2ind(size(db),(1:length(db))',mtb);
    tabb = [tabb;nb(indb,1) nb(indb,2) k*ones(sum(indb),1) xtb(mb) ytb(mb) db(mb)];
end
    

x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
% Can easily map points on [0,N+1] to [a,b]
lx = @(nx) (x1-x0)/(Nx+1)*nx+x0;
ly = @(ny) (y1-y0)/(Ny+1)*ny+y0;


