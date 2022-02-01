%% 1D NKS 
% The krylov is the outer loop and newton is the inner loop
clear; close all

a = 0;  b = 1;      % domain
kap = 1;            % curvature
ga = 0;             % u(a) = ga
gb = 1;             % u(b) = gb
u_exact =@(x) -sqrt(1-x.^2)+1; % exact viscosity solution

% N = global mesh size 
N = 13;
h = 1/(N+1);
xx_global = a+h:h:b-h;

% delta = how many points beyond center point 
delta = 1;

% size of unknowns on each subdomain
Nloc = (N+1)/2 - 1 + delta; %
xx_1 = xx_global(1:Nloc);
xx_2 = xx_global(end-Nloc+1:end);

% understand which nodes are sent to the neighbor
u1send = (N+1)/2-delta;
u2send = 2*delta;

%% Step 1 - Do one solve using given data on the exterior 

% Do local solve on Omega1 with given BC on left end, 0 on interface
U1 = zeros(Nloc,1);
U1 = fsolve( @(u)makeF(u,ga,0,h,kap),U1); 

% Do local solve on Omega2 with given BC on right end, 0 on interface
U2 = zeros(Nloc,1);
U2 = fsolve( @(u)makeF(u,0,gb,h,kap),U2); 


%% Step 2 - Setup Orthodir 

max_iter = 30;
conv_iter = max_iter;
tol = 1e-16;

dim = 2; % we only send 1 Dirichlet data 

% mu is the desired output
% mu(1) contains the data supplied to Omega1 
% mu(2) contains the data supplied to Omega2

% Need to initialize 
mu   = zeros(dim,1);         
res  = zeros(dim,1);
APj  = zeros(dim,max_iter);
P    = zeros(dim,max_iter);
beta = zeros(dim,max_iter);

% initialize with residue and P 
res(1) = U2(u2send); % data supplied to Omega1
res(2) = U1(u1send); % data supplied to Omega2
P(:,1) = res;

if(0)
    figure(1), hold on
    plot(xx_global,u_exact(xx_global),'ko-','LineWidth',1)
    plot(xx_1,U1,'bo-','LineWidth',1.5)
    plot(xx_2,U2,'ro-','LineWidth',1.5)
    plot(xx_1(u1send),U1(u1send),'ro')
    plot(xx_2(u2send),U2(u2send),'bo')
end


%% Step 3 - The orthodir iteration

ff = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

for j = 1:max_iter
    
    %--------------Compute alpha_j
    APj(:,j) = A_projection(P(:,j),Nloc,h,u1send,u2send,kap,ga,gb);
    
    alpha = dot(res,APj(:,j))/dot(APj(:,j),APj(:,j));
    
    if isnan(alpha)
        conv_iter = j;
        break
    elseif (abs(alpha)<tol)
        conv_iter = j;
        break
    end
    
    %------ Compute mu 
    mu = mu + alpha*P(:,j);
    
    %------ Compute residual
    res = res - alpha*APj(:,j);
    
    %----- Compute Beta_ij
    temp = A_projection(APj(:,j),Nloc,h,u1send,u2send,kap,ga,gb);
    
    for i = 1:j
        beta(i,j) = -dot(temp,APj(:,i))/dot(APj(:,i),APj(:,i));
    end
    
    % ----- Compute next basis
    P(:,j+1) = APj(:,j);
    for i = 1:j
        P(:,j+1) = P(:,j+1) + beta(i,j)*P(:,i);
    end
    
    if(1)
        U1 = zeros(Nloc,1);
        U1 = fsolve( @(u)makeF(u,ga,mu(1),h,kap),U1);
        U2 = zeros(Nloc,1);
        U2 = fsolve( @(u)makeF(u,mu(2),gb,h,kap),U2);
        
        hold on
        plot(xx_global,u_exact(xx_global),'ko-','LineWidth',1)
        
        plot([a,xx_1,xx_1(end)+h],[ga;U1;mu(1)],'bo-','LineWidth',1.5)
        plot([a,xx_1(end)+h],[ga;mu(1)],'bo','LineWidth',1.5,'MarkerFaceColor','b')
        
        plot([xx_2(1)-h,xx_2,b],[mu(2);U2;gb],'ro-','LineWidth',1.5)
        plot([xx_2(1)-h,b],[mu(2);gb],'ro','LineWidth',1.5,'MarkerFaceColor','r')

        plot(xx_1(u1send),U1(u1send),'rx')
        plot(xx_2(u2send),U2(u2send),'bx')

        ylim([-1,1])
        xlim([0,1])
        title(sprintf('j = %3d',j))
        
        % Draw plot 
        drawnow
        
        % Capture the plot as an image
        frame = getframe(ff);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if j == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        
        cla % clear the figure
    end
    
end



%% Step 4 - One last solve with external data and output of orthodir 

% Do local solve on Omega1 with given BC on left end, P(1) on interface
U1 = zeros(Nloc,1);
U1 = fsolve( @(u)makeF(u,ga,P(1,conv_iter),h,kap),U1); 

% Do local solve on Omega2 with given BC on right end, 0 on interface
U2 = zeros(Nloc,1);
U2 = fsolve( @(u)makeF(u,P(2,conv_iter),gb,h,kap),U2); 

if(0) 
   figure(3)
   plot(xx_1,U1,xx_2,U2)
end


%% Subroutines 

function F = makeF(u,ga,gb,h,kap) % system to find root of
N = length(u);
F = zeros(N,1);
ih2 = 1/h^2;

F(1) = ih2*(u(2)-2*u(1)+ga)-kap*(1+ih2*max([u(1)-ga,u(1)-u(2),0])^2)^(3/2); % at the left endpoint

for j = 2:N-1
    F(j) = ih2*(u(j+1)-2*u(j)+u(j-1))-kap*(1+ih2*max([u(j)-u(j-1),u(j)-u(j+1),0])^2)^(3/2); % at the middle points
end

F(N) = ih2*(gb-2*u(N)+u(N-1))-kap*(1+ih2*max([u(N)-u(N-1),u(N)-gb,0])^2)^(3/2); % at the right endpoint
end


function [Ap] = A_projection(p,Nloc,h,u1send,u2send,kap,ga,gb) % can be parallelized

% p(1) has data for Omega1, p(2) has data for Omega2
U1 = zeros(Nloc,1);
U2 = zeros(Nloc,1);

zeroD = 1; 
if(zeroD)
    U1 = fsolve( @(u)makeF(u,0,p(1),h,kap),U1); 
    U2 = fsolve( @(u)makeF(u,p(2),0,h,kap),U2); 
else
    U1 = fsolve( @(u)makeF(u,ga,p(1),h,kap),U1); 
    U2 = fsolve( @(u)makeF(u,p(2),gb,h,kap),U2); 
end

% construct A*p
Ap = p;
if(zeroD)
    Ap(1) = p(1) + U2(u2send);
    Ap(2) = p(2) + U1(u1send);
else
    Ap(1) =  U2(u2send);
    Ap(2) =  U1(u1send);
end

end



function s = combine(s1,s2,N,delta)

s = zeros(N,1);

% transfer solution from s1 and s2 to s
nn = (N+1)/2-1; % number of non-overlapping interior nodes 
s(1:nn) = s1(1:nn);
s(end-nn+1:end) = s2(end-nn+1:end); 


mid = (N+1)/2;
s(mid) = (s1(end-delta+1)+s2(delta))/2; 
end