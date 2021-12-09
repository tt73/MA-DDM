%% 1D Newton-Krylov-Schwarz 

a = 0;  b = 1;      % domain
kap = 1;            % curvature
ga = 0;             % u(a) = ga
gb = 1;             % u(b) = gb
u_exact =@(x) -sqrt(1-x.^2)+1; % exact viscosity solution

N = 51;
h = 1/(N+1);
delta = 1;
Nloc = (N+1)/2 + delta - 1;
x_global = h:h:1-h;
U_newton= zeros(N,1);


N_iter = 100;

for l = 1:N_iter
    
    % construct DDM jacobian system
    [J1,J2] = makeJs(U_newton,ga,gb,h,kap,delta); % should be parallel
    
    [b1,b2] = makeRHS(U_newton,ga,gb,h,kap,delta); % should be parallel 
    
    % use orthodir to compute local solutions of s 
    [s1,s2,conv_iter] = run_orthodir(N_iter,1e-12,Nloc,delta,h,J1,J2,b1,b2);
    
    % piece together s1 and s2 to make s 
    s = combine(s1,s2,N,delta);
    
    % update newton iteration
    U_newton = U_newton + s;
end

plot(x_global,U_newton)



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

function [J1,J2] = makeJs(u,ga,gb,h,kap,delta)  % jacobian
N = length(u); 
Nloc = (N+1)/2 + delta - 1; 
ih2 = h^-2; % inverse h^2 

J1 = gallery('tridiag',Nloc,ih2,0,ih2);
J2 = J1;

% left domain
for i = 1:Nloc
    if (i==1)
        m1 = max([u(i)-ga, -u(i+1)+u(i), 0]);
    elseif(i==Nloc)
        m1 = max([u(i)-u(i-1),-u(i+1)+u(i), 0]);
    else
        m1 = max([u(i)-u(i-1),-u(i+1)+u(i), 0]);
    end
    J1(i,i) = -2*ih2 - 3*ih2*kap*(1+ih2*m1^2)^(1/2)*m1;
end

% right domain
shift = N-Nloc;
for i = 1:Nloc
    if (i==1)
        m2 = max([u(i+shift)-u(i-1+shift),-u(i+1+shift)+u(i+shift), 0]);
    elseif(i==Nloc)
        m2 = max([u(i+shift)-u(i-1+shift),-gb+u(i+shift), 0]);
    else
        m2 = max([u(i+shift)-u(i-1+shift),-u(i+1+shift)+u(i+shift), 0]);
    end
    J2(i,i) = -2*ih2 - 3*ih2*kap*(1+ih2*m2^2)^(1/2)*m2;
end
end

function [b1, b2]= makeRHS(u,ga,gb,h,kap,delta) % rhs of oddm jacobi system 
F = makeF(u,ga,gb,h,kap);
N = length(u);
Nloc = (N+1)/2 + delta - 1; 

F1 = F(1:Nloc);
F2 = F(end-Nloc+1:end);
b1 = -F1;
b2 = -F2;
end


function [Ap] = A_projection(p,Nloc,h,delta,D1,D2) % can be parallelized

% construct RHS without the static part
g12 = zeros(Nloc,1);
g21 = g12;
g12(Nloc) = g12(Nloc) + p(1);
g21(1)    = g21(1) + p(2);

% solve system
u1 = D1\g12;
u2 = D2\g21;

% construct A*p
Ap = p;
Ap(1) = Ap(1) + u2(2*delta+2)/h^2;
Ap(2) = Ap(2) + u1(Nloc-(2*delta+1))/h^2;
end

function [u1,u2,conv_iter] = run_orthodir(N_iters,tol,Nloc,delta,h,J1,J2,g1,g2)

% init orthodir
APj  = zeros(2,N_iters);
P = zeros(2,N_iters);
beta = zeros(2,N_iters);
mu = zeros(2,1);
res = zeros(2,1);

% first computation
u1 = J1\g1;
u2 = J2\g2;
res(1) = -u2(2*delta+2)/h^2;
res(2) = -u1(Nloc-(2*delta+1))/h^2;
P(:,1) = res;

for j = 1:N_iters
    
    %--------------Compute alpha_j
    APj(:,j) = A_projection(P(:,j),Nloc,h,delta,J1,J2);
    alpha = dot(res,APj(:,j))/dot(APj(:,j),APj(:,j));
    
    if isnan(alpha)
        conv_iter = j;
        break
    elseif (abs(alpha)<tol)
        conv_iter = j;
        break
    end
    
    %------ Compute mu_slow a solution
    mu =  mu + alpha*P(:,j);
    
    %------ Compute residual
    res =  res - alpha*APj(:,j);
    
    %----- Compute Beta_ij
    temp = A_projection(APj(:,j),Nloc,h,delta,J1,J2);
    
    for i = 1:j
        beta(i,j) = -dot(temp,APj(:,i))/dot(APj(:,i),APj(:,i));
    end
    
    % ----- Compute next basis
    P(:,j+1) = APj(:,j);
    for i = 1:j
        P(:,j+1) = P(:,j+1) + beta(i,j)*P(:,i);
    end
end

% final computation
g12 = zeros(Nloc,1); g12(Nloc) = mu(1);
g21 = zeros(Nloc,1);    g21(1) = mu(2);

% do one last solve with the static part 
u1 = J1\(g1+g12);
u2 = J2\(g2+g21);
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