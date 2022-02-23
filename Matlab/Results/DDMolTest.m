clear
addpath('Subroutines Ver1')
% Demonstrate basic convergence of DDM iterations

% 1:Errors = {1:ress, 2:err_exact, 3:err_direct, 4:min_err };
% 2:Iterations = {1:conv_iter, 2:newt_iters, 3:direct_iters};
% 3:Times = {1:ddm_time, 2:direct_time};
% 4:Solutions = {1:u_global, 2:uSoln, 3:uTrue, 4:F, 5:uBdry};
% 5:Domains = {1:Dom1 ,2:Dom2, 3:Points, 4:Interior};
% 6:Parameters = {1:h, 2:depth, 3:tol};

% Fix parameters
method = 2;
tol = 1e-6;
N1 = 2^5+1;
OL1 = [ceil(1+(N1-1)/2*(0:.08:.72))' ceil(1+(N1-1)/2*(0:.08:.72))' ceil(1+(N1-1)/2*(0:.08:.72))';];
Outputs1 = cell(10,3);
for choice = 1:3
    for ol = 1:10
        overlap = ceil(1+(N1-1)/2*(ol-1)*.08);
        [Errors,Iterations,Times,Solutions,Domains,Parameters] = SchwarzSolver(N1,overlap,choice,method,tol);   
        Outputs1{ol,choice} = {Errors,Iterations,Times,Solutions,Domains,Parameters};
    end
end
%%
DDM_Iters1 = zeros(10,3);
DDM_Newts1 = zeros(10,3);
DDM_Time1 = zeros(10,3);
Sub_Size1 = zeros(10,1);
Dir_Size1 = zeros(10,1);
OLP1 = OL1/((N1-1)/2);

for i = 1:3
    for j = 1:10
        DDM_Iters1(j,i) = Outputs1{j,i}{2}{1};
        DDM_Newts1(j,i) = sum(Outputs1{j,i}{2}{2});
        DDM_Time1(j,i) = Outputs1{j,i}{3}{1};
        Sub_Size1(j) = length(Outputs1{j,i}{5}{1}.l2g);
        Dir_Size1(j) = length(Outputs1{j,i}{5}{3});
    end
end

Examples = {'Smooth Radial','Degnerate','Blow-Up'};
figure
hold on
plot(OLP1(:,1),DDM_Iters1(:,1),'-o','LineWidth',1)
plot(OLP1(:,1),DDM_Iters1(:,2),'-*','LineWidth',1)
plot(OLP1(:,1),DDM_Iters1(:,3),'-s','LineWidth',1)
title(sprintf('DDM Solver: Resolution vs Error \n Overlap = 10%%, tol = 1e-6'),'FontSize',20)
legend(Examples,"Location","best","FontSize",11)
xlabel('log10(h)','FontSize',14)
ylabel('log10(||u-u_h||_{\infty})','FontSize',14)
%%
% Fix parameters
method = 2;
tol = 1e-6;
N2 = 2^6+1;
OL2 = [ceil(1+(N2-1)/2*(0:.08:.72))' ceil(1+(N2-1)/2*(0:.08:.72))' ceil(1+(N2-1)/2*(0:.08:.72))';];
Outputs2 = cell(10,3);
for choice = 1:3
    for ol = 1:10
        overlap = ceil(1+(N2-1)/2*(ol-1)*.08);
        [Errors,Iterations,Times,Solutions,Domains,Parameters] = SchwarzSolver(N2,overlap,choice,method,tol);   
        Outputs2{ol,choice} = {Errors,Iterations,Times,Solutions,Domains,Parameters};
    end
end
%%
DDM_Iters2 = zeros(10,3);
DDM_Newts2 = zeros(10,3);
DDM_Time2 = zeros(10,3);
Sub_Size2 = zeros(10,1);
Dir_Size2 = zeros(10,1);
OLP2 = OL2/((N2-1)/2);

for i = 1:3
    for j = 1:10
        DDM_Iters2(j,i) = Outputs2{j,i}{2}{1};
        DDM_Newts2(j,i) = sum(Outputs2{j,i}{2}{2});
        DDM_Time2(j,i) = Outputs2{j,i}{3}{1};
        Sub_Size2(j) = length(Outputs2{j,i}{5}{1}.l2g);
        Dir_Size2(j) = length(Outputs2{j,i}{5}{3});
    end
end


figure
hold on
plot(OLP2(:,1),DDM_Iters2(:,1),'-o','LineWidth',1)
plot(OLP2(:,1),DDM_Iters2(:,2),'-*','LineWidth',1)
plot(OLP2(:,1),DDM_Iters2(:,3),'-s','LineWidth',1)
title(sprintf('DDM Solver: Overlap vs DDM Iterations \n N = %d, tol = 1e-6',N2),'FontSize',20)
legend(Examples,"Location","best","FontSize",11)
xlabel('Percent Overlap','FontSize',14)
ylabel('DDM Iterations','FontSize',14)

figure
hold on
plot(OLP2(:,1),DDM_Newts2(:,1),'-o','LineWidth',1)
plot(OLP2(:,1),DDM_Newts2(:,2),'-*','LineWidth',1)
plot(OLP2(:,1),DDM_Newts2(:,3),'-s','LineWidth',1)
title(sprintf('DDM Solver: Overlap vs Newton Iterations \n N = %d, tol = 1e-6',N2),'FontSize',20)
legend(Examples,"Location","best","FontSize",11)
xlabel('Percent Overlap','FontSize',14)
ylabel('Total Newton Iterations','FontSize',14)

figure
hold on
plot(OLP2(:,1),DDM_Time2(:,1),'-o','LineWidth',1)
plot(OLP2(:,1),DDM_Time2(:,2),'-*','LineWidth',1)
plot(OLP2(:,1),DDM_Time2(:,3),'-s','LineWidth',1)
title('DDM TIME','FontSize',20)
legend(Examples,"Location","best","FontSize",11)
xlabel('Percent Overlap','FontSize',14)
ylabel('Total Newton Iterations','FontSize',14)

figure
hold on
plot(OLP1(:,1),Sub_Size1./Dir_Size1,'-o','LineWidth',1)
plot(OLP2(:,1),Sub_Size2./Dir_Size2,'-*','LineWidth',1)
title('DDM Solver: Overlap vs Subdomains','FontSize',20)
legend({'N = 33','N = 65'},"Location","best","FontSize",11)
xlabel('Percent Overlap','FontSize',14)
ylabel('Relative Subdomain Size','FontSize',14)