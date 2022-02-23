clear
addpath('Subroutines Ver1')
% Demonstrate basic convergence of DDM iterations

% 1:Errors = {1:ress, 2:err_exact, 3:err_direct, 4:min_err };
% 2:Iterations = {1:conv_iter, 2:newt_iters, 3:direct_iters};
% 3:Times = {1:ddm_time, 2:direct_time};
% 4:Solutions = {1:u_global, 2:uSoln, 3:uTrue, 4:F, 5:uBdry};
% 5:Domains = {1:Dom1 ,2:Dom2, 3:Points, 4:Interior};
% 6:Parameters = {1:h, 2:depth, 3:tol};

NN = [(2.^((1:5)+2)+1)' (2.^((1:5)+2)+1)' (2.^((1:5)+2)+1)'];
% Fix parameters
method = 3;
tol = 1e-6;
Outputs = cell(6,3);
for choice = 1:3
    for n = 1:5
        N = 2^(2+n)+1;
        overlap = ceil((N-1)/20);
        [Errors,Iterations,Times,Solutions,Domains,Parameters] = SchwarzSolver(N,overlap,choice,method,tol);   
        Outputs{n,choice} = {Errors,Iterations,Times,Solutions,Domains,Parameters};
    end
end
%%
H = zeros(5,3);
Err = H;
DDM_Iter = H;
DDM_Time = H;
DDM_Newts = H;
Dir_Time = H;
Dir_Count = H;
Grid_Size = H;
for i = 1:3
    for j = 1:5
        H(j,i) = Outputs{j,i}{6}{1};
        DDM_Iter(j,i) = Outputs{j,i}{2}{1};
        Err(j,i) = Outputs{j,i}{1}{2}(DDM_Iter(j,i)-1);
        DDM_Newts(j,i) = sum(Outputs{j,i}{2}{2});
        DDM_Time(j,i) = Outputs{j,i}{3}{1};
        Dir_Time(j,i) = Outputs{j,i}{3}{2};
        Dir_Count(j,i) = Outputs{j,i}{2}{3};

        Grid_Size(j,i) = length(Outputs{j,i}{5}{3});
    end
end

Examples = {'Smooth Radial','Degnerate','Blow-Up'};



figure
hold on
plot(log10(Grid_Size(:,1)),DDM_Newts(:,1),'-o','LineWidth',1.5)
plot(log10(Grid_Size(:,2)),DDM_Newts(:,2),'-*','LineWidth',1.5)
plot(log10(Grid_Size(:,3)),DDM_Newts(:,3),'-s','LineWidth',1.5)
title('Single Step Newton' ,'FontSize',20)
legend(Examples,"Location","best","FontSize",11)
xlabel('log10(Total Grid Points)','FontSize',14)
ylabel('Total Newton Steps','FontSize',14)

