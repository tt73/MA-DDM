clear
addpath('Subroutines Ver1')
% Demonstrate basic convergence of DDM iterations for each representative
% example

% 1:Errors = {1:ress, 2:err_exact, 3:err_direct, 4:min_err };
% 2:Iterations = {1:conv_iter, 2:newt_iters, 3:direct_iters};
% 3:Times = {1:ddm_time, 2:direct_time};
% 4:Solutions = {1:u_global, 2:uSoln, 3:uTrue, 4:F, 5:uBdry};
% 5:Domains = {1:Dom1 ,2:Dom2, 3:Points, 4:Interior};
% 6:Parameters = {1:h, 2:depth, 3:tol};

% Fix parameters
N = 2^5+1;
overlap = ceil((N-1)/20);
tol = 1e-6;
choice = 1;

method = 1;
    [Errors1,Iterations1,Times1,Solutions1,Domains1,Parameters1] = SchwarzSolver(N,overlap,choice,method,tol);   
    conv_iter1 = Iterations1{1};
    ress1 = Errors1{1};
    err_exact1 = Errors1{2};
    min_err = Errors1{4};
    depth = Parameters1{2};
method = 2;
    [Errors2,Iterations2,Times2,Solutions2,Domains2,Parameters2] = SchwarzSolver(N,overlap,choice,method,tol);   
    conv_iter2 = Iterations2{1};
    ress2 = Errors2{1};
    err_exact2 = Errors2{2};

    figure
    subplot(2,1,1), hold on
    title('Standard Initialization')
    plot(1:conv_iter1,log10(ress1(1:conv_iter1)),'g--','linewidth',1.5)
    plot(1:conv_iter1,log10(err_exact1(1:conv_iter1)),'b:','linewidth',1.5)
    yline(log10(min_err),'k--')
    legend('log10(residue)','log10(err exact)','min err')
    xlabel('Iterations')
subplot(2,1,2), hold on
    title('Last Iteration Initialization')
    plot(1:conv_iter2,log10(ress2(1:conv_iter2)),'g--','linewidth',1.5)
    plot(1:conv_iter2,log10(err_exact2(1:conv_iter2)),'b:','linewidth',1.5)
    yline(log10(min_err),'k--')
    legend('log10(residue)','log10(err exact)','min err')
    xlabel('Iterations')
