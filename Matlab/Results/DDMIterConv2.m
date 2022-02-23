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
N = 2^6+1;
overlap = ceil((N-1)/20);
method = 3;
Outputs = cell(1,3);
tol = 1e-6;

examples = {'Radial Smooth','Degenate','Blow-Up' };
for choice = 1
    [Errors,Iterations,Times,Solutions,Domains,Parameters] = SchwarzSolver(N,overlap,choice,method,tol);   
    Outputs{choice} = {Errors,Iterations,Times,Solutions,Domains,Parameters};

    conv_iter = Iterations{1};
    ress = Errors{1};
    err_exact = Errors{2};
    min_err = Errors{4};
    depth = Parameters{2};

    figure, hold on, title(sprintf(['Error over Iteration with \n N = %d, depth = %d, overlap = %d, tol = 1e-6 \n Example = ' examples{choice}],N,depth,overlap),'FontSize',20)
    plot(1:conv_iter,log10(ress(1:conv_iter)),'g--','linewidth',2)
    plot(1:conv_iter,log10(err_exact(1:conv_iter)),'b:','linewidth',2)
    yline(log10(min_err),'k--','LineWidth',2)
    legend('log10(residue)','log10(err exact)','min err',"Location","best","FontSize",12)
    xlabel('Iterations','FontSize',18)
    
end
