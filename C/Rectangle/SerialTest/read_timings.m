ksp_types="cg groppcg pipecg pipecgrr pipelcg pipeprcg pipecg2 nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls";
pc_types="none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat cp redistribute svd gamg kaczmarz hmg lmvm";
K_list = split(ksp_types);
P_list = split(pc_types);

nk = length(K_list);
np = length(P_list);

%% data to get 

iters  = zeros(np,nk);
errors = zeros(np,nk);
times  = zeros(np,nk); 

%% read file

fid = fopen('timing.out');
tline = fgetl(fid);
K_count = 1; 
P_count = 0;
while ischar(tline)
    tline = fgetl(fid);
    if(tline==-1)
        break
    end
    if(startsWith(tline,'ksp'))
        P_count = P_count+1; 
        tline = fgetl(fid);
        if(strcmp(tline,""))
            iters(P_count,K_count) = nan;
            errors(P_count,K_count) = nan;
            times(P_count,K_count) = nan;
        else
            iters(P_count,K_count) = str2double(tline);
            tline = fgetl(fid);
            errors(P_count,K_count) = str2double(tline);
            tline = fgetl(fid);
            times(P_count,K_count) = str2double(tline);
        end
    end
    if (P_count == np) 
        P_count = 0;
        K_count = K_count + 1;
    end
end

%% heatmap

figure 
heatmap(K_list,P_list,iters);
title("Iterations")

figure 
heatmap(K_list,P_list,errors);
title("Errors")

minerr = min(errors(:,1));


for k = 1:nk
    for p = 1:np
        if (iters(p,k) == 0) 
            times(p,k) = nan;
        elseif(errors(p,k)~=minerr) 
            times(p,k) = nan;
        end
    end
end
times(find(P_list == "mat"),find(K_list == "gcr")) = nan; % blank out gcr/mat 
times(find(P_list == "mg"),find(K_list == "richardson")) = nan; 


figure 
heatmap(K_list,P_list,times,'colormethod','median');
title("Times")

%% Best for each PC 

disp('Best pc for each ksp:')
for k = 1:nk 
    [mtime,ind] = min(times(:,k));
    fprintf('ksp = %10s, best pc = %12s, best time = %f\n',K_list(k),P_list(ind),mtime);
end

%% Notes 

% 


