% We want to solve on an N by N grid which has resolution h. The goal is to
% construct an initial guess by solving a coarse problem with resolution
% hc = k*h. The coarse grid is Nc by Nc, where Nc < N is to be determined.

% The inputs are N and k, and the outpud is Nc. 

% In petsc, there is a restriction on the interpolation between meshsizes.
% We must have (N-1)/(Nc-1) be an integer. That puts severe restrictions 



test_cases = [ 
  200, 4;    % 200-1 = 199 is a prime number 
  199, 4; 
    6, 4;
   16, 4; 
   16, 2;
   15, 3;
   31, 4;
   101, 4;
   4,  2; 
   6,  2; 
   7,  2;
];

for i = 1:length(test_cases)
    N = test_cases(i,1);
    k = test_cases(i,2);
    fprintf("\n\nTest %d: N = %d, k = %d\n",i,N,k)
    [b,Nc] = find_best(N,k); 
    fprintf("outut: best = %d, which leads to Nc = %d\n",b,Nc)
    if (Nc == N)
        fprintf("Interpolation not possible\n")
    elseif (b==k)
        fprintf("Exact coarse solve on 4h is possible\n")
    else
        fprintf("Ratio between h and coarse h_c is %f\n",(N+1)/(Nc-1))
        
    end
end


function [best,Nc] = find_best(N,k)

    frac = (N-1)/k;
    if (floor(frac)==frac) 
        Nc = frac + 1;
        best = k;
        return 
    end

    dist = 100;
    best = 0;
    Nc = N;
    for i = 1: k+15
        frac = (N-1)/i;
        if (floor(frac) == frac) 
            new_dist = abs(k-i);
            
            if (dist > new_dist)
                dist = new_dist;
                best = i;
                Nc = frac + 1;
           
            else 
%                 fprintf("Candidate i = %d with Nc = %d is too far\n",i,frac+1)
            end
  
        end
    end

end
