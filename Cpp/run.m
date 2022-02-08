%% List problem parameters here

% Parameters needed to generate grid
x0 = -1; x1 = 1;
y0 = -1; y1 = 1;
N = 2^7+1;
h = (x1-x0)/(N+1);

% requirement: overlap + depth - 1 <= (N-1)/2
depth = ceil(h^(-1/3));
overlap = ceil((N-1)/20);
if (overlap + depth - 1 > (N-1)/2)
   error("overlap + depth exceeds mesh size")
end
% for around 10% overlap, choose overlap ~ (N-1)/20 

% choose F
choice = 1;
switch(choice)
   case 1
      DirBC = @(x,y) (exp((x.^2+y.^2)/2));
      contF = @(x,y) ((1+x.^2+y.^2).*exp(x.^2+y.^2));
   case 2
      pos = @(x) max(x,0);
      DirBC = @(x,y) (.5*pos(((x-.5).^2+(y-.5).^2).^.5-.2).^2);
      contF = @(x,y) (pos(1-.2./sqrt((x-.5).^2+(y-.5).^2)));
   case 3
      DirBC = @(x,y) (-sqrt(2-(x.^2+y.^2)));
      contF = @(x,y) (2*(2-(x.^2+y.^2)).^-2);
end


%% Call NewtonSchwarz 

% figure out how to send parameters to C++ program 

% 
                
