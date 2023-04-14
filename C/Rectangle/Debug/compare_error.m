
N = 200; 

h = 2/N;

load_exact4
load_u4

u4 = h*reshape(u,N,N)
u4_exact = reshape(u_exact,N,N)

load_exact8
load_u8

u8 = reshape(u,N,N)
u8_exact = reshape(u_exact,N,N)


figure 
subplot(121) 
surf(u4) 
subplot(122)
surf(u8)