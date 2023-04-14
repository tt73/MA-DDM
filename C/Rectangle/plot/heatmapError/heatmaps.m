N = 400; 

load_exact2
u2e = reshape(u_exact,N,N);
load_u2
u2 = reshape(u,N,N);
e2 = u2e - u2;

load_exact4
u4e = reshape(u_exact,N,N);
load_u4
u4 = reshape(u,N,N);
e4 = u4e - u4;

load_exact6
u6e = reshape(u_exact,N,N);
load_u6
u6 = reshape(u,N,N);
e6 = u6e - u6;

load_exact8
u8e = reshape(u_exact,N,N);
load_u8
u8 = reshape(u,N,N);
e8 = u8e - u8;

load_exact9
u9e = reshape(u_exact,N,N);
load_u9
u9 = reshape(u,N,N);
e9 = u9e - u9; 

%% 
figure 

subplot(151)
imagesc(abs(e2)) 
title('N_d=2')
colorbar('location','southoutside')

subplot(152)
imagesc(abs(e4))
title('N_d=4')
colorbar('location','southoutside')

subplot(153)
imagesc(abs(e6))
title('N_d=6')
colorbar('location','southoutside')

subplot(154) 
imagesc(abs(e8))
title('N_d=8')
colorbar('location','southoutside')

subplot(155) 
imagesc(abs(e9))
title('N_d=9')
colorbar('location','southoutside')


