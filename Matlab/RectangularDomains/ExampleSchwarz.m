% This is a script that does runs schwarz alternating method with 2
% overlapping subdomains.  

addpath('Subroutines')

x0 = -1; x1 = 1; 
y0 = -1; y1 = 1; 
N = 2^3; 
h = (x1-x0)/N;
depth = 1; 
% Parameters needed to generate grid


