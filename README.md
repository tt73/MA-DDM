# MA-DDM
Domain decomposition method for the Monge-Ampere by <br />
> Tadanaga Takahashi - tt73@njit.edu <br />
> Jake Bursca - jb327@njit.edu <br />

## Introduction
This is a repository for the implementation of an experimental parallel domain decomposition method applied to the Monge-Ampere equation. The Monge-Ampere equation is a 2nd order fully non-linear "elliptic" partial differential equation. Literature on domain decomposition methods for non-linear problems is sparse, especially for fully non-linear equations. We implemented a Schwarz Alternating Procedure driven by a wide-stencil Newton's method for the local solve. The method is empirically stable for 3 test cases. One test tests for a degenerate solution. 

## Model 
In this project we are solving equations of the form <br />
![equation](https://bit.ly/3eJbnJs" align="center" border="0" alt="det(D^2u(x))=f(x)")
<br />

$ \det( D^2 u ) = f(x) $
![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

## Method 


## Todo List 
- [x] 1D Matlab Implementation 
- [x] 2D Matlab Jacobi Implementation on a rectangle
- [ ] 2D Matlad Krylov Implementation on a rectangle  
- [ ] 2D Solver in C++ on a rectangle with N vertical strips 
- [ ] 2D Solver in C++ for other domain splitting 
