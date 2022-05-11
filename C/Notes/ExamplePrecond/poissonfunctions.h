#ifndef POISSONFUNCTIONS_H_
#define POISSONFUNCTIONS_H_

/*
In 3D these functions approximate the residual of, and the Jacobian of, the
(slightly-generalized) Poisson equation
    - cx u_xx - cy u_yy - cz u_zz = f(x,y,z)
with Dirichlet boundary conditions  u = g(x,y,z)  on a domain
Omega = (0,Lx) x (0,Ly) x (0,Lz) discretized using a DMDA structured grid.
In 1D and 2D they do the same thing with the obvious reductions of dimension
in these equations.  (The domain is an interval, rectangle, or rectangular
solid.)  All of these function work with equally-spaced structured grids.  The
dimensions hx, hy, hz of the rectangular cells can have any positive values.

These functions promote code reuse and serve as canonical examples.  They
are used in ch6/fish.c, ch6/minimal.c, and ch12/obstacle.c.

The functions PoissonXDFunctionLocal(), X=1,2,3, compute residuals and are
designed as call-backs:

  ierr = DMDASNESSetFunctionLocal(dmda,INSERT_VALUES,
             (DMDASNESFunction)PoissonXDFunctionLocal,&user); CHKERRQ(ierr);

A rough estimate is made of the number of flops in the functions
PoissonXDFunctionLocal().

The PoissonXDJacobianLocal() functions are call-backs which assemble Jacobians
for the same problems:

  ierr = DMDASNESSetJacobianLocal(dmda,
             (DMDASNESJacobian)PoissonXDJacobianLocal,&user); CHKERRQ(ierr);

The matrices A are normalized so that if cells are square (h = hx = hy = hz)
then A / h^d approximates the Laplacian in d dimensions.  This is the way
the rows would be scaled in a Galerkin FEM scheme.  (The entries are O(1)
only if d=2.)  The Dirichlet boundary conditions are approximated using
diagonal Jacobian entries with the same values as the diagonal entries for
points in the interior.  Thus these Jacobian matrices have constant diagonal.
*/

// warning: the user is in charge of setting up ALL of this content!
//STARTDECLARE
typedef struct {
    // domain dimensions
    PetscReal Lx, Ly, Lz;
    // coefficients in  - cx u_xx - cy u_yy - cz u_zz = f
    PetscReal cx, cy, cz;
    // right-hand-side f(x,y,z)
    PetscReal (*f_rhs)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
    // Dirichlet boundary condition g(x,y,z)
    PetscReal (*g_bdry)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
    // additional context; see example usage in ch7/minimal.c
    void   *addctx;
} PoissonCtx;

PetscErrorCode Poisson1DFunctionLocal(DMDALocalInfo *info,
    PetscReal *au, PetscReal *aF, PoissonCtx *user);

PetscErrorCode Poisson2DFunctionLocal(DMDALocalInfo *info,
    PetscReal **au, PetscReal **aF, PoissonCtx *user);

PetscErrorCode Poisson3DFunctionLocal(DMDALocalInfo *info,
    PetscReal ***au, PetscReal ***aF, PoissonCtx *user);
//ENDDECLARE

/* This generates a tridiagonal sparse matrix.  If cx=1 then it has 2 on the
diagonal and -1 or zero in off-diagonal positions.  For example,
    ./fish -fsh_dim 1 -mat_view ::ascii_dense -da_refine N                */
PetscErrorCode Poisson1DJacobianLocal(DMDALocalInfo *info, PetscReal *au,
                                      Mat J, Mat Jpre, PoissonCtx *user);

/* If h = hx = hy and h = L/(m-1) then this generates a 2m-1 bandwidth
sparse matrix.  If cx=cy=1 then it has 4 on the diagonal and -1 or zero in
off-diagonal positions.  For example,
    ./fish -fsh_dim 2 -mat_view :foo.m:ascii_matlab -da_refine N
produces a matrix which can be read into Matlab/Octave.                   */
PetscErrorCode Poisson2DJacobianLocal(DMDALocalInfo *info, PetscReal **au,
                                      Mat J, Mat Jpre, PoissonCtx *user);

/* If h = hx = hy = hz and h = L/(m-1) then this generates a 2m^2-1 bandwidth
sparse matrix.  For example,
    ./fish -fsh_dim 3 -mat_view :foo.m:ascii_matlab -da_refine N          */
PetscErrorCode Poisson3DJacobianLocal(DMDALocalInfo *info, PetscReal ***au,
                                      Mat J, Mat Jpre, PoissonCtx *user);

/* The following function generates an initial iterate using either
  * zero
  * a random function (white noise; *no* smoothness)
In addition, one can initialize either using the boundary function g for
the boundary locations in the initial state, or not.                      */

typedef enum {ZEROS, RANDOM} InitialType;

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry,
                            Vec u, PoissonCtx *user);

#endif

