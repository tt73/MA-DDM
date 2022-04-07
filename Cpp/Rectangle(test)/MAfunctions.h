#ifndef MAFUNCTIONS_H_
#define MAFUNCTIONS_H_

/*
In 3D these functions approximate the residual of, and the Jacobian of, the
(slightly-generalized) MA equation
    det(D^2(u(x,y,z))) = f(x,y,z)
with Dirichlet boundary conditions  u = g(x,y,z)  on a domain
Omega = (0,Lx) x (0,Ly) x (0,Lz) discretized using a DMDA structured grid.
In 1D and 2D they do the same thing with the obvious reductions of dimension
in these equations.  (The domain is an interval, rectangle, or rectangular
solid.)  All of these function work with equally-spaced structured grids.  The
dimensions hx, hy, hz of the rectangular cells can have any positive values.

The functions MAXDFunctionLocal(), X=1,2,3, compute residuals and are
designed as call-backs:

  ierr = DMDASNESSetFunctionLocal(dmda,INSERT_VALUES,
             (DMDASNESFunction)MAXDFunctionLocal,&user); CHKERRQ(ierr);

A rough estimate is made of the number of flops in the functions
MAXDFunctionLocal().

The MAXDJacobianLocal() functions are call-backs which assemble Jacobians
for the same problems:

  ierr = DMDASNESSetJacobianLocal(dmda,
             (DMDASNESJacobian)MAXDJacobianLocal,&user); CHKERRQ(ierr);

The matrices A are normalized so that if cells are square (h = hx = hy = hz)
then A / h^d approximates the Laplacian in d dimensions.  This is the way
the rows would be scaled in a Galerkin FEM scheme.  (The entries are O(1)
only if d=2.)  The Dirichlet boundary conditions are approximated using
diagonal Jacobian entries with the same values as the diagonal entries for
points in the interior.  Thus these Jacobian matrices have constant diagonal.
*/
typedef struct {
   // domain dimensions
   PetscReal Lx, Ly, Lz;
   // epsilon is the regularization term
   PetscReal epsilon;
   // right-hand-side f(x,y,z)
   PetscReal (*f_rhs)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
   // Dirichlet boundary condition g(x,y,z)
   PetscReal (*g_bdry)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
   // additional context; see example usage in ch7/minimal.c
   void   *addctx;


} MACtx;

PetscErrorCode MA1DFunctionLocal(DMDALocalInfo *info,
    PetscReal *au, PetscReal *aF, MACtx *user);

PetscErrorCode MA2DFunctionLocal(DMDALocalInfo *info,
    PetscReal **au, PetscReal **aF, MACtx *user);

PetscErrorCode MA3DFunctionLocal(DMDALocalInfo *info,
    PetscReal ***au, PetscReal ***aF, MACtx *user);
//ENDDECLARE

PetscErrorCode MA1DJacobianLocal(DMDALocalInfo *info, PetscReal *au, Mat J, Mat Jpre, MACtx *user);

PetscErrorCode MA2DJacobianLocal(DMDALocalInfo *info, PetscReal **au, Mat J, Mat Jpre, MACtx *user);

PetscErrorCode MA3DJacobianLocal(DMDALocalInfo *info, PetscReal ***au, Mat J, Mat Jpre, MACtx *user);

/* The following function generates an initial iterate using either
  * zero
  * a random function (white noise; *no* smoothness)
In addition, one can initialize either using the boundary function g for
the boundary locations in the initial state, or not.                      */
typedef enum {ZEROS, RANDOM} InitialType;

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry, Vec u, MACtx *user);

#endif