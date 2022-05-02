#ifndef MAFUNCTIONS_H_
#define MAFUNCTIONS_H_
/*

*/
typedef struct {
   // Stencil width
   PetscInt  width;
   // domain dimensions
   PetscReal Lx, Ly, Lz;
   // epsilon is the regularization term
   PetscReal epsilon;
   // quadrature weights
   PetscReal *weights;
   // Forward stencil directions in i and j (only 2D)
   PetscInt  *Si, *Sj;
   // right-hand-side f(x,y,z)
   PetscReal (*f_rhs)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
   // Dirichlet boundary condition g(x,y,z)
   PetscReal (*g_bdry)(PetscReal x, PetscReal y, PetscReal z, void *ctx);
   // Debugging flag 
   PetscBool debug; 
   // additional context; s
   void   *addctx;
} MACtx;

PetscErrorCode MA1DFunctionLocal(DMDALocalInfo *info, PetscReal *au, PetscReal *aF, MACtx *user);

PetscErrorCode MA2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, MACtx *user);

PetscErrorCode MA3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au, PetscReal ***aF, MACtx *user);

PetscErrorCode MA1DJacobianLocal(DMDALocalInfo *info, PetscReal *au, Mat J, Mat Jpre, MACtx *user);

PetscErrorCode MA2DJacobianLocal(DMDALocalInfo *info, PetscReal **au, Mat J, Mat Jpre, MACtx *user);

PetscErrorCode MA3DJacobianLocal(DMDALocalInfo *info, PetscReal ***au, Mat J, Mat Jpre, MACtx *user);

/* The following function generates an initial iterate using either
  * zero
  * a random function (white noise; *no* smoothness)
In addition, one can initialize either using the boundary function g for
the boundary locations in the initial state, or not.                      */
typedef enum {ZEROS, RANDOM, POISSON} InitialType;

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry, Vec u, MACtx *user);

PetscErrorCode ComputeWeights(PetscInt width, PetscInt order, MACtx *user);

PetscErrorCode ComputeFwdStencilDirs(PetscInt width, MACtx *user);

PetscErrorCode ComputeProjectionIndeces(PetscReal *di, PetscReal *dj, PetscInt i, PetscInt j, PetscInt Si, PetscInt Sj, PetscInt Nx, PetscInt Ny);

PetscErrorCode ApproxDetD2u(PetscReal *DetD2u, PetscInt dim, PetscReal *SDD, MACtx *user);

PetscErrorCode ComputeSDD(DMDALocalInfo *info, PetscReal **au, MACtx *user, PetscInt i, PetscInt j, PetscReal x, PetscReal y, PetscReal *SDD, PetscReal *hFwd, PetscReal *hBak);

PetscErrorCode PrintProjection(DM da, MACtx *user);

#endif