/*
   This is a test code to solve the Monge-Ampere on a 2D rectangular domain
   with the stencil of width 1.

   This is the serial version, but most of it is ready to be
   done in parallel.
*/
static char help[] = "Create a wide-stencil grid in 2D.\n\n";

#include <petsc.h>
#include "MAfunctions.h"


// functions simply to put u_exact()=g_bdry() into a grid
// these are irritatingly-dimension-dependent inside ...
extern PetscErrorCode Form1DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode Form2DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode Form3DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode u_exact_1Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_2Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_3Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_1Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_2Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_3Dex10(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_1Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_2Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_3Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_1Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_2Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_3Dex11(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_1Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_2Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode u_exact_3Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_1Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_2Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);
extern PetscErrorCode f_rhs_3Dex12(PetscReal, PetscReal, PetscReal, void*, PetscReal*);

/*
   Solution    u(x) = exp(|x|^2/2), for x in Rn
   RHS: Det(D^2u(x)) = (1+|x|^2)*exp(n/2*|x|^2)
*/
PetscErrorCode u_exact_1Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal((x*x)/2.0);
   return 0;
}
PetscErrorCode u_exact_2Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal((x*x + y*y)/2.0);
   return 0;
}
PetscErrorCode u_exact_3Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal((x*x + y*y + z*z)/2.0);
   return 0;
}
// right-hand-side functions
PetscErrorCode f_rhs_1Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x)*PetscExpReal(x*x*0.5);
   return 0;
}
PetscErrorCode f_rhs_2Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x + y*y)*PetscExpReal(x*x + y*y);
   return 0;
}
PetscErrorCode f_rhs_3Dex10(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x + y*y + z*z)*PetscExpReal((x*x + y*y + z*z)*1.5);
   return 0;
}


/*
   Solution    u(x) = 1/2 max(|x-0.5| - 0.2,0)        , for x in Rn
   RHS: Det(D^2u(x)) = max(1-0.2/|x-0.5|, 0)
*/
PetscErrorCode u_exact_1Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscAbsReal(x-0.5)-0.2,0),2);
   return 0;
}
PetscErrorCode u_exact_2Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5))-0.2,0),2);
   return 0;
}
PetscErrorCode u_exact_3Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)+PetscSqr(z-0.5))-0.2,0),2);
   return 0;
}
// RHS
PetscErrorCode f_rhs_1Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) { // not sure if this is right
   if (x > 0.7 || x < 0.3){
      f[0] = 1.0;
   } else {
      f[0] = 0.0;
   }
   return 0;
}
PetscErrorCode f_rhs_2Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = PetscMax(1-0.2/PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)),0);
   return 0;
}
PetscErrorCode f_rhs_3Dex11(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = PetscMax(1-0.2/PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)),0);
   return 0;
}

/*
   Solution    u(x) = -sqrt(2 - |x|^2),              for x in Rn
   RHS: Det(D^2u(x)) = 2/(2- |x|^2)^p(n)
*/
PetscErrorCode u_exact_1Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = -PetscSqrtReal(2.0 - x*x);
   return 0;
}
PetscErrorCode u_exact_2Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = -PetscSqrtReal(2.0 - x*x - y*y);
   return 0;
}
PetscErrorCode u_exact_3Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = -PetscSqrtReal(2.0 - x*x - y*y - z*z);
   return 0;
}

// RHS
PetscErrorCode f_rhs_1Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2-x*x,1.5);
   return 0;
}
PetscErrorCode f_rhs_2Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2-x*x-y*y,2.0);
   return 0;
}
PetscErrorCode f_rhs_3Dex12(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2-x*x-y*y-z*z,2.5);
   return 0;
}


//STARTPTRARRAYS
// arrays of pointers to functions
static DMDASNESFunction residual_ptr[3]
    = {(DMDASNESFunction)&MA1DFunctionLocal,
       (DMDASNESFunction)&MA2DFunctionLocal,
       (DMDASNESFunction)&MA3DFunctionLocal};

static DMDASNESJacobian jacobian_ptr[3]
    = {(DMDASNESJacobian)&MA1DJacobianLocal,
       (DMDASNESJacobian)&MA2DJacobianLocal,
       (DMDASNESJacobian)&MA3DJacobianLocal};

typedef PetscErrorCode (*ExactFcnVec)(DMDALocalInfo*,Vec,MACtx*);

static ExactFcnVec getuexact_ptr[3]
    = {&Form1DUExact, &Form2DUExact, &Form3DUExact};
//ENDPTRARRAYS

typedef enum {ex10, ex11, ex12} ProblemType;
static const char* ProblemTypes[] = {"ex10","ex11","ex12",
                                     "ProblemType", "", NULL};

// more arrays of pointers to functions:   ..._ptr[DIMS][PROBLEMS]
typedef PetscErrorCode (*PointwiseFcn)(PetscReal,PetscReal,PetscReal,void*,PetscReal*);

static PointwiseFcn g_bdry_ptr[3][3]
    = {{&u_exact_1Dex10, &u_exact_1Dex11, &u_exact_1Dex12},
       {&u_exact_2Dex10, &u_exact_2Dex11, &u_exact_2Dex12},
       {&u_exact_3Dex10, &u_exact_3Dex11, &u_exact_3Dex12}};

static PointwiseFcn f_rhs_ptr[3][3]
    = {{&f_rhs_1Dex10, &f_rhs_1Dex11, &f_rhs_1Dex12},
       {&f_rhs_2Dex10, &f_rhs_2Dex11, &f_rhs_2Dex12},
       {&f_rhs_3Dex10, &f_rhs_3Dex11, &f_rhs_3Dex12}};

static const char* InitialTypes[] = {"zeros","random",
                                     "InitialType", "", NULL};


int main(int argc,char **args) {
   PetscErrorCode ierr;
   DM             da, da_after;
   SNES           snes;
   KSP            ksp;
   Vec            u_initial, u, u_exact, err;
   MACtx          user;
   DMDALocalInfo  info;
   PetscReal      errinf, normconst2h, err2h;
   char           gridstr[99];
   ExactFcnVec    getuexact;

   // fish defaults:
   PetscInt       dim = 2;                  // 2D
   ProblemType    problem = ex10;        // manufactured problem using exp()
   InitialType    initial = ZEROS;          // set u=0 for initial iterate
   PetscBool      debug = PETSC_FALSE,gonboundary = PETSC_TRUE; // initial iterate has u=g on boundary

   PetscInt s = 1;  // Stencil width
   PetscInt N = 2;  // We want an interior domain that is N by N

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;

   // get options and configure context
   user.epsilon = 1.0/((N+1)*(N+1));
   user.Lx = 1.0;
   user.Ly = 1.0;
   user.Lz = 1.0;
   // command args for this programa need a 't1_' prefix
   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"t1_", "options for test1.c", ""); CHKERRQ(ierr);
   ierr = PetscOptionsInt("-dim","dimension of problem (=1,2,3 only)","test1.c",dim,&dim,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-initial_gonboundary","set initial iterate to have correct boundary values",
      "test1.c",gonboundary,&gonboundary,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-debug","print out extra info","test1.c",debug,&debug,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-initial_type","type of initial iterate",
      "test1.c",InitialTypes,(PetscEnum)initial,(PetscEnum*)&initial,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Lx", "set Lx in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
      "test1.c",user.Lx,&user.Lx,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Ly","set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
      "test1.c",user.Ly,&user.Ly,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Lz","set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
      "test1.c",user.Lz,&user.Lz,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-problem","problem type; determines exact solution and RHS",
      "test1.c",ProblemTypes,(PetscEnum)problem,(PetscEnum*)&problem,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsEnd(); CHKERRQ(ierr);

   user.g_bdry = g_bdry_ptr[dim-1][problem];
   user.f_rhs = f_rhs_ptr[dim-1][problem];
   user.debug = debug;

   /*
      Only doing dim = 1 and 2.
      Since this is a serial code, the only important information is:
       1. A N+2 by N+2 square mesh is created.
       2. The information is saved to `da`

      Note: The reason why we do N+2 is because we want two extra layers
            of structured nodes for the boundary leaving the desired
            N by N interior nodes. So h = 1/(N+1).
   */
   switch (dim) {
      case 1:
         ierr = DMDACreate1d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE,N+2,1,1, NULL, &da); CHKERRQ(ierr);
         break;
      case 2:
         ierr = DMDACreate2d(PETSC_COMM_WORLD,   // MPI not important
                           DM_BOUNDARY_NONE,  // periodicity in  x
                           DM_BOUNDARY_NONE,  // periodicity in y
                           DMDA_STENCIL_STAR, // stencil type: box vs star
                           N+2,N+2,           // mesh size in x & y directions
                           PETSC_DECIDE,PETSC_DECIDE, // local mesh size
                           1,                 // degree of freedom
                           s,                 // stencil width
                           NULL,NULL,         // not important
                           &da); CHKERRQ(ierr);
         break;
      case 3:
         ierr = DMDACreate3d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_STAR,
               3,3,3,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
               1,1,NULL,NULL,NULL,&da); CHKERRQ(ierr);
         break;
      default:
         SETERRQ(PETSC_COMM_SELF,1,"invalid dim for DMDA creation\n");
   }

   /*
      The line below allows you to make changes to the default settings on
      the DM object. For this code, the only important ones are:
       -da_grid_x <m>, which changes the number of cols of points to m
       -da_grid_y <n>, which changes the number of rows of points to n
       -da_refine <s>, which multiplies the default N+2 x N+2 grid by a factor of s+1
   */
   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMDASetUniformCoordinates(da,0.0,user.Lx,0.0,user.Ly,0.0,user.Lz); CHKERRQ(ierr);
   ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);

   // set SNES call-backs
   /*
      SetFuntionLocal - link the custom residue function to `da`
      SetJacobianLocal - link the custom Jacobian function to `da`
      I am not going to write an explicit Jacobian.
   */
   ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
   ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
   ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)(residual_ptr[dim-1]),&user); CHKERRQ(ierr);
   ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)(jacobian_ptr[dim-1]),&user); CHKERRQ(ierr);

   ierr = SNESSetType(snes,SNESNASM); CHKERRQ(ierr);
   ierr = SNESNASMSetType(snes,PC_ASM_RESTRICT);
   ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,500,PETSC_DEFAULT);
   ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
   ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
   ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

   // set initial iterate and then solve
   ierr = DMGetGlobalVector(da,&u_initial); CHKERRQ(ierr);
   ierr = InitialState(da, initial, gonboundary, u_initial, &user); CHKERRQ(ierr);
   ierr = SNESSolve(snes,NULL,u_initial); CHKERRQ(ierr);

   // Get the numerical and exact solution as Vecs
   ierr = SNESGetDM(snes,&da_after); // let DM da_after hold the data management info (possibly different than original da)
   ierr = SNESGetSolution(snes,&u); // let Vec u hold the solution
   ierr = DMDAGetLocalInfo(da_after,&info); // retrieve local process info from DA
   ierr = DMCreateGlobalVector(da_after,&u_exact);
   getuexact = getuexact_ptr[dim-1];
   ierr = (*getuexact)(&info,u_exact,&user); // compute exact solution using g(x,y)
   ierr = VecDuplicate(u_exact,&err);        // allocated mem for `err`
   ierr = VecCopy(u_exact,err);              // make a copy of exact solution to `err`

   // Compute the error
   ierr = VecAYPX(err,-1.0,u); CHKERRQ(ierr);   // err <- u + (-1.0)*uexact
   ierr = VecNorm(err,NORM_INFINITY,&errinf); CHKERRQ(ierr);
   ierr = VecNorm(err,NORM_2,&err2h); CHKERRQ(ierr);

   switch (dim) {
      case 1:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1));
         snprintf(gridstr,99,"%d point 1D",info.mx);
         break;
      case 2:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1)*(info.my-1));
         snprintf(gridstr,99,"%d x %d point 2D",info.mx,info.my);
         break;
      case 3:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1)*(info.my-1)*(info.mz-1));
         snprintf(gridstr,99,"%d x %d x %d point 3D",info.mx,info.my,info.mz);
         break;
      default:
         SETERRQ(PETSC_COMM_SELF,4,"invalid dim value in final report\n");
   }
   err2h /= normconst2h; // like continuous L2
   ierr = PetscPrintf(PETSC_COMM_WORLD, "problem %s on %s grid:\n"
               "  error |u-uexact|_inf = %.3e, |u-uexact|_h = %.3e\n",
               ProblemTypes[problem],gridstr,errinf,err2h); CHKERRQ(ierr);

   // Print out the solution to another file
   PetscViewer viewer;  // Viewer object fascilitates printing out solution
   PetscViewerCreate(PETSC_COMM_WORLD, &viewer); // initialize the viewer object
   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_u.m",&viewer);  // set the file name
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for u
   PetscObjectSetName((PetscObject)u,"u");
   VecView(u,viewer);

   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_exact.m",&viewer);  // set the file name
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for u
   PetscObjectSetName((PetscObject)u_exact,"u_exact");
   VecView(u_exact,viewer);

   ierr = DMDestroy(&da); CHKERRQ(ierr);
   ierr = VecDestroy(&err); CHKERRQ(ierr);
   ierr = VecDestroy(&u_exact); CHKERRQ(ierr);
   ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
   ierr = SNESDestroy(&snes); CHKERRQ(ierr);
   return PetscFinalize();
}


PetscErrorCode Form1DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i;
   PetscReal  xmax[1], xmin[1], hx, x, *au,temp;
   ierr = DMGetBoundingBox(info->da,xmin,xmax); CHKERRQ(ierr);
   hx = (xmax[0] - xmin[0]) / (info->mx - 1);
   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (i=info->xs; i<info->xs+info->xm; i++) {
      x = xmin[0] + i*hx;
      user->g_bdry(x,0.0,0.0,user,&temp);
      au[i] = temp;
   }
   ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
   return 0;
}

PetscErrorCode Form2DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i, j;
   PetscReal  xymin[2], xymax[2], hx, hy, x, y, **au,temp;
   ierr = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
   hx = (xymax[0] - xymin[0]) / (info->mx - 1);
   hy = (xymax[1] - xymin[1]) / (info->my - 1); // mx = my = N+2

   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (j=info->ys; j<info->ys+info->ym; j++) {
      y = xymin[1] + j*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = xymin[0] + i*hx;
         user->g_bdry(x,y,0.0,user,&temp);
         au[j][i] = temp;
      }
   }
   ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
   return 0;
}

PetscErrorCode Form3DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt  i, j, k;
   PetscReal xyzmin[3], xyzmax[3], hx, hy, hz, x, y, z, ***au, temp;
   ierr = DMGetBoundingBox(info->da,xyzmin,xyzmax); CHKERRQ(ierr);
   hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
   hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
   hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (k=info->zs; k<info->zs+info->zm; k++) {
      z = xyzmin[2] + k*hz;
      for (j=info->ys; j<info->ys+info->ym; j++) {
         y = xyzmin[1] + j*hy;
         for (i=info->xs; i<info->xs+info->xm; i++) {
            x = xyzmin[0] + i*hx;
            user->g_bdry(x,y,z,user,&temp);
            au[k][j][i] = temp;
         }
      }
   }
   ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
   return 0;
}