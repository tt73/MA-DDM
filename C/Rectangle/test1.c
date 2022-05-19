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
   f[0] = PetscMax(1-0.2/PetscAbsReal(x-0.5),0);
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

static ExactFcnVec getuexact_ptr[3] = {&Form1DUExact, &Form2DUExact, &Form3DUExact};


typedef enum {ex10, ex11, ex12} ProblemType;
static const char* ProblemTypes[] = {"ex10","ex11","ex12","ProblemType","", NULL};

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

static const char* InitialTypes[] = {"zeros","random","corner","pyramid","InitialType","", NULL};


int main(int argc,char **args) {
   PetscErrorCode ierr;
   DM             da, da_after;
   SNES           snes;
   KSP            ksp;
   PC             pc;
   SNESLineSearch ls;
   DMDALocalInfo  info;
   Vec            u_initial,u,u_exact,err;
   MACtx          user; // see header file
   ExactFcnVec    getuexact;
   InitialType    initial;
   ProblemType    problem;           // manufactured problem using exp()
   PetscBool      debug,set_N,set_eps,set_width;
   PetscReal      h_eff,hx,hy,eps,errinf,normconst2h,err2h;
   char           gridstr[99];
   PetscInt       dim,width,N,Nx,Ny,order;

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;
   /*  Default params  - - - - - - - - - - - - - - - - - - - - - - - - - - -

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   N = Nx = Ny = 4; // # of interior points
   dim         = 2; // PDE space dimension
   width       = 1; // stencil width
   eps         = 0.25;  // epsilon = (hd)^2
   order       = 2; // guadrature order
   initial     = CORNER;
   problem     = ex10;
   user.Lx     = 1.0;
   user.Ly     = 1.0;
   user.Lz     = 1.0;
   debug       = PETSC_FALSE;
   // Get command args
   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"t1_", "options for test1.c", ""); CHKERRQ(ierr);
   ierr = PetscOptionsInt("-dim","dimension of problem (=1,2,3 only)","test1.c",dim,&dim,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-N","make the interior N by N","test1.c",N,&N,&set_N);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-Nx","number of interior nodes horizontally","test1.c",Nx,&Nx,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-Ny","number of interior nodes vertically","test1.c",Ny,&Ny,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-order","order of quadrature (default is 2)","test1.c",order,&order,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-debug","print out extra info","test1.c",debug,&debug,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-init_type","type of initial iterate","test1.c",InitialTypes,(PetscEnum)initial,(PetscEnum*)&initial,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Lx","set Lx in domain ([-Lx,Lx] x [-Ly,Ly] x [-Lz,Lz], etc.)","test1.c",user.Lx,&user.Lx,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Ly","set Ly in domain ([-Lx,Lx] x [-Ly,Ly] x [-Lz,Lz], etc.)","test1.c",user.Ly,&user.Ly,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-Lz","set Ly in domain ([-Lx,Lx] x [-Ly,Ly] x [-Lz,Lz], etc.)","test1.c",user.Lz,&user.Lz,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-problem","problem type; determines exact solution and RHS","test1.c",ProblemTypes,(PetscEnum)problem,(PetscEnum*)&problem,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsInt("-width","stencil width","test1.c",width,&width,&set_width);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-eps","regularization constant epsilon","test1.c",eps,&eps,&set_eps);CHKERRQ(ierr);
   ierr = PetscOptionsEnd(); CHKERRQ(ierr);
   /* Update dependent params
      It's optimal to choose depth = ceil(h^(-1/3)) where h is the stepsize.
      However, you can choose during runtime with "-t1_width <d>", where <d> is an integer.
      Epsilon is the regularization constant. We choose epsilon = (hd)^2 by default.
      We may change epsilon during runtime with "-t1_eps <f>", where <f> is a float
   */
   if (set_N) {
      Nx = Ny = N;
   }
   hx    = 2.0*user.Lx/(Nx+1);
   hy    = 2.0*user.Ly/(Ny+1);
   h_eff = PetscMax(hx,hy);
   if (!set_width) {
      width = PetscCeilReal(PetscPowReal(h_eff,-0.333));
   }
   h_eff *= width;
   if (!set_eps) {
      eps = h_eff*h_eff; // epsilon = (h*d)^2
   }
   user.debug   = debug;
   user.width   = width;
   user.epsilon = eps;
   user.g_bdry  = g_bdry_ptr[dim-1][problem];
   user.f_rhs   = f_rhs_ptr[dim-1][problem];
   /* DMDACreate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Set up a system of indexed nodes configured in a 1d, 2d, or 3d lattice.
      The dimensions are free to choose and the distribution of nodes is
      handled automatically.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   switch (dim) {
      case 1:
         ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Nx,1,1,NULL,&da); CHKERRQ(ierr);
         break;
      case 2:
         ierr = DMDACreate2d(PETSC_COMM_WORLD,   // MPI not important
                        DM_BOUNDARY_NONE,        // no periodicity in x
                        DM_BOUNDARY_NONE,        // no periodicity in y
                        width==1?DMDA_STENCIL_STAR:DMDA_STENCIL_BOX, // star means + sign stencil, box means stencil has diagonals
                        Nx,Ny,                     // mesh size in x & y directions
                        PETSC_DECIDE,PETSC_DECIDE, // local mesh size
                        1,                         // degree of freedom
                        width,                     // stencil width
                        NULL,NULL,                 // not important
                        &da); CHKERRQ(ierr);
         break;
      case 3:
         SETERRQ(PETSC_COMM_SELF,3,"dim = 3 not yet supported\n");
         ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,3,3,3,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,width,NULL,NULL,NULL,&da); CHKERRQ(ierr);
         break;
      default:
         SETERRQ(PETSC_COMM_SELF,1,"invalid dim for DMDA creation\n");
   }
   /* DM setup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      The code is set up to mesh a domain [-Lx,Lx] x [-Ly,Ly]
      with uniformly distributed points. Nx and Ny are the number
      of interior points. You can adjust the overlapping nodes
      with -da_overlap <int>.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = DMDASetOverlap(da,width,width,width);
   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMDASetUniformCoordinates(da,-user.Lx+hx,user.Lx-hx,-user.Ly+hy,user.Ly-hy,-user.Lz,user.Lz); CHKERRQ(ierr);
   ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
   /* SNES setup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      By default, we use Nonlinear Additive Schwarz method (NASM) for the
      nonlinear solver. On the local subdomains, we use one step of basic
      newton method. For the Jacobian solve, we use GMRES with multigrid
      preconditioning.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
   ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
   ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)(residual_ptr[dim-1]),&user); CHKERRQ(ierr);
   ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)(jacobian_ptr[dim-1]),&user); CHKERRQ(ierr);
   ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
   ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
   ierr = KSPSetType(ksp,KSPGMRES); CHKERRQ(ierr);
   // ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
   ierr = SNESSetType(snes,SNESNASM); CHKERRQ(ierr);
   ierr = SNESNASMSetType(snes,PC_ASM_RESTRICT);
   /* Linesearch - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = SNESGetLineSearch(snes,&ls);
   ierr = SNESLineSearchSetType(ls,SNESLINESEARCHBASIC);
   ierr = SNESLineSearchSetTolerances(ls,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1);
   /* Wide-stencil params - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Compute forward stencil directions for the determinant
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ComputeFwdStencilDirs(width,&user); // Stencil directions based on the L1 norm
   ComputeWeights(width,order,&user);  // Quadrature weights
   /* Solve - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Set up the inital guess on each sub-domain

   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = DMGetGlobalVector(da,&u_initial); CHKERRQ(ierr);
   ierr = InitialState(da,initial,u_initial,&user); CHKERRQ(ierr);
   ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
   ierr = SNESSolve(snes,NULL,u_initial); CHKERRQ(ierr);
   /* Error info - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Get the solution from DA.
      Get the exact solution with g.
      Compute the norm of the error.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = SNESGetDM(snes,&da_after); // let DM da_after hold the data management info (possibly different than original da)
   ierr = SNESGetSolution(snes,&u); // let Vec u hold the solution
   ierr = DMDAGetLocalInfo(da_after,&info); // retrieve local process info from DA
   ierr = DMCreateGlobalVector(da_after,&u_exact);
   getuexact = getuexact_ptr[dim-1];
   ierr = (*getuexact)(&info,u_exact,&user); // compute exact solution using g(x,y)
   ierr = VecDuplicate(u_exact,&err);        // allocated mem for `err`
   ierr = VecCopy(u_exact,err);              // make a copy of exact solution to `err`
   ierr = VecAYPX(err,-1.0,u); CHKERRQ(ierr);   // err <- u + (-1.0)*uexact
   ierr = VecNorm(err,NORM_INFINITY,&errinf); CHKERRQ(ierr);
   ierr = VecNorm(err,NORM_2,&err2h); CHKERRQ(ierr);
   /* Print Message - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Problem type, grid size, stencil width, epsilon, and error
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
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
   ierr = PetscPrintf(PETSC_COMM_WORLD, "problem %s on %s grid with d = %d, and eps = %.3f:\n"
               "  error |u-uexact|_inf = %.3e, |u-uexact|_h = %.3e\n",
               ProblemTypes[problem],gridstr,info.sw,user.epsilon,errinf,err2h); CHKERRQ(ierr);

   /* Debugging Info - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      To print out debugging info, add the option -t1_debug.
      The goal is to get info about the subdomains.
      Unfortunately, da and da_after lose local info by this point.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if (debug) {
      PetscInt MM, NN, mm, nn, dof, ss;
      DMDAGetInfo(da_after,&dim,&MM,&NN,NULL,&mm,&nn,NULL,&dof,&ss,NULL,NULL,NULL,NULL);
      PetscPrintf(PETSC_COMM_WORLD,"-- Grid is %d by %d, processors divided in %d by %d format.\n",MM,NN,mm,nn);
      PetscPrintf(PETSC_COMM_WORLD,"-- Dof = %d, Stencil = %d\n",dof,ss);
      PetscInt ox, oy;
      DMDAGetOverlap(da_after,&ox,&oy,NULL);
      PetscPrintf(PETSC_COMM_WORLD,"-- Overlap in x: %d, Overlap in y: %d\n",ox,oy);
      // PetscMPIInt    rank, size;
      // DMDALocalInfo  info;
      // MPI_Comm_size(PETSC_COMM_WORLD,&size);
      // MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      // for (int i=0; i<size; i++) {
      //    if (i==rank) {
      //       DMDAGetLocalInfo(da_after,&info);
      //       printf("-- Rank %d: x indeces [%d to %d], y indeces [%d to %d]\n",rank,info.xs,i<info.xs+info.xm,info.ys,i<info.ys+info.ym);
      //    }
      // }
   }
   /* MATLAB  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create MATLAB files load_u.m and load_exact.m which loads
      the numerical and exact solutions into a workspace.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   PetscViewer viewer;  // Viewer object fascilitates printing out solution
   ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); // initialize the viewer object
   ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_u.m",&viewer);  // set the file name
   ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for u
   ierr = PetscObjectSetName((PetscObject)u,"u");
   ierr = VecView(u,viewer);
   ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_exact.m",&viewer);  // set the file name
   ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for u
   ierr = PetscObjectSetName((PetscObject)u_exact,"u_exact");
   ierr = VecView(u_exact,viewer);
   // Free memory
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
   PetscReal  Lx, hx, x, *au,temp;
   Lx = user->Lx;
   hx = 2.0*Lx/(PetscReal)(info->mx-1);
   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (i=info->xs; i<info->xs+info->xm; i++) {
      x = -Lx + (i+1)*hx;
      user->g_bdry(x,0.0,0.0,user,&temp);
      au[i] = temp;
   }
   ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
   return 0;
}

PetscErrorCode Form2DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i, j;
   PetscReal  Lx, Ly, hx, hy, x, y, **au,temp;
   Lx = user->Lx; Ly = user->Ly;
   hx = 2.0*Lx/(PetscReal)(info->mx + 1);
   hy = 2.0*Ly/(PetscReal)(info->my + 1);
   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (j=info->ys; j<info->ys+info->ym; j++) {
      y = -Ly + (j+1)*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = -Lx + (i+1)*hx;
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