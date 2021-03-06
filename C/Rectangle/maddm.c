/*
   This is a code to solve the Monge-Ampere on a 2D rectangular domain
   discretized into a structured grid. It is set to solve the problem
   with nonlinear additive schwarz + Newton's method.

   Compile this code with the makefile. Run with `./maddm`.
   You can type  `./maddm -help | grep maddm` to see the options for this code.
*/
static char help[] = "Solve the Monge-Ampere equation for one of three example problems in 1D or 2D using a wide-stencil discretization scheme on a rectangular grid.\n\n";

#include <petsc.h>
#include "MAfunctions.h"

extern PetscErrorCode Form1DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode Form2DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode Form3DUExact(DMDALocalInfo*, Vec, MACtx*);
extern PetscErrorCode ComputeRHS(DMDALocalInfo*, Vec, MACtx*);

/* Problem 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Smooth example.
   Solution    u(x) = exp(|x|^2/2), for x in Rn
   RHS: Det(D^2u(x)) = (1+|x|^2)*exp(n/2*|x|^2)
   Defualt domain: [-1,-1]^2
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode u_exact_1D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal(x*x/2.0);
   return 0;
}
PetscErrorCode u_exact_2D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal((x*x + y*y)/2.0);
   return 0;
}
PetscErrorCode u_exact_3D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscExpReal((x*x + y*y + z*z)/2.0);
   return 0;
}
PetscErrorCode f_rhs_1D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x)*PetscExpReal(x*x/2.0);
   return 0;
}
PetscErrorCode f_rhs_2D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x + y*y)*PetscExpReal(x*x + y*y);
   return 0;
}
PetscErrorCode f_rhs_3D_ex1(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (1.0 + x*x + y*y + z*z)*PetscExpReal((x*x + y*y + z*z)*1.5);
   return 0;
}

/* Problem 2 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Fully degenerate example.
   Solution    u(x) = 1/2 max(|x-0.5| - 0.2,0)^2        , for x in Rn
   RHS: Det(D^2u(x)) = max(1-0.2/|x-0.5|, 0)
   Defualt domain: [0,-1]^2
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode u_exact_1D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscAbsReal(x-0.5)-0.2,0),2);
   return 0;
}
PetscErrorCode u_exact_2D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5))-0.2,0),2);
   return 0;
}
PetscErrorCode u_exact_3D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.5*PetscPowReal(PetscMax(PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)+PetscSqr(z-0.5))-0.2,0),2);
   return 0;
}
PetscErrorCode f_rhs_1D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = (x<=0.3 || x>=0.7) ? 1.0 : 0.0;
   return 0;
}
PetscErrorCode f_rhs_2D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = PetscMax(1-0.2/PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)),0);
   return 0;
}
PetscErrorCode f_rhs_3D_ex2(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = PetscMax(1-0.2/PetscSqrtReal(PetscSqr(x-0.5)+PetscSqr(y-0.5)),0);
   return 0;
}

/* Problem 3 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Gradient blow-up on the boundary.
   Solution    u(x) = -sqrt(2 - |x|^2),              for x in Rn
   RHS: Det(D^2u(x)) = 2/(2- |x|^2)^p(n)
   Defualt domain: [0,-1]^2
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode u_exact_1D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = -PetscSqrtReal(2.0 - x*x);
   return 0;
}
PetscErrorCode u_exact_2D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   PetscReal temp = -PetscSqrtReal(2.0 - x*x - y*y);
   u[0] = PetscIsInfOrNanScalar(temp) ? 0.0 : temp;
   return 0;
}
PetscErrorCode u_exact_3D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = -PetscSqrtReal(2.0 - x*x - y*y - z*z);
   return 0;
}
PetscErrorCode f_rhs_1D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2.0-x*x,1.5);
   return 0;
}
PetscErrorCode f_rhs_2D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2.0-x*x-y*y,2.0);
   return 0;
}
PetscErrorCode f_rhs_3D_ex3(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 2.0/PetscPowReal(2.0-x*x-y*y-z*z,2.5);
   return 0;
}

/* Problem 4 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Semi-degenerate polynomial. ?? = (0.6,0.4)
   Solution    u(x) = (????x)^2,              for x in R2
   RHS: 0
   Defualt domain: [-1,-1]^2
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode u_exact_1D_ex4(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = 0.6*0.6*x*x;
   return 0;
}
PetscErrorCode u_exact_2D_ex4(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscPowReal(0.6*x + 0.4*y,2.0);
   return 0;
}
PetscErrorCode u_exact_3D_ex4(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * u) {
   u[0] = PetscPowReal(0.6*x + 0.4*y + z, 2.0);
   return 0;
}
PetscErrorCode f_rhs_ex4(PetscReal x, PetscReal y, PetscReal z, void *ctx, PetscReal * f) {
   f[0] = 0;
   return 0;
}

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

typedef enum {ex1, ex2, ex3, ex4} ProblemType;
static const char* ProblemTypes[] = {"ex1","ex2","ex3","ex4","ProblemType","", NULL};

// more arrays of pointers to functions:   ..._ptr[DIMS][PROBLEMS]
typedef PetscErrorCode (*PointwiseFcn)(PetscReal,PetscReal,PetscReal,void*,PetscReal*);

static PointwiseFcn g_bdry_ptr[3][4]
    = {{&u_exact_1D_ex1, &u_exact_1D_ex2, &u_exact_1D_ex3, &u_exact_1D_ex4},
       {&u_exact_2D_ex1, &u_exact_2D_ex2, &u_exact_2D_ex3, &u_exact_2D_ex4},
       {&u_exact_3D_ex1, &u_exact_3D_ex2, &u_exact_3D_ex3, &u_exact_3D_ex3}};

static PointwiseFcn f_rhs_ptr[3][4]
    = {{&f_rhs_1D_ex1, &f_rhs_1D_ex2, &f_rhs_1D_ex3, &f_rhs_ex4},
       {&f_rhs_2D_ex1, &f_rhs_2D_ex2, &f_rhs_2D_ex3, &f_rhs_ex4},
       {&f_rhs_3D_ex1, &f_rhs_3D_ex2, &f_rhs_3D_ex3, &f_rhs_ex4}};

static const char* InitialTypes[] = {"zeros","random","corner","pyramid","InitialType","", NULL};


int main(int argc,char **args) {
   PetscErrorCode ierr;
   DM             da,da_after;
   SNES           snes,subsnes;
   SNESLineSearch subls;
   KSP            subksp;
   PC             subpc;
   DMDALocalInfo  info;
   Vec            u_initial,u,u_exact,err;
   MACtx          user; // monge-ampere user context. see header file
   ExactFcnVec    getuexact;
   InitialType    initial;
   ProblemType    problem;
   PetscBool      debug,set_N,set_eps,set_width,printSol,fast,mixed,htn,sin,aspin;
   PetscReal      h_eff,hx,hy,eps,errinf,normconst2h,err2h,op;
   char           gridstr[99];
   PetscInt       dim,width,N,Nx,Ny,order,its;
   PetscLogDouble t1,t2;
   PetscMPIInt    rank,size;

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;
   /*  Get parameters  - - - - - - - - - - - - - - - - - - - - - - - - - - -
      The defualt parameters are initialized below.
      They can all be changed during runtime.
      By default, we solve ex1 on [-1,1]^2.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   N = Nx = Ny = 4; // # of interior points
   dim         = 2; // PDE space dimension
   width       = 1; // stencil width
   eps         = 0.25;  // epsilon = (hd)^2
   order       = 2;     // guadrature order
   op          = 0.1;   // domain overlap percentage
   initial     = ZEROS; // initial guess for outer iterative method
   problem     = ex1;   // choose ex1, ex2, or ex3
   user.xmin   = -1.0; user.xmax = 1.0; // x limits [-1, 1]
   user.ymin   = -1.0; user.ymax = 1.0; // y limits [-1, 1]
   user.zmin   = -1.0; user.zmax = 1.0; // z limits [-1, 1]
   fast        = PETSC_FALSE; // fast might be unstable so leave it off by default
   debug       = PETSC_FALSE; // option to print out extra info
   printSol    = PETSC_FALSE; // option to generate MATLAB solution scripts
   mixed       = PETSC_FALSE; // option to use different solver on each subdomain
   htn         = PETSC_FALSE; // option to use high-tolerance Newton
   sin         = PETSC_FALSE; // option to use single-iteration Newton
   aspin       = PETSC_FALSE; // option to use additive schwarz preconditioned inexact newton

   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","options for maddm.c",""); CHKERRQ(ierr);
   ierr = PetscOptionsInt("-dim","dimension of problem (=1,2,3 only)","maddm.c",dim,&dim,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-N","make the interior N by N","maddm.c",N,&N,&set_N);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-Nx","number of interior nodes horizontally","maddm.c",Nx,&Nx,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-Ny","number of interior nodes vertically","maddm.c",Ny,&Ny,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-order","order of quadrature (default is 2)","maddm.c",order,&order,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsInt("-width","stencil width for MA discretization","maddm.c",width,&width,&set_width);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-fast","use settings that are tuned to under-solve the local problems","maddm.c",fast,&fast,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-debug","print out extra info","maddm.c",debug,&debug,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-sol","generate MATLAB solution files","maddm.c",printSol,&printSol,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-mixed","sub-index the local domains and use a different solver on the last subdomain","maddm.c",mixed,&mixed,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-htn","option to use high-tolerance Newton (HTN)","maddm.c",htn,&htn,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-sin","option to use single-iteration Newton (HTN)","maddm.c",sin,&sin,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsBool("-aspin","option to use additive schwarz preconditioned inexact newton (ASPIN)","maddm.c",aspin,&aspin,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-init_type","type of initial iterate","maddm.c",InitialTypes,(PetscEnum)initial,(PetscEnum*)&initial,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsReal("-eps","regularization constant epsilon","maddm.c",eps,&eps,&set_eps);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-op","domain overlap percentage (0.0 to 1.0)","maddm.c",op,&op,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnum("-problem","problem type; (ex1, ex2, ex3, ex4)","maddm.c",ProblemTypes,(PetscEnum)problem,(PetscEnum*)&problem,NULL); CHKERRQ(ierr);
   if (problem==ex2 || problem==ex3) { // default limits for ex2 & ex3
      user.xmin   = 0.0; user.xmax = 1.0; // x limits [0, 1]
      user.ymin   = 0.0; user.ymax = 1.0; // y limits [0, 1]
      user.zmin   = 0.0; user.zmax = 1.0; // z limits [0, 1]
   }
   ierr = PetscOptionsReal("-xmin","set limit of domain ([xmin,1] x [-1,1] x [-1,1])","maddm.c",user.xmin,&user.xmin,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-xmax","set limit of domain ([-1,xmax] x [-1,1] x [-1,1])","maddm.c",user.xmax,&user.xmax,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-ymin","set limit of domain ([-1,1] x [ymin,1] x [-1,1])","maddm.c",user.ymin,&user.ymin,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-ymax","set limit of domain ([-1,1] x [-1,ymax] x [-1,1])","maddm.c",user.ymax,&user.ymax,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-zmin","set limit of domain ([-1,1] x [-1,1] x [zmin,1])","maddm.c",user.zmin,&user.zmin,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsReal("-zmax","set limit of domain ([-1,1] x [-1,1] x [-1,zmax])","maddm.c",user.zmax,&user.zmax,NULL);CHKERRQ(ierr);
   ierr = PetscOptionsEnd(); CHKERRQ(ierr);
   /* Update dependent params
      It's optimal to choose depth = ceil(h^(-1/3)) where h is the stepsize.
      However, you can choose during runtime with "-t1_width <d>", where <d> is an integer.
      Epsilon is the regularization constant. We choose epsilon = (hd)^2 by default.
      We may change epsilon during runtime with "-t1_eps <f>", where <f> is a float
   */
   if (set_N)
      Nx = Ny = N;
   hx    = (user.xmax - user.xmin)/(Nx+1);
   hy    = (user.ymax - user.ymin)/(Ny+1);
   h_eff = PetscMax(hx,hy);
   if (!set_width)
      width = PetscCeilReal(PetscPowReal(h_eff,-0.333));
   h_eff *= width;
   if (!set_eps)
      eps = h_eff*h_eff; // epsilon = (h*d)^2
   user.debug   = debug;
   user.width   = width;
   user.epsilon = eps;
   user.g_bdry  = g_bdry_ptr[dim-1][problem];
   user.f_rhs   = f_rhs_ptr[dim-1][problem];
   /* DMDACreate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Set up a system of indexed nodes configured in a 1d, 2d, or 3d lattice.
      The dimensions are free to choose and the distribution of nodes of the
      grid into subdomains is handled automatically, hence the `PETSC_DECIDE`.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   switch (dim) {
      case 1:
         ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Nx,1,1,NULL,&da); CHKERRQ(ierr);
         break;
      case 2:
         ierr = DMDACreate2d(PETSC_COMM_WORLD,     //
                        DM_BOUNDARY_NONE,          // no periodicity in x
                        DM_BOUNDARY_NONE,          // no periodicity in y
                        width==1?DMDA_STENCIL_STAR:DMDA_STENCIL_BOX, // star = cardinal directions, box = more general
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
      The code is set up to mesh a domain [xmin,xmax] x [ymin,ymax]
      with uniformly distributed points. Nx and Ny are the number
      of interior points.
      You can adjust the overlapping nodes directly with -da_overlap <int>
      or adjust the percentage with -op <float> where float ranges from 0 to 1.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   PetscInt mm, nn, olx, oly;
   ierr = DMSetUp(da); CHKERRQ(ierr); // initialize the grid distribution, then get mm and nn
   ierr = DMDAGetInfo(da,NULL,NULL,NULL,NULL,&mm,&nn,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);
   olx = PetscCeilReal(op*Nx/(PetscReal)mm); // x overlap
   oly = PetscCeilReal(op*Ny/(PetscReal)nn); // y overlap
   ierr = DMDASetOverlap(da,olx,oly,olx);
   ierr = DMSetFromOptions(da); CHKERRQ(ierr); // the nodes
   ierr = DMDASetUniformCoordinates(da,user.xmin+hx,user.xmax-hx,user.ymin+hy,user.ymax-hy,user.zmin,user.zmax); CHKERRQ(ierr);
   ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
   /* SNES setup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      By default, we use Nonlinear Additive Schwarz method (NASM) for the
      nonlinear solver. On the local subdomains, we use one step of basic
      newton method. For the Jacobian solve, we use GMRES with a SSOR
      preconditioner.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
   ierr = SNESSetType(snes,SNESNASM); CHKERRQ(ierr);
   ierr = SNESNASMSetType(snes,PC_ASM_RESTRICT); CHKERRQ(ierr);
   ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
   ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)(residual_ptr[dim-1]),&user); CHKERRQ(ierr);
   ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)(jacobian_ptr[dim-1]),&user); CHKERRQ(ierr);
   /* Subdomain Solve - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      The following lines of code extracts the SNES, LS, KSP, and PC
      objects from the local subsolver given to the current rank.

      By default, we set the local solver to use the deflated GMRES and use an
      efficient SSOR preconditioner. This was experimentally found to be a
      fast and reliable linear solver combo.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = SNESSetUp(snes); CHKERRQ(ierr); // initialize subdomains
   MPI_Comm_size(PETSC_COMM_WORLD,&size); // get # of processors
   MPI_Comm_rank(PETSC_COMM_WORLD,&rank); // get index of current procs
   SNESNASMGetSNES(snes,0,&subsnes);      // get local SNES
   SNESGetLineSearch(subsnes,&subls);     // get local linesearch
   SNESGetKSP(subsnes,&subksp);           // get local KSP
   KSPGetPC(subksp,&subpc);               // get local PC
   KSPSetType(subksp,KSPDGMRES);  // rtol = 1e-5 by default
   PCSetType(subpc,PCEISENSTAT);  // fast accurate linear solver combo

   /* Serial Code - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      If the code is run serially, there is only 1 subdomain.
      In this case, NASM is equivalent to Newton's method.
      The default SNES iterations is 50. We increase it to avoid needless restarting.
      You need to add -sub_snes_monitor to view the residues of the method.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if (size==1) { // only 1 subdomain - increase defualt iterations which is only 50
      SNESSetTolerances(subsnes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,10000,PETSC_DEFAULT);
   }

   /* Tuning Local the Solve - - - - - - - - - - - - - - - - - - - - - - - - -
      In order to get speedup with NASM, we need to under-solve the local
      Newton, linesearch, and GMRES options.

      There's two recommended options: the high-tolerance Newton (HTN) and
      the single-iteration Newton (SIN). The HTN changes the Newton tolerance
      and the GMRES tolerance both to 1e-1. The SIN caps the max number of
      Newton iterations to one and sets the GMRES tolerance to 1e-2.

      For both methods, we set the linesearch to 2nd order BT. The 3rd order
      BT does a while loop to meet a convergence criterion which can be costly.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if (sin) {
      SNESSetTolerances(subsnes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
      KSPSetTolerances(subksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
      SNESLineSearchSetType(subls,SNESLINESEARCHBT);
      SNESLineSearchSetOrder(subls,2);
   } else if(htn) {
      SNESSetTolerances(subsnes,PETSC_DEFAULT,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
      KSPSetTolerances(subksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
      SNESLineSearchSetType(subls,SNESLINESEARCHBT);
      SNESLineSearchSetOrder(subls,2);
   }


   // } else if (rank<size-1) { // default settings
   //    // newton rres < 0.01
   //    // gmres rres < 0.01
   //    // ls rres < 0.01
   //    SNESSetTolerances(subsnes,PETSC_DEFAULT,1.e-2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //    KSPSetTolerances(subksp,1.e-2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //    SNESLineSearchSetTolerances(subls,PETSC_DEFAULT,PETSC_DEFAULT,1.e-2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   // }


   // if ((rank==size-1) && mixed) { // final domain
   //    // SNESSetType(subsnes,SNESFAS); CHKERRQ(ierr);
   //    SNESSetType(subsnes,SNESNEWTONLS);
   //    if (fast){
   //       // do 1 iteration
   //       SNESSetTolerances(subsnes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
   //    } else {
   //       SNESSetTolerances(subsnes,PETSC_DEFAULT,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //    }
   //    KSPSetType(subksp,KSPGMRES);
   //    PCSetType(subpc,PCEISENSTAT);
   // } else  {
   //    if (fast) {
   //       // 1 iteration only + L2 linesearch
   //       SNESLineSearchSetType(subls,SNESLINESEARCHL2); // secant L2 linesearch
   //       SNESSetTolerances(subsnes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,1,PETSC_DEFAULT);
   //       KSPSetTolerances(subksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //       SNESLineSearchSetTolerances(subls,PETSC_DEFAULT,PETSC_DEFAULT,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //    } else if (size > 1) {
   //       // newton rres < 0.1
   //       // gmres res < 0.1
   //       SNESSetTolerances(subsnes,PETSC_DEFAULT,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //       KSPSetTolerances(subksp,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //       SNESLineSearchSetTolerances(subls,PETSC_DEFAULT,PETSC_DEFAULT,1.e-1,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
   //    }
   //    KSPSetType(subksp,KSPGMRES);
   //    PCSetType(subpc,PCEISENSTAT);
   // }

   /* ASPIN - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Additive Schwarz Preconditioned Inexact Newton is a completely different
      from NASM. It's using Newton as the global method and using NASM as
      a preconditioner. The default settings are terrible so I made these
      options available.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if(aspin) {
      ierr = SNESSetType(snes,SNESASPIN); CHKERRQ(ierr);
      ierr = SNESGetNPC(snes,&subsnes);
      // ierr = SNESGetLineSearch(subsnes,&subls); // get local linesearch, default is cubic BT
      // ierr = SNESGetKSP(subsnes,&subksp);       // get local KSP, default is preonly
      // ierr = KSPGetPC(subksp,&subpc);           // get local PC, default is lu
      // ierr = KSPSetType(subksp,KSPDGMRES);      // change KSP to DGMRES
      // ierr = PCSetType(subpc,PCEISENSTAT);      // change PC to eisenstat
      // ierr = SNESLineSearchSetOrder(subls,2);   // change BT order to quadratic
      // ierr = SNESSetNPC(snes,subsnes);          // apply the changes

      PetscOptionsSetValue(NULL,"-npc_sub_ksp_type","fgmres");
      PetscOptionsSetValue(NULL,"-npc_sub_pc_type","ilu");
      PetscOptionsSetValue(NULL,"-npc_sub_snes_linesearch_order","2");
      // PetscOptionsSetValue(NULL,"-npc_sub_snes_rtol","1e-1");
      // PetscOptionsSetValue(NULL,"-npc_sub_ksp_rtol","1e-1");
      // PetscOptionsSetValue(NULL,"-npc_sub_snes_max_it","1");
   }

   if (mixed) { // each subdomain gets their own indexed prefix
      char prefix[10];
      sprintf(prefix,"sub_%d_",rank);
      SNESSetOptionsPrefix(subsnes,prefix);
   }
   ierr = SNESSetFromOptions(subsnes); CHKERRQ(ierr);
   ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
   /* Wide-stencil params - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Compute forward stencil directions for the determinant
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ComputeFwdStencilDirs(width,&user); // Stencil directions based on the L1 norm
   ComputeWeights(width,order,&user);  // Quadrature weights
   /* Solve - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Set up the inital guess on each sub-domain.
      The runtime is calculated for the SNESSolve only.
      The number of outer SNES iterations is also obtained.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   ierr = DMGetGlobalVector(da,&u_initial); CHKERRQ(ierr);
   ierr = InitialState(da,initial,u_initial,&user); CHKERRQ(ierr);
   ierr = PetscTime(&t1); CHKERRQ(ierr);
   ierr = SNESSolve(snes,NULL,u_initial); CHKERRQ(ierr);
   ierr = PetscTime(&t2); CHKERRQ(ierr);
   ierr = SNESGetIterationNumber(snes,&its); CHKERRQ(ierr); // its = # of iterations
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
      Problem type, grid size, stencil width, epsilon, iterations, and error
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   switch (dim) {
      case 1:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1));
         snprintf(gridstr,99,"domain [%.2f,%.2f] with %d point 1D",user.xmin,user.xmax,info.mx);
         break;
      case 2:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1)*(info.my-1));
         snprintf(gridstr,99,"domain [%.2f,%.2f]x[%.2f,%.2f] with %dx%d point 2D",user.xmin,user.xmax,user.ymin,user.ymax,info.mx,info.my);
         break;
      case 3:
         normconst2h = PetscSqrtReal((PetscReal)(info.mx-1)*(info.my-1)*(info.mz-1));
         snprintf(gridstr,99,"domain [%.2f,%.2f]x[%.2f,%.2f]x[%.2f,%.2f] %dx%dx%d point 3D",user.xmin,user.xmax,user.ymin,user.ymax,user.zmin,user.zmax,info.mx,info.my,info.mz);
         break;
      default:
         SETERRQ(PETSC_COMM_SELF,4,"invalid dim value in final report\n");
   }
   err2h /= normconst2h; // like continuous L2
   ierr = PetscPrintf(PETSC_COMM_WORLD,
               "*Problem: %s on %s grid\n"
               "*Params:  Nd = %d, width = %d, eps = %.3f, op = %.2f\n"
               "*Error:   |u-uexact|_inf = %.3e, |u-uexact|_h = %.3e\n"
               "*WTime:   %.6f\n"
               "*Iters:   %d\n",
               ProblemTypes[problem],gridstr,size,info.sw,user.epsilon,op,errinf,err2h,t2-t1,its); CHKERRQ(ierr);

   /* Debugging Info - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      To print out debugging info, add the option -t1_debug.
      The goal is to get info about the subdomains.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if (debug) {
      PetscInt MM, NN, mm, nn, dof, ss;
      DMDAGetInfo(da_after,&dim,&MM,&NN,NULL,&mm,&nn,NULL,&dof,&ss,NULL,NULL,NULL,NULL);
      PetscPrintf(PETSC_COMM_WORLD,"-- Grid is %d by %d, processors, and the domain is divided into %d stacks of %d.\n",MM,NN,mm,nn);
      PetscPrintf(PETSC_COMM_WORLD,"-- Dof = %d, Stencil = %d\n",dof,ss);
      PetscInt ox, oy;
      DMDAGetOverlap(da_after,&ox,&oy,NULL);
      PetscPrintf(PETSC_COMM_WORLD,"-- Overlap in x: %d, Overlap in y: %d\n",ox,oy);
      DMDALocalInfo  info;
      MPI_Barrier(PETSC_COMM_WORLD);
      for (int i=0; i<size; i++) {
         if (i==rank) {
            DMDAGetLocalInfo(da,&info);
            printf("-- Rank %2d: x indeces [%2d (%2d) to %2d (%2d)], y indeces [%2d (%2d) to %2d (%2d)] (ghost)\n",rank,info.xs,info.gxs,info.xs+info.xm-1,info.gxs+info.gxm-1,info.ys,info.gys,info.ys+info.ym-1,info.gys+info.gym-1);
         }
      }
   }
   /* MATLAB  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Create MATLAB files load_u.m and load_exact.m which loads
      the numerical and exact solutions into a workspace.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
   if (printSol) {
      PetscViewer viewer; // Viewer object fascilitates printing out solution
      ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_u.m",&viewer);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
      ierr = PetscObjectSetName((PetscObject)u,"u");
      ierr = VecView(u,viewer);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_exact.m",&viewer);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
      ierr = PetscObjectSetName((PetscObject)u_exact,"u_exact");
      ierr = VecView(u_exact,viewer);
      if (debug) { // extras for sanity check
         // print out intial guess in matlab format
         InitialState(da,initial,u_initial,&user);
         ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_initial.m",&viewer);
         ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
         ierr = PetscObjectSetName((PetscObject)u_initial,"u_initial");
         ierr = VecView(u_initial,viewer);
         Vec  RHS; // print out the source term or the RHS
         ierr = VecDuplicate(u_exact,&RHS);
         ierr = ComputeRHS(&info,RHS,&user);
         ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_rhs.m",&viewer);
         ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
         ierr = PetscObjectSetName((PetscObject)RHS,"rhs");
         ierr = VecView(RHS,viewer);
         ierr = VecDestroy(&RHS);
      }
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
   }
   // Free memory
   ierr = DMDestroy(&da); CHKERRQ(ierr); // destroying da also destroys snes
   ierr = VecDestroy(&err); CHKERRQ(ierr);
   ierr = VecDestroy(&u); CHKERRQ(ierr);
   ierr = VecDestroy(&u_exact); CHKERRQ(ierr);
   ierr = VecDestroy(&u_initial); CHKERRQ(ierr);
   return PetscFinalize();
}

PetscErrorCode Form1DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i;
   PetscReal  hx,x,*au,temp;
   hx = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   ierr = DMDAVecGetArray(info->da,u,&au); CHKERRQ(ierr);
   for (i=info->xs; i<info->xs+info->xm; i++) {
      x = user->xmin + (i+1)*hx;
      user->g_bdry(x,0.0,0.0,user,&temp);
      au[i] = temp;
   }
   ierr = DMDAVecRestoreArray(info->da,u,&au); CHKERRQ(ierr);
   return 0;
}

PetscErrorCode Form2DUExact(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i,j;
   PetscReal  hx,hy,x,y,**au,temp;
   hx = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   hy = (user->ymax-user->ymin)/(PetscReal)(info->my+1);
   ierr = DMDAVecGetArray(info->da,u,&au); CHKERRQ(ierr);
   for (j=info->ys; j<info->ys+info->ym; j++) {
      y = user->ymin + (j+1)*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = user->xmin + (i+1)*hx;
         user->g_bdry(x,y,0.0,user,&temp);
         au[j][i] = temp;
      }
   }
   ierr = DMDAVecRestoreArray(info->da,u,&au); CHKERRQ(ierr);
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
   ierr = DMDAVecGetArray(info->da, u, &au); CHKERRQ(ierr);
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
   ierr = DMDAVecRestoreArray(info->da, u, &au); CHKERRQ(ierr);
   return 0;
}

// Compute the RHS and store it in u. This exists as a sanity check.
PetscErrorCode ComputeRHS(DMDALocalInfo *info, Vec u, MACtx* user) {
   PetscErrorCode ierr;
   PetscInt   i, j;
   PetscReal  Lx, Ly, hx, hy, x, y, **au,temp;
   Lx = user->xmax - user->xmin;
   Ly = user->ymax - user->ymin;
   hx = Lx/(PetscReal)(info->mx + 1);
   hy = Ly/(PetscReal)(info->my + 1);
   ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
   for (j=info->ys; j<info->ys+info->ym; j++) {
      y = user->ymin + (j+1)*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = user->xmin + (i+1)*hx;
         user->f_rhs(x,y,0.0,user,&temp);
         au[j][i] = temp;
      }
   }
   ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
   return 0;
}