static char help[] =
"Solves structured-grid Poisson problem in 1D, 2D, 3D.  Option prefix fsh_.\n"
"Equation is\n"
"    - cx u_xx - cy u_yy - cz u_zz = f,\n"
"subject to Dirichlet boundary conditions.  Solves three different problems\n"
"where exact solution is known.  Uses DMDA and SNES.  Equation is put in form\n"
"F(u) = - grad^2 u - f.  Call-backs fully-rediscretize for the supplied grid.\n"
"Defaults to 2D, a SNESType of KSPONLY, and a KSPType of CG.\n\n";

#include <petsc.h>
#include "poissonfunctions.h"

// exact solutions  u(x,y),  for boundary condition and error calculation

static PetscReal u_exact_1Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x);
}

static PetscReal u_exact_2Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x) * y*y *(y*y - 1.0);
}

static PetscReal u_exact_3Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x*x * (1.0 - x*x) * y*y * (y*y - 1.0) * z*z * (z*z - 1.0);
}

static PetscReal u_exact_1Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - PetscExpReal(x);
}

static PetscReal u_exact_2Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - x * PetscExpReal(y);
}

static PetscReal u_exact_3Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return - x * PetscExpReal(y + z);
}

static PetscReal zero(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return 0.0;
}

// right-hand-side functions  f(x,y) = - laplacian u

static PetscReal f_rhs_1Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    PoissonCtx* user = (PoissonCtx*)ctx;
    return user->cx * 12.0 * x*x - 2.0;
}

static PetscReal f_rhs_2Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    PoissonCtx* user = (PoissonCtx*)ctx;
    PetscReal   aa, bb, ddaa, ddbb;
    aa = x*x * (1.0 - x*x);
    bb = y*y * (y*y - 1.0);
    ddaa = 2.0 * (1.0 - 6.0 * x*x);
    ddbb = 2.0 * (6.0 * y*y - 1.0);
    return - (user->cx * ddaa * bb + user->cy * aa * ddbb);
}

static PetscReal f_rhs_3Dmanupoly(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    PoissonCtx* user = (PoissonCtx*)ctx;
    PetscReal   aa, bb, cc, ddaa, ddbb, ddcc;
    aa = x*x * (1.0 - x*x);
    bb = y*y * (y*y - 1.0);
    cc = z*z * (z*z - 1.0);
    ddaa = 2.0 * (1.0 - 6.0 * x*x);
    ddbb = 2.0 * (6.0 * y*y - 1.0);
    ddcc = 2.0 * (6.0 * z*z - 1.0);
    return - (user->cx * ddaa * bb * cc + user->cy * aa * ddbb * cc + user->cz * aa * bb * ddcc);
}

static PetscReal f_rhs_1Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return PetscExpReal(x);
}

static PetscReal f_rhs_2Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return x * PetscExpReal(y);  // note  f = - (u_xx + u_yy) = - u
}

static PetscReal f_rhs_3Dmanuexp(PetscReal x, PetscReal y, PetscReal z, void *ctx) {
    return 2.0 * x * PetscExpReal(y + z);  // note  f = - laplacian u = - 2 u
}

// functions simply to put u_exact()=g_bdry() into a grid
// these are irritatingly-dimension-dependent inside ...
extern PetscErrorCode Form1DUExact(DMDALocalInfo*, Vec, PoissonCtx*);
extern PetscErrorCode Form2DUExact(DMDALocalInfo*, Vec, PoissonCtx*);
extern PetscErrorCode Form3DUExact(DMDALocalInfo*, Vec, PoissonCtx*);

//STARTPTRARRAYS
// arrays of pointers to functions
static DMDASNESFunction residual_ptr[3]
    = {(DMDASNESFunction)&Poisson1DFunctionLocal,
       (DMDASNESFunction)&Poisson2DFunctionLocal,
       (DMDASNESFunction)&Poisson3DFunctionLocal};

static DMDASNESJacobian jacobian_ptr[3]
    = {(DMDASNESJacobian)&Poisson1DJacobianLocal,
       (DMDASNESJacobian)&Poisson2DJacobianLocal,
       (DMDASNESJacobian)&Poisson3DJacobianLocal};

typedef PetscErrorCode (*ExactFcnVec)(DMDALocalInfo*,Vec,PoissonCtx*);

static ExactFcnVec getuexact_ptr[3]
    = {&Form1DUExact, &Form2DUExact, &Form3DUExact};
//ENDPTRARRAYS

typedef enum {MANUPOLY, MANUEXP, ZERO} ProblemType;
static const char* ProblemTypes[] = {"manupoly","manuexp","zero",
                                     "ProblemType", "", NULL};

// more arrays of pointers to functions:   ..._ptr[DIMS][PROBLEMS]
typedef PetscReal (*PointwiseFcn)(PetscReal,PetscReal,PetscReal,void*);

static PointwiseFcn g_bdry_ptr[3][3]
    = {{&u_exact_1Dmanupoly, &u_exact_1Dmanuexp, &zero},
       {&u_exact_2Dmanupoly, &u_exact_2Dmanuexp, &zero},
       {&u_exact_3Dmanupoly, &u_exact_3Dmanuexp, &zero}};

static PointwiseFcn f_rhs_ptr[3][3]
    = {{&f_rhs_1Dmanupoly, &f_rhs_1Dmanuexp, &zero},
       {&f_rhs_2Dmanupoly, &f_rhs_2Dmanuexp, &zero},
       {&f_rhs_3Dmanupoly, &f_rhs_3Dmanuexp, &zero}};

static const char* InitialTypes[] = {"zeros","random",
                                     "InitialType", "", NULL};

int main(int argc,char **argv) {
    PetscErrorCode ierr;
    DM             da, da_after;
    SNES           snes;
    KSP            ksp;
    Vec            u_initial, u, u_exact;
    PoissonCtx     user;
    DMDALocalInfo  info;
    PetscReal      errinf, normconst2h, err2h;
    char           gridstr[99];
    ExactFcnVec    getuexact;

    // fish defaults:
    PetscInt       dim = 2;                  // 2D
    ProblemType    problem = MANUEXP;        // manufactured problem using exp()
    InitialType    initial = ZEROS;          // set u=0 for initial iterate
    PetscBool      gonboundary = PETSC_TRUE; // initial iterate has u=g on boundary

    ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;

    // get options and configure context
    user.Lx = 1.0;
    user.Ly = 1.0;
    user.Lz = 1.0;
    user.cx = 1.0;
    user.cy = 1.0;
    user.cz = 1.0;
    ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"fsh_", "options for fish.c", ""); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cx",
         "set coefficient of x term u_xx in equation",
         "fish.c",user.cx,&user.cx,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cy",
         "set coefficient of y term u_yy in equation",
         "fish.c",user.cy,&user.cy,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cz",
         "set coefficient of z term u_zz in equation",
         "fish.c",user.cz,&user.cz,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsInt("-dim",
         "dimension of problem (=1,2,3 only)",
         "fish.c",dim,&dim,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-initial_gonboundary",
         "set initial iterate to have correct boundary values",
         "fish.c",gonboundary,&gonboundary,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsEnum("-initial_type",
         "type of initial iterate",
         "fish.c",InitialTypes,(PetscEnum)initial,(PetscEnum*)&initial,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Lx",
         "set Lx in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
         "fish.c",user.Lx,&user.Lx,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Ly",
         "set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
         "fish.c",user.Ly,&user.Ly,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Lz",
         "set Ly in domain ([0,Lx] x [0,Ly] x [0,Lz], etc.)",
         "fish.c",user.Lz,&user.Lz,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsEnum("-problem",
         "problem type; determines exact solution and RHS",
         "fish.c",ProblemTypes,(PetscEnum)problem,(PetscEnum*)&problem,NULL); CHKERRQ(ierr);
    ierr = PetscOptionsEnd(); CHKERRQ(ierr);
    user.g_bdry = g_bdry_ptr[dim-1][problem];
    user.f_rhs = f_rhs_ptr[dim-1][problem];
    if ( user.cx <= 0.0 || user.cy <= 0.0 || user.cz <= 0.0 ) {
        SETERRQ(PETSC_COMM_SELF,2,"positivity required for coefficients cx,cy,cz\n");
    }
    if ((problem == MANUEXP) && ( user.cx != 1.0 || user.cy != 1.0 || user.cz != 1.0)) {
        SETERRQ(PETSC_COMM_SELF,3,"cx=cy=cz=1 required for problem MANUEXP\n");
    }

//STARTCREATE
    // create DMDA in chosen dimension
    switch (dim) {
        case 1:
            ierr = DMDACreate1d(PETSC_COMM_WORLD,
                DM_BOUNDARY_NONE,3,1,1, NULL, &da); CHKERRQ(ierr);
            break;
        case 2:
            ierr = DMDACreate2d(PETSC_COMM_WORLD,
                DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,
                3,3,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da); CHKERRQ(ierr);
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
    ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);  // call BEFORE SetUniformCoordinates
    ierr = DMDASetUniformCoordinates(da,0.0,user.Lx,0.0,user.Ly,0.0,user.Lz); CHKERRQ(ierr);

    // set SNES call-backs
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
    ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
    ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)(residual_ptr[dim-1]),&user); CHKERRQ(ierr);
    ierr = DMDASNESSetJacobianLocal(da,
             (DMDASNESJacobian)(jacobian_ptr[dim-1]),&user); CHKERRQ(ierr);

    // default to KSPONLY+CG because problem is linear and SPD
    ierr = SNESSetType(snes,SNESKSPONLY); CHKERRQ(ierr);
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    // set initial iterate and then solve
    ierr = DMGetGlobalVector(da,&u_initial); CHKERRQ(ierr);
    ierr = InitialState(da, initial, gonboundary, u_initial, &user); CHKERRQ(ierr);
    ierr = SNESSolve(snes,NULL,u_initial); CHKERRQ(ierr);
//ENDCREATE

//STARTGETSOLUTION
    // -snes_grid_sequence could change grid resolution
    ierr = DMRestoreGlobalVector(da,&u_initial); CHKERRQ(ierr);
    ierr = DMDestroy(&da); CHKERRQ(ierr);

    // evaluate error and report
    ierr = SNESGetSolution(snes,&u); CHKERRQ(ierr);  // SNES owns u; do not destroy it
    ierr = SNESGetDM(snes,&da_after); CHKERRQ(ierr); // SNES owns da_after; do not destroy it
    ierr = DMDAGetLocalInfo(da_after,&info); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da_after,&u_exact); CHKERRQ(ierr);
    getuexact = getuexact_ptr[dim-1];
    ierr = (*getuexact)(&info,u_exact,&user); CHKERRQ(ierr);
    ierr = VecAXPY(u,-1.0,u_exact); CHKERRQ(ierr);   // u <- u + (-1.0) uexact
    ierr = VecDestroy(&u_exact); CHKERRQ(ierr);      // no longer needed
    ierr = VecNorm(u,NORM_INFINITY,&errinf); CHKERRQ(ierr);
    ierr = VecNorm(u,NORM_2,&err2h); CHKERRQ(ierr);
//ENDGETSOLUTION

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
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                "problem %s on %s grid:\n"
                "  error |u-uexact|_inf = %.3e, |u-uexact|_h = %.3e\n",
                ProblemTypes[problem],gridstr,errinf,err2h); CHKERRQ(ierr);

    // destroy what we explicitly Created
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    return PetscFinalize();
}

PetscErrorCode Form1DUExact(DMDALocalInfo *info, Vec u, PoissonCtx* user) {
  PetscErrorCode ierr;
  PetscInt   i;
  PetscReal  xmax[1], xmin[1], hx, x, *au;
  ierr = DMGetBoundingBox(info->da,xmin,xmax); CHKERRQ(ierr);
  hx = (xmax[0] - xmin[0]) / (info->mx - 1);
  ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
  for (i=info->xs; i<info->xs+info->xm; i++) {
      x = xmin[0] + i * hx;
      au[i] = user->g_bdry(x,0.0,0.0,user);
  }
  ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
  return 0;
}

PetscErrorCode Form2DUExact(DMDALocalInfo *info, Vec u, PoissonCtx* user) {
    PetscErrorCode ierr;
    PetscInt   i, j;
    PetscReal  xymin[2], xymax[2], hx, hy, x, y, **au;
    ierr = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
    hx = (xymax[0] - xymin[0]) / (info->mx - 1);
    hy = (xymax[1] - xymin[1]) / (info->my - 1);
    ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
    for (j=info->ys; j<info->ys+info->ym; j++) {
        y = xymin[1] + j * hy;
        for (i=info->xs; i<info->xs+info->xm; i++) {
            x = xymin[0] + i * hx;
            au[j][i] = user->g_bdry(x,y,0.0,user);
        }
    }
    ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
    return 0;
}

PetscErrorCode Form3DUExact(DMDALocalInfo *info, Vec u, PoissonCtx* user) {
    PetscErrorCode ierr;
    PetscInt  i, j, k;
    PetscReal xyzmin[3], xyzmax[3], hx, hy, hz, x, y, z, ***au;
    ierr = DMGetBoundingBox(info->da,xyzmin,xyzmax); CHKERRQ(ierr);
    hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
    hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
    hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
    ierr = DMDAVecGetArray(info->da, u, &au);CHKERRQ(ierr);
    for (k=info->zs; k<info->zs+info->zm; k++) {
        z = xyzmin[2] + k * hz;
        for (j=info->ys; j<info->ys+info->ym; j++) {
            y = xyzmin[1] + j * hy;
            for (i=info->xs; i<info->xs+info->xm; i++) {
                x = xyzmin[0] + i * hx;
                au[k][j][i] = user->g_bdry(x,y,z,user);
            }
        }
    }
    ierr = DMDAVecRestoreArray(info->da, u, &au);CHKERRQ(ierr);
    return 0;
}

