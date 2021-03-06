static char help[] = "1D reaction-diffusion problem with DMDA and SNES.  Option prefix -rct_.\n\n";

/*
   The 1D grid is broken up into sub-intervals automatically with DMDA.
   The main motivation to use DMDA with SNES is to use something called "coloring"
   which works well with Finite Diff (a type of Jacobian-free) newton's method.
   This code also allows the user to use the Jacobian free Newton-Krylov method.
   You can even choose which Krylov method you want to use.
   It's recommended to use conjugate gradient since the Jacobian is symmetric (?).


   Run with `-snes_fd_color` to use Finite Diff with the node coloring.
   Run with `-snes_mf` to use the matrix free Newton Krylov method.
   The Krylov is unconditioned and apparently not that good.
   Run with `-snes_mf -ksp_type cg -ksp

   This code comes with a special feature, an inexact Jacobian.
   Using this inexact Jacobian directly leads to more Newton iterations.
   To compare, run   `./diffeq1d -rct_noRinJ -snes_monitor` vs `./diffeq1d -snes_monitor`.
   However, this inexact Jacobian is actually a good preconditioner for the Jacobian-free Newton Krylov method.
   Add the option `-snes_mf_operator` to use the preconditioned NK method.
   To see why its good run `./diffeq1d -snes_monitor -rct_noRinJ -snes_mf_operator -ksp_converged_reason -da_refine 5`.
   The `ksp_converged_reason` tells us how the Krylov performed at each Newton iteration.

*/
#include <petsc.h>

typedef struct {
    PetscReal  rho, M, alpha, beta;
    PetscBool  noRinJ;
} AppCtx;

extern PetscReal f_source(PetscReal);
extern PetscErrorCode InitialAndExact(DMDALocalInfo*, PetscReal*, PetscReal*, AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, PetscReal*, PetscReal*, AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*, PetscReal*, Mat, Mat, AppCtx*);

//STARTMAIN
int main(int argc,char **args) {
   PetscErrorCode ierr;
   DM            da;
   SNES          snes;
   AppCtx        user;
   Vec           u, uexact;           // u and uexact are petsc vectors
   PetscReal     errnorm, *au, *auex; // au and auex are C arrays
   DMDALocalInfo info;

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;
   user.rho   = 10.0;
   user.M     = PetscSqr(user.rho / 12.0);
   user.alpha = user.M;
   user.beta  = 16.0 * user.M;
   user.noRinJ = PETSC_FALSE;


   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"rct_","options for reaction",""); CHKERRQ(ierr);
   ierr = PetscOptionsBool("-noRinJ","do not include R(u) term in Jacobian",
         "diffeq1d.c",user.noRinJ,&(user.noRinJ),NULL); CHKERRQ(ierr);
   ierr = PetscOptionsEnd(); CHKERRQ(ierr);

   // 1d DMDA
   // ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,9,1,1,NULL,&da); CHKERRQ(ierr);
   ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,9,1,1,NULL,&da); CHKERRQ(ierr);

   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);

   // Any vector representing a function defined over the domain needs
   // to be created in association with the `da` object.
   ierr = DMCreateGlobalVector(da,&u); CHKERRQ(ierr);
   ierr = VecDuplicate(u,&uexact); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da,u,&au); CHKERRQ(ierr);

   // unsure what these are
   ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
   ierr = DMDAVecGetArray(da,uexact,&auex); CHKERRQ(ierr);

   // this is a custom function
   ierr = InitialAndExact(&info,au,auex,&user); CHKERRQ(ierr);

   ierr = DMDAVecRestoreArray(da,u,&au); CHKERRQ(ierr);
   ierr = DMDAVecRestoreArray(da,uexact,&auex); CHKERRQ(ierr);

   // SNES is the nonlinear solver object
   // The snes object and the da object need to be associated
   ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
   ierr = SNESSetDM(snes,da); CHKERRQ(ierr);

   // This is the syntax to associate user-defined Jacobians and Residue functions to the snes
   // Instead of being linked directly to snes, its linked to the da.
   ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
   ierr = DMDASNESSetJacobianLocal(da,(DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
   ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

   // call the solver
   ierr = SNESSolve(snes,NULL,u); CHKERRQ(ierr);

   ierr = VecAXPY(u,-1.0,uexact); CHKERRQ(ierr);    // u <- u + (-1.0) uexact
   ierr = VecNorm(u,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"on %d point grid:  |u-u_exact|_inf = %g\n",info.mx,errnorm); CHKERRQ(ierr);

   VecDestroy(&u);  VecDestroy(&uexact);
   SNESDestroy(&snes);  DMDestroy(&da);
   return PetscFinalize();
}
//ENDMAIN

PetscReal f_source(PetscReal x) {
   return 0.0;
}


/*

*/
PetscErrorCode InitialAndExact(DMDALocalInfo *info, PetscReal *u0, PetscReal *uex, AppCtx *user) {
   PetscInt   i;
   PetscReal  h = 1.0 / (info->mx-1), x;
   for (i=info->xs; i<info->xs+info->xm; i++) {
      x = h * i;
      u0[i]  = user->alpha * (1.0 - x) + user->beta * x;
      uex[i] = user->M * PetscPowReal(x + 1.0,4.0);
   }
   return 0;
}

//STARTFUNCTIONS
PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, PetscReal *u, PetscReal *FF, AppCtx *user) {
   PetscInt   i;
   PetscReal  h = 1.0 / (info->mx-1), x, R;
   for (i=info->xs; i<info->xs+info->xm; i++) {
      if (i == 0) {
         FF[i] = u[i] - user->alpha;
      } else if (i == info->mx-1) {
         FF[i] = u[i] - user->beta;
      } else {  // interior location
         if (i == 1) {
               FF[i] = - u[i+1] + 2.0 * u[i] - user->alpha;
         } else if (i == info->mx-2) {
               FF[i] = - user->beta + 2.0 * u[i] - u[i-1];
         } else {
               FF[i] = - u[i+1] + 2.0 * u[i] - u[i-1];
         }
         R = - user->rho * PetscSqrtReal(u[i]);
         x = i * h;
         FF[i] -= h*h * (R + f_source(x));
      }
   }
   return 0;
}

PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, PetscReal *u, Mat J, Mat P, AppCtx *user) {
   PetscErrorCode ierr;
   PetscInt   i, col[3];
   PetscReal  h = 1.0 / (info->mx-1), dRdu, v[3];
   for (i=info->xs; i<info->xs+info->xm; i++) {
      if ((i == 0) | (i == info->mx-1)) {
         v[0] = 1.0;
         ierr = MatSetValues(P,1,&i,1,&i,v,INSERT_VALUES); CHKERRQ(ierr);
      } else {
         col[0] = i;
         v[0] = 2.0;
         if (!user->noRinJ) {
               dRdu = - (user->rho / 2.0) / PetscSqrtReal(u[i]);
               v[0] -= h*h * dRdu;
         }
         col[1] = i-1;   v[1] = (i > 1) ? - 1.0 : 0.0;
         col[2] = i+1;   v[2] = (i < info->mx-2) ? - 1.0 : 0.0;
         ierr = MatSetValues(P,1,&i,3,col,v,INSERT_VALUES); CHKERRQ(ierr);
      }
   }
   ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   if (J != P) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   }
   return 0;
}
//ENDFUNCTIONS