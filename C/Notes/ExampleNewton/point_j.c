static char help[] = "Newton's method for a two-variable system.  Implements analytical Jacobian and a struct to hold a parameter.\n";

#include <petsc.h>

/*
   Again, we find the point (x,y) where the two curves: y = 1/b*exp(b*x) and
   x^2 +y^2 + 1 intersect using Newton's method. This time, the Jacobian of
   the system is explicitly coded in.

   Running `./point_j -snes_test_jacobian` tells petsc to see if the Finite Diff and
   the hand-coded Jacobian match up.

*/




// This is custom struct. AppCtx stands for "application context".
// It's used to pass bunch of parameters to the Jacobian in one neat struct.
// Actually, in this example, there's only 1 parameter.
typedef struct {
  PetscReal  b;
} AppCtx;

extern PetscErrorCode FormFunction(SNES, Vec, Vec, void*);
extern PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void*);

//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  SNES   snes;         // nonlinear solver
  Vec    x,r;          // solution, residual vectors
  Mat    J;
  AppCtx user;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.b = 2.0;

  ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,2); CHKERRQ(ierr);
  ierr = VecSetFromOptions(x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r); CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&J); CHKERRQ(ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,2,2); CHKERRQ(ierr);
  ierr = MatSetFromOptions(J); CHKERRQ(ierr);
  ierr = MatSetUp(J); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,r,FormFunction,&user);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,FormJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  ierr = VecSet(x,1.0);CHKERRQ(ierr);            // initial iterate
  ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  SNESDestroy(&snes);  MatDestroy(&J);  VecDestroy(&x);  VecDestroy(&r);
  return PetscFinalize();
}
//ENDMAIN

//STARTFORM
PetscErrorCode FormFunction(SNES snes, Vec x, Vec F, void *ctx) {
   PetscErrorCode ierr;
   AppCtx          *user = (AppCtx*)ctx;
   const PetscReal b = user->b, *ax;
   PetscReal       *aF;

   ierr = VecGetArrayRead(x,&ax);CHKERRQ(ierr);
   ierr = VecGetArray(F,&aF);CHKERRQ(ierr);
   aF[0] = (1.0 / b) * PetscExpReal(b * ax[0]) - ax[1];
   aF[1] = ax[0] * ax[0] + ax[1] * ax[1] - 1.0;
   ierr = VecRestoreArrayRead(x,&ax);CHKERRQ(ierr);
   ierr = VecRestoreArray(F,&aF);CHKERRQ(ierr);
   return 0;
}

/*
This function contains the rule to construct the Jacobian. It has to obey a specific structure.
   The 1st arg is always the SNES object.
   The 2nd arg is the current iterate x_k
   The 3rd arg is the Jacobian itself
   The 4th arg is the material used to precondition the Jacobian system.
   The 5th arg is anything that parameters for the jacobian.
The user is responsible for constructing P.
*/
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat P, void *ctx) {
   PetscErrorCode ierr;
   AppCtx           *user = (AppCtx*)ctx;
   const PetscReal  b = user->b, *ax;
   PetscReal        v[4];
   PetscInt         row[2] = {0,1}, col[2] = {0,1};

   ierr = VecGetArrayRead(x,&ax); CHKERRQ(ierr);
   v[0] = PetscExpReal(b * ax[0]);  v[1] = -1.0;
   v[2] = 2.0 * ax[0];              v[3] = 2.0 * ax[1];
   ierr = VecRestoreArrayRead(x,&ax); CHKERRQ(ierr);
   ierr = MatSetValues(P,2,row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);
   ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   if (J != P) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   }
   return 0;
}
//ENDFORM

