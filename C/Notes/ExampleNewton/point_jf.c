//STARTWHOLE
static char help[] = "Newton's method for a two-variable system.\n"
"This code is Jacobian free.  Run with -snes_fd or -snes_mf. This is a serial code.\n\n";
/*
   Find the point (x,y) where the two curves: y = 1/b*exp(b*x) and
   x^2 +y^2 + 1 intersect using Newton's method. Let b = 2.

   All we need to do is set up the problem as F(x,y) = [0;0].
   F is defined in `FormFunction`. We don't provide a Jacobian.

   The goal is to learn the syntax of the SNES.

   What's going on under the hood:
    - start with initial guess (1,1) which is quite close
    - forward diff to approx Jacobian, step is chosen internally
    - line search is used to scale the newton steps
    - iterates until relative or absolute residue reaches 10^-8
*/


#include <petsc.h>

extern PetscErrorCode FormFunction(SNES, Vec, Vec, void*);

int main(int argc,char **argv) {
   PetscErrorCode ierr;
   SNES  snes;          // scalable nonlinear equation solver
   Vec   x, r;          // solution, residual vectors

   ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;

   ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
   ierr = VecSetSizes(x,PETSC_DECIDE,2); CHKERRQ(ierr);
   ierr = VecSetFromOptions(x); CHKERRQ(ierr);
   ierr = VecSet(x,1.0); CHKERRQ(ierr);         // initial iterate
   ierr = VecDuplicate(x,&r); CHKERRQ(ierr);

   // initialize the SNES object
   ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);

   // Associate SNES object with F and residue vector
   ierr = SNESSetFunction(snes,r,FormFunction,NULL); CHKERRQ(ierr);

   // Apply command arg options
   /*
      run this code with `-snes_monitor` to see residue
      run this code with `-snes_rtol 1e-14` to lower deflt tol
   */
   ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

   // Call the non-linear solver
   ierr = SNESSolve(snes,NULL,x); CHKERRQ(ierr);

   // Print solution
   ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

   // Clean up
   SNESDestroy(&snes);  VecDestroy(&x);  VecDestroy(&r);
   return PetscFinalize();
}

// This is the function that we want to find the root of.
// The template for this function is very specific.
// You can see that `ctx` is never used but it needs to be there. It stands for "context".
PetscErrorCode FormFunction(SNES snes, Vec x, Vec F, void *ctx) {
   PetscErrorCode ierr;
   const PetscReal  b = 2.0, *ax;
   PetscReal        *aF;
   /*
      F is a 2 by 1 array which represents the residue so
         F[0] = 1/b * exp(b * x) - y
         F[1] = x^2 + y^2 - 1
      The vector x is the input so really we have
         x[0] = x
         x[1] = y
   */
   ierr = VecGetArrayRead(x,&ax);CHKERRQ(ierr);
   ierr = VecGetArray(F,&aF);CHKERRQ(ierr);
   aF[0] = (1.0 / b) * PetscExpReal(b * ax[0]) - ax[1];
   aF[1] = ax[0] * ax[0] + ax[1] * ax[1] - 1.0;
   ierr = VecRestoreArrayRead(x,&ax);CHKERRQ(ierr);
   ierr = VecRestoreArray(F,&aF);CHKERRQ(ierr);
   return 0;
}
//ENDWHOLE