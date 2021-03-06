//STARTWHOLE
static char help[] = "Solve a tridiagonal system of arbitrary size.\n"
"Option prefix = tri_.\n";

#include <petsc.h>

int main(int argc,char **args) {
   PetscErrorCode ierr;
   Vec         x, b, xexact;
   Mat         A;
   KSP         ksp;
   PetscInt    m = 4, i, Istart, Iend, j[3];
   PetscReal   v[3], xval, errnorm;

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;

   // Set up args
   /*
      m is the size of the system
      m is 4 by default
      m can be changed during runtime
   */
   ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"tri_","options for tri",""); CHKERRQ(ierr);
   ierr = PetscOptionsInt("-m","dimension of linear system","tri.c",m,&m,NULL); CHKERRQ(ierr);
   ierr = PetscOptionsEnd(); CHKERRQ(ierr);

   // Set up vectors
   /*
      x is a vector of size m
      x is distributed automatically across some number P processors
      P is chosen during runtime, and P=1 by default
      b is a clone of x
      xexact is another clone
   */
   ierr = VecCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
   ierr = VecSetSizes(x,PETSC_DECIDE,m); CHKERRQ(ierr);
   ierr = VecSetFromOptions(x); CHKERRQ(ierr);
   ierr = VecDuplicate(x,&b); CHKERRQ(ierr);
   ierr = VecDuplicate(x,&xexact); CHKERRQ(ierr);

   // Set up A matrix
   /*
      There's a lot going on here...
      A is a m by m tridiagonal system with entries -1,3,-1
      The entries by rows. Only the few non-zero entries are specified.
      The matrix is assembled at the end in parallel.
   */
   ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
   ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m,m); CHKERRQ(ierr);
   ierr = MatSetOptionsPrefix(A,"a_"); CHKERRQ(ierr);
   ierr = MatSetFromOptions(A); CHKERRQ(ierr);
   ierr = MatSetUp(A); CHKERRQ(ierr);
   ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);
   for (i=Istart; i<Iend; i++) {
      if (i == 0) {
         v[0] = 3.0;  v[1] = -1.0;
         j[0] = 0;    j[1] = 1;
         ierr = MatSetValues(A,1,&i,2,j,v,INSERT_VALUES); CHKERRQ(ierr);
      } else {
         v[0] = -1.0;  v[1] = 3.0;  v[2] = -1.0;
         j[0] = i-1;   j[1] = i;    j[2] = i+1;
         if (i == m-1) {
               ierr = MatSetValues(A,1,&i,2,j,v,INSERT_VALUES); CHKERRQ(ierr);
         } else {
               ierr = MatSetValues(A,1,&i,3,j,v,INSERT_VALUES); CHKERRQ(ierr);
         }
      }
   }
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

   // Set up exact solution
   /*
      the exact solution is [1,2,3,...,m]^T
      the exact solution is assembled in parallel
      the RHS is computed with Mat-Vec multiplication.
   */
   for (i=Istart; i<Iend; i++) {
      xval = PetscExpReal(PetscCosReal((double)i));
      ierr = VecSetValues(xexact,1,&i,&xval,INSERT_VALUES); CHKERRQ(ierr);
   }
   ierr = VecAssemblyBegin(xexact); CHKERRQ(ierr);
   ierr = VecAssemblyEnd(xexact); CHKERRQ(ierr);
   ierr = MatMult(A,xexact,b); CHKERRQ(ierr);

   // set up the solver
   ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
   ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);
   ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

   // solve linear system Ax=b
   ierr = KSPSolve(ksp,b,x); CHKERRQ(ierr);

   // compute error
   ierr = VecAXPY(x,-1.0,xexact); CHKERRQ(ierr);
   ierr = VecNorm(x,NORM_2,&errnorm); CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,
   "error for m = %d system is |x-xexact|_2 = %.1e\n",m,errnorm); CHKERRQ(ierr);

   KSPDestroy(&ksp);  MatDestroy(&A);
   VecDestroy(&x);  VecDestroy(&b);  VecDestroy(&xexact);
   return PetscFinalize();
}
//ENDWHOLE
