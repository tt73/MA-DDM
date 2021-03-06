/*
   This is a sample code from Bueler's Petsc book.
   It solves the 2D Poisson problem on a square with a Krylov method.

   Goals:
   * Understand how the DMDACreate2d function works and all of its features
   * Learn how to interface PETSc with Matlab

   Compile this code by typing `make clr` and then `make poisson`.
   Run this code by typing: `./poisson`.
   Run this code with Np threads with an M by N grid by typing
      `mpiexec -n Np ./poisson -da_grid_x M -da_grid_y N`
*/
static char help[] = "A structured-grid Poisson solver using DMDA+KSP.\n\n";

#include <petsc.h>

extern PetscErrorCode formMatrix(DM, Mat);
extern PetscErrorCode formExact(DM, Vec);
extern PetscErrorCode formRHS(DM, Vec);

//STARTMAIN
int main(int argc,char **args) {
   Mat           A;           // Problem: Au=b
   Vec           b,u,uexact;  //
   DM            da;          // DataManager  - distributed array
   KSP           ksp;         // Krylov object fascilitates solving Au=b
   PetscViewer   viewer;      // Viewer object fascilitates printing out A, u, b
   PetscReal     errnorm;     // errnorm = |u-uexact|
   DMDALocalInfo info;        // Info contains parallel info which is different for each thread
   PetscErrorCode ierr;       // standard error handling variable

   ierr = PetscInitialize(&argc,&args,NULL,help); if (ierr) return ierr;

   // change default 9x9 size using -da_grid_x M -da_grid_y N
   ierr = DMDACreate2d(PETSC_COMM_WORLD,
                       DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                       DMDA_STENCIL_STAR,
                       9,9,
                       PETSC_DECIDE,PETSC_DECIDE,
                       1,1,NULL,NULL,
                       &da);
   CHKERRQ(ierr);

   // Initialize and apply the given settings to `da` and `A`
   ierr = DMSetFromOptions(da); CHKERRQ(ierr);
   ierr = DMSetUp(da); CHKERRQ(ierr);
   ierr = DMCreateMatrix(da,&A); CHKERRQ(ierr);
   ierr = MatSetFromOptions(A); CHKERRQ(ierr);

   // Initialize RHS b, approx solution u, exact solution uexact
   ierr = DMCreateGlobalVector(da,&b); CHKERRQ(ierr);
   ierr = VecDuplicate(b,&u); CHKERRQ(ierr);
   ierr = VecDuplicate(b,&uexact); CHKERRQ(ierr);

   // fill vectors and assemble linear system
   ierr = formExact(da,uexact); CHKERRQ(ierr);
   ierr = formRHS(da,b); CHKERRQ(ierr);
   ierr = formMatrix(da,A); CHKERRQ(ierr);

   // create and solve the linear system
   ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
   ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);
   ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
   ierr = KSPSolve(ksp,b,u); CHKERRQ(ierr);

   // report on grid and numerical error
   ierr = VecAXPY(u,-1.0,uexact); CHKERRQ(ierr);    // u <- u + (-1.0) uxact
   ierr = VecNorm(u,NORM_INFINITY,&errnorm); CHKERRQ(ierr);
   ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,
               "on %d x %d grid:  error |u-uexact|_inf = %g\n",
               info.mx,info.my,errnorm); CHKERRQ(ierr);


   // Print sparsity pattern
   // MatView(A,PETSC_VIEWER_DRAW_WORLD); // this doesn't work, need X windows

   // Set matrix and vector print format
   PetscViewerCreate(PETSC_COMM_WORLD, &viewer); // initialize the viewer object

   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_mat.m",&viewer);  // set viewer to print out to "load_mat.m"
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // set viewer to print in matlab syntax
   PetscObjectSetName((PetscObject)A,"A");                  // set the variable name of A to A in matlab
   MatView(A,viewer);                                       // print out the matrix

   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_b.m",&viewer);  // these commands are the same as the ones above
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for b
   PetscObjectSetName((PetscObject)b,"b");
   VecView(b,viewer);

   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"load_u.m",&viewer);  // these commands are the same as the ones above
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // except its for u
   PetscObjectSetName((PetscObject)u,"u");
   VecView(u,viewer);

   // Done. Clean up
   VecDestroy(&u);  VecDestroy(&uexact);  VecDestroy(&b);
   MatDestroy(&A);  KSPDestroy(&ksp);     DMDestroy(&da);
   PetscViewerDestroy(&viewer); // the viewer need to be destroyed as well
   return PetscFinalize();
}
//ENDMAIN

//STARTMATRIX
PetscErrorCode formMatrix(DM da, Mat A) {
   PetscErrorCode ierr;
   DMDALocalInfo  info;
   MatStencil     row, col[5];  // stencil variable - define i & j
   PetscReal      hx, hy, v[5];
   PetscInt       i, j, ncols;

   // Each thread gets its own unique values stored in `info`
   ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

   // info.mx = # of local horizontal nodes, info.my # of local vertical nodes (ends included)
   hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);

   // Each thread loops over its own range of (ys < j < ym),(xs < i < xm)
   for (j = info.ys; j < info.ys+info.ym; j++) {
      for (i = info.xs; i < info.xs+info.xm; i++) {

         /*
            The entries of A get filled with a subroutine called
            `MatSetValuesStencil`. You need to create `MatStencil`
            variables. This example is a 5-point star stencil.
         */
         row.j = j;           // row of A corresponding to (x_i,y_j)
         row.i = i;
         col[0].j = j;        // diagonal entry
         col[0].i = i;
         ncols = 1;   // there's at least 1 col to affect
         if (i==0 || i==info.mx-1 || j==0 || j==info.my-1) {
               v[0] = 1.0;      // on boundary: trivial equation
         } else {
               v[0] = 2*(hy/hx + hx/hy); // interior: build a row
               if (i-1 > 0) {
                  col[ncols].j = j;
                  col[ncols].i = i-1;
                  v[ncols++] = -hy/hx; // west
               }
               if (i+1 < info.mx-1) {
                  col[ncols].j = j;
                  col[ncols].i = i+1;
                  v[ncols++] = -hy/hx; // east
               }
               if (j-1 > 0) {
                  col[ncols].j = j-1;
                  col[ncols].i = i;
                  v[ncols++] = -hx/hy; // south
               }
               if (j+1 < info.my-1) {
                  col[ncols].j = j+1;
                  col[ncols].i = i;
                  v[ncols++] = -hx/hy; // north
               }
         }
         // insert values for 1 row, and ncol columns (which is usually 3, but can be 2 or 1)
         ierr = MatSetValuesStencil(A,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
      }
   }
   // this is the standard assembly procedure to create the matrix in parallel
   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   return 0;
}
//ENDMATRIX


//STARTEXACT
/*
   We are solving for a poisson problem with homogeneous BC.
   The exact solution is
      U(x,y) = (x^2 - x^4)*(y^2-y^4)
   on [0,1]x[0,1].
*/
PetscErrorCode formExact(DM da, Vec uexact) {
   PetscErrorCode ierr;
   PetscInt       i, j;
   PetscReal      hx, hy, x, y, **auexact;
   DMDALocalInfo  info;

   // Each thread gets unique info variable
   ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
   hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);

   /*
      The parallel vector assembly process sandwiched between
      GetArray and RestoreArray.
   */
   ierr = DMDAVecGetArray(da, uexact, &auexact);CHKERRQ(ierr);
   for (j = info.ys; j < info.ys+info.ym; j++) {
      y = j * hy;
      for (i = info.xs; i < info.xs+info.xm; i++) {
         x = i * hx;
         auexact[j][i] = x*x * (1.0 - x*x) * y*y * (y*y - 1.0);
      }
   }
   ierr = DMDAVecRestoreArray(da, uexact, &auexact);CHKERRQ(ierr);


   return 0;
}


PetscErrorCode formRHS(DM da, Vec b) {
   PetscErrorCode ierr;
   PetscInt       i, j;
   PetscReal      hx, hy, x, y, f, **ab;
   DMDALocalInfo  info;

   ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
   hx = 1.0/(info.mx-1);  hy = 1.0/(info.my-1);
   ierr = DMDAVecGetArray(da, b, &ab);CHKERRQ(ierr);
   for (j=info.ys; j<info.ys+info.ym; j++) {
      y = j * hy;
      for (i=info.xs; i<info.xs+info.xm; i++) {
         x = i * hx;
         if (i==0 || i==info.mx-1 || j==0 || j==info.my-1) {
               ab[j][i] = 0.0;  // on boundary: 1*u = 0
         } else {
               f = 2.0 * ( (1.0 - 6.0*x*x) * y*y * (1.0 - y*y)
                  + (1.0 - 6.0*y*y) * x*x * (1.0 - x*x) );
               ab[j][i] = hx * hy * f;
         }
      }
   }
   ierr = DMDAVecRestoreArray(da, b, &ab); CHKERRQ(ierr);
   return 0;
}
//ENDEXACT

