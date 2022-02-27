
static char help[] = "Tests VecView() contour plotting for 2d DMDAs.\n\n";

   /*


   VecView() on DMDA vectors first puts the Vec elements into global natural ordering before printing (or plotting)
   them. In 2d 5 by 2 DMDA this means the numbering is

      5  6   7   8   9
      0  1   2   3   4

   Now the default split across 2 processors with the DM  is (by rank)

      0  0   0  1   1
      0  0   0  1   1

   So the global PETSc ordering is

      3  4  5   8  9
      0  1  2   6  7

   Use the options
      -da_grid_x <nx> - number of grid points in x direction, if M < 0
      -da_grid_y <ny> - number of grid points in y direction, if N < 0
      -da_processors_x <MX> number of processors in x directio
      -da_processors_y <MY> number of processors in x direction
   */

#include <petscdm.h>
#include <petscdmda.h>

int main(int argc,char **argv)
{
   PetscMPIInt      rank;
   PetscInt         M = 10,N = 8;
   PetscErrorCode   ierr;
   PetscBool        flg = PETSC_FALSE;
   DM               da;
   PetscViewer      viewer;
   Vec              local,global;
   PetscScalar      value;
   DMBoundaryType   bx    = DM_BOUNDARY_NONE,by = DM_BOUNDARY_NONE;
   DMDAStencilType  stype = DMDA_STENCIL_BOX;

   ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;
   //   ierr = PetscViewerDrawOpen(PETSC_COMM_WORLD,0,"",300,0,300,300,&viewer);CHKERRQ(ierr);
   ierr = PetscViewerCreate(PETSC_COMM_WORLD, &viewer); // initialize the viewer object
      PetscViewerASCIIOpen(PETSC_COMM_WORLD,"view2d.out",&viewer);  //
   PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); // set viewer to print in matlab syntax


   ierr = PetscOptionsGetBool(NULL,NULL,"-star_stencil",&flg,NULL);CHKERRQ(ierr);
   if (flg) stype = DMDA_STENCIL_STAR;

   /* Create distributed array and get vectors */
   ierr = DMDACreate2d(PETSC_COMM_WORLD,bx,by,stype,M,N,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,&da);CHKERRQ(ierr);
   ierr = DMSetFromOptions(da);CHKERRQ(ierr);
   ierr = DMSetUp(da);CHKERRQ(ierr);
   ierr = DMCreateGlobalVector(da,&global);CHKERRQ(ierr);
   ierr = DMCreateLocalVector(da,&local);CHKERRQ(ierr);

   // The global vector is filled with -3
   // The local vectors are distributed
   value = -3.0;
   ierr  = VecSet(global,value);CHKERRQ(ierr);
   ierr  = DMGlobalToLocalBegin(da,global,INSERT_VALUES,local);CHKERRQ(ierr);
   ierr  = DMGlobalToLocalEnd(da,global,INSERT_VALUES,local);CHKERRQ(ierr);

   // Rank is the processer number 0, 1, 2, ..., Np-1
   // local value = 1, 2, 3, ..., Np
   ierr  = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRMPI(ierr);
   value = rank+1;

   ierr  = VecScale(local,value);CHKERRQ(ierr);
   ierr  = DMLocalToGlobalBegin(da,local,ADD_VALUES,global);CHKERRQ(ierr);
   ierr  = DMLocalToGlobalEnd(da,local,ADD_VALUES,global);CHKERRQ(ierr);

   flg  = PETSC_FALSE;
   ierr = PetscOptionsGetBool(NULL,NULL, "-view_global", &flg,NULL);CHKERRQ(ierr);
   if (flg) { /* view global vector in natural ordering */
      ierr = VecView(global,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   }
   ierr = DMView(da,viewer);CHKERRQ(ierr);
   // ierr = VecView(global,viewer);CHKERRQ(ierr);

   /* Free memory */
   ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
   ierr = VecDestroy(&local);CHKERRQ(ierr);
   ierr = VecDestroy(&global);CHKERRQ(ierr);
   ierr = DMDestroy(&da);CHKERRQ(ierr);
   ierr = PetscFinalize();
   return ierr;
}

   /*TEST

      test:
         requires: x
         nsize: 2
         args: -nox

   TEST*/
