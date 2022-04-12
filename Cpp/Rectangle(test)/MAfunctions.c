#include <petsc.h>
#include "MAfunctions.h"

/*
   Note on DMDALocalInfo:
   PetscInt  mx,my,mz;    global number of grid points in each direction
   PetscInt  xs,ys,zs;    starting point of this processor, excluding ghosts
   PetscInt  xm,ym,zm;    number of grid points on this processor, excluding ghosts
   PetscInt  gxs,gys,gzs;    starting point of this processor including ghosts
   PetscInt  gxm,gym,gzm;    number of grid points on this processor including ghosts
*/

/*
   This is the residue function for the 1D Monge-Ampere. Which is actually just
      det(D^2u) = u''.
   So we wish to find the root of
      F(u) = {f - u''  in the domain
             {u - g    on the boundary
   Since this is 1D, the ouptut F is represented as an array and so is the input u.
*/
PetscErrorCode MA1DFunctionLocal(DMDALocalInfo *info, PetscReal *u, PetscReal *F, MACtx *user) {
   PetscInt   i;
   PetscReal  xmax[1], xmin[1], h, x, ue, uw;

   DMGetBoundingBox(info->da,xmin,xmax);
   h = (xmax[0] - xmin[0])/(info->mx + 1);
   for (i = info->xs; i < info->xs + info->xm; i++) {
      x = xmin[0] + (1+i)*h;
      ue = (i == info->mx-1) ? user->g_bdry(x+h,0.0,0.0,user) : u[i+1];
      uw = (i == 0)          ? user->g_bdry(x-h,0.0,0.0,user) : u[i-1];
      F[i] = (-uw + 2.0*u[i] - ue)/(h*h) + user->f_rhs(x,0.0,0.0,user);
   }
   return 0;
}


/*
   The equation we want to find the root of is F(u). F is the discretized version of the nonlinear system:
      F(u) =  det(D^2(u)) - f
   The 2nd arg **au is a 2D array representing the discretization of u i.e. au[i][j] ~ u(x_i,y_j)
   The 3rd arg **aF is a 2D array representing the discritization of F i.e. aF[i][j] ~ F(u(x_i,y_j))

   We use the approximation
                  ⎛2d + 1            ⎞-2
                  ⎜ ____             ⎟
                  ⎜ ╲                ⎟
                 2⎜  ╲       w_k     ⎟
     aF[i][j] = π ⎜  ╱   ──────────  ⎟  + min(min(D_k^2u, ε)) - f[i][j]
                  ⎜ ╱  max(D_k^2u, ε)⎟     k
                  ⎜ ‾‾‾‾               ⎟
                  ⎝ k = 0            ⎠
   where d is the stencil width, w_k is the quad weight, and D_k^2u is the kth 2nd directional derivative.
   See write-up for more details.
*/
PetscErrorCode MA2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, MACtx *user) {
   PetscErrorCode ierr;
   PetscInt   i, j, k;
   PetscReal  xymin[2], xymax[2], hx, hy, x, y;
   PetscReal  left, right; // left and right terms in the MA operator approximation
   PetscReal  *SDD; // second directional deriv
   PetscReal  *uFwd, *uBak; // u in the the forward and backward position for each direction k
   PetscReal  *hFwd, *hBak; // magnitude of step in forward and backward dirs for each direction k
   PetscInt   width, M; // stencil width, and total number of directions

   // get info from DA
   ierr = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
   hx   = (xymax[0] - xymin[0])/(info->mx + 1);
   hy   = (xymax[1] - xymin[1])/(info->my + 1);
   // allocate mem for derivative stuff
   width = info->sw;
   M     = 2*width + 1;
   PetscMalloc1(M,&SDD);          //
   PetscMalloc2(M,&uFwd,M,&uBak); //
   PetscMalloc2(M,&hFwd,M,&hBak); //
   // begin loop over all local interior nodes
   for (j = info->ys; j < info->ys + info->ym; j++) {
      y = xymin[1] + (j+1)*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = xymin[0] + (i+1)*hx;
         // Compute the SDDs
         if (width==1) {
            // width 1 case is hardcoded because it is simple
            uFwd[0] = (i == info->mx-1) ? user->g_bdry(x+hx,y,0.0,user) : au[j][i+1]; // east
            uBak[0] = (i == 0)          ? user->g_bdry(x-hx,y,0.0,user) : au[j][i-1]; // west
            uFwd[1] = (j == info->my-1) ? user->g_bdry(x,y+hy,0.0,user) : au[j+1][i]; // north
            uBak[1] = (j == 0)          ? user->g_bdry(x,y-hy,0.0,user) : au[j-1][i]; // south
            uFwd[2] = (i == 0)          ? user->g_bdry(x-hx,y,0.0,user) : au[j][i-1]; // west
            uBak[2] = (i == info->mx-1) ? user->g_bdry(x+hx,y,0.0,user) : au[j][i+1]; // east
            // 2nd dir. deriv
            SDD[0] = (uFwd[0] - 2.0*au[j][i] + uBak[0])/(hx*hx); // horizontal centered-diff
            SDD[1] = (uFwd[1] - 2.0*au[j][i] + uBak[1])/(hy*hy); // vertical centered-diff
            SDD[2] = (uFwd[2] - 2.0*au[j][i] + uBak[2])/(hx*hx); // horizontal centered-diff again
         } else if (width==2) { //
            // get fwd & bak u for north direction
            uFwd[width] = (j+2 > info->my-1) ? user->g_bdry(x,xymax[1],0.0,user) : au[j+2][i]; // north
            uBak[width] = (j-2 < 0)          ? user->g_bdry(x,xymin[1],0.0,user) : au[j-2][i]; // south
            hFwd[width] = (j+2 > info->my-1) ? hy : 2*hy; // north
            hBak[width] = (j-2 < 0)          ? hy : 2*hy; // south
            // get fwd & bak u for east and west
            uFwd[0]   = (i+2 > info->mx-1) ? user->g_bdry(xymax[0],y,0.0,user) : au[j][i+2];  // east
            uBak[0]   = (i-2 < 0)          ? user->g_bdry(xymin[0],y,0.0,user) : au[j][i-2];  // west
            uFwd[M-1] = uBak[0];
            uBak[M-1] = uFwd[0];
            hFwd[0]   = (i+2 > info->mx-1) ? hx : 2*hx;  // east
            hBak[0]   = (i-2 < 0)          ? hx : 2*hx;  // west
            hFwd[M-1] = hBak[0];
            hBak[M-1] = hFwd[0];
            // NE direction
            uFwd[1] = (i<info->mx-1 && j<info->my-1)? au[j+1][i+1] : user->g_bdry(x+hx,y+hy,0.0,user);
            hFwd[1] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // NW direction
            uFwd[M-2] = (i>0 && j<info->my-1)? au[j+1][i-1] : user->g_bdry(x-hx,y+hy,0.0,user);
            hFwd[M-2] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // SW
            uBak[1] = (i>0 && j>0)? au[j-1][i-1] : user->g_bdry(x-hx,y-hy,0.0,user);
            hBak[1] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // SE
            uBak[M-2] = (i<info->mx-1 && j>0)? au[j-1][i+1] : user->g_bdry(x+hx,y-hy,0.0,user);
            hBak[M-2] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // Compute SDD
            for (k=0; k<2*width+1; k++) {
               // Use formula for generalized centered difference
               SDD[k] = (hBak[k]*uFwd[k] - (hBak[k]+hFwd[k])*au[j][i] + hFwd[k]*uBak[k])/(hBak[k]*hFwd[k]*(hBak[k]+hFwd[k]));
            }
         } else { // Figure out how to compute SDD for general width

            // loop over all directions k
               // determmine uf ub
               // uf = u[j+dy][i+dx] if interior, else
               // use table(i,j,k) to get xp yp uf = g(xp,yp)

            // loop again
               // compute SDD

         }
         // calculate the left term in the MA operator
         left = 0;
         for (k=0; k<3; k++) {
            left = left + user->weights[k]/PetscMax(SDD[k],user->epsilon);
         }
         left = PetscPowReal(left,-2.0);
         left = left * PETSC_PI*PETSC_PI;
         // calculate the right term in the MA operator
         right = user->epsilon;
         for (k=0; k<3; k++) {
            if (right > SDD[k]){ // find minimum
               right = SDD[k];
            }
         }
         aF[j][i] = (left + right) - user->f_rhs(x,y,0.0,user);
      }
   }
   PetscFree(SDD);
   PetscFree2(uFwd,uBak);
   PetscFree2(hFwd,hBak);
   return 0;
}

// Place holder for 3D residue function
PetscErrorCode MA3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au, PetscReal ***aF, MACtx *user) {
   // PetscErrorCode ierr;
   // PetscInt   i, j, k;
   // PetscReal  xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag,
   //          x, y, z, ue, uw, un, us, uu, ud;
   // ierr = DMGetBoundingBox(info->da,xyzmin,xyzmax); CHKERRQ(ierr);
   // hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
   // hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
   // hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
   // dvol = hx * hy * hz;
   // scx = user->cx * dvol / (hx*hx);
   // scy = user->cy * dvol / (hy*hy);
   // scz = user->cz * dvol / (hz*hz);
   // scdiag = 2.0 * (scx + scy + scz);
   // for (k = info->zs; k < info->zs + info->zm; k++) {
   //    z = xyzmin[2] + k * hz;
   //    for (j = info->ys; j < info->ys + info->ym; j++) {
   //       y = xyzmin[1] + j * hy;
   //       for (i = info->xs; i < info->xs + info->xm; i++) {
   //          x = xyzmin[0] + i * hx;
   //          if (   i==0 || i==info->mx-1
   //             || j==0 || j==info->my-1
   //             || k==0 || k==info->mz-1) {
   //             aF[k][j][i] = au[k][j][i] - user->g_bdry(x,y,z,user);
   //             aF[k][j][i] *= scdiag;
   //          } else {
   //             ue = (i+1 == info->mx-1) ? user->g_bdry(x+hx,y,z,user) : au[k][j][i+1];
   //             uw = (i-1 == 0)          ? user->g_bdry(x-hx,y,z,user) : au[k][j][i-1];
   //             un = (j+1 == info->my-1) ? user->g_bdry(x,y+hy,z,user) : au[k][j+1][i];
   //             us = (j-1 == 0)          ? user->g_bdry(x,y-hy,z,user) : au[k][j-1][i];
   //             uu = (k+1 == info->mz-1) ? user->g_bdry(x,y,z+hz,user) : au[k+1][j][i];
   //             ud = (k-1 == 0)          ? user->g_bdry(x,y,z-hz,user) : au[k-1][j][i];
   //             aF[k][j][i] = scdiag * au[k][j][i]
   //                - scx * (uw + ue) - scy * (us + un) - scz * (uu + ud)
   //                - dvol * user->f_rhs(x,y,z,user);
   //          }
   //       }
   //    }
   // }
   // ierr = PetscLogFlops(14.0*info->xm*info->ym*info->zm);CHKERRQ(ierr);
   return 0;
}

/*
   The 1D Jacobian

   The au is the input of size mx.
   The output J is a matrix of size mx by mx.
   The values of J are filled exactly how you would expect...
      For each row, fill in 3 values
*/
PetscErrorCode MA1DJacobianLocal(DMDALocalInfo *info, PetscScalar *au, Mat J, Mat Jpre, MACtx *user) {
   PetscErrorCode  ierr;
   PetscInt     i,ncols;
   PetscReal    xmin[1], xmax[1], h, v[3];
   MatStencil   col[3],row;

   ierr = DMGetBoundingBox(info->da,xmin,xmax); CHKERRQ(ierr);
   h = (xmax[0]-xmin[0])/(info->mx+1);
   for (i = info->xs; i < info->xs+info->xm; i++) { // loop over each row of J (mx by mx)
      row.i = i;
      col[0].i = i;
      ncols = 1;

      v[0] = 2.0/(h*h); // middle J_{i,j} = 2/h^2
      if (i-1 > 0) {
         col[ncols].i = i-1;
         v[ncols++] = -1.0/(h*h); // left J_{i-1,j} = -1/h^2
      }
      if (i+1 < info->mx-1) {
         col[ncols].i = i+1;
         v[ncols++] = -1.0/(h*h); // right J_{i+1,j} = -1/h^2
      }
      // Insert up to 3 values of the Jacobian per row
      ierr = MatSetValuesStencil(Jpre,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
   }

   ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   if (J != Jpre) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   }
   return 0;
}

/*
   The 2D Jacobian

   The input au is a 2D array of size my by mx.
   The ouptut J is a matrix of size mx*my by mx*my.
   The J is constructed implicitly using stencil values.
*/
PetscErrorCode MA2DJacobianLocal(DMDALocalInfo *info, PetscScalar **au, Mat J, Mat Jpre, MACtx *user) {
   PetscErrorCode  ierr;
   PetscInt     i,j,k,di,dj,ncols,width,min_k,count,Nk,Ns; // Nk = # of directions, Ns = # of stencil pts
   PetscReal    xymin[2], xymax[2], x, y, hx, hy, h2;
   PetscReal    *v;
   MatStencil   *col;
   MatStencil   row;
   PetscReal    common;
   PetscReal    *uFwd, *uBak; // u in the the forward and backward position for each direction k
   PetscReal    *hFwd, *hBak; // magnitude of step in forward and backward dirs for each direction k
   PetscReal    *SDD;   // second directional deriv
   PetscBool    *SGTE;  // stands for "SDD[k] greather than epsilon"
   PetscReal    **dDt; // Jacobian of SDD operator for each direction and for the fwd, center, and bak

   // Retrieve info from DMDA
   ierr  = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
   hx    = (xymax[0] - xymin[0]) / (info->mx + 1);
   hy    = (xymax[1] - xymin[1]) / (info->my + 1);
   h2    = hx*hy;
   width = info->sw;
   Nk    = 2*width + 1; // number of angular forward directions
   Ns    = 4*width + 1; // total number of stencil points
   // Allocate memory
   PetscMalloc2(Ns,&v,Ns,&col);
   PetscCalloc2(Nk,&SDD,Nk,&SGTE);
   PetscCalloc2(Nk,&uFwd,Nk,&uBak);
   PetscCalloc2(Nk,&hFwd,Nk,&hBak);
   PetscMalloc1(Ns,&dDt);
   for (i=0; i<Ns; i++) {      // dDt is a 2D array, so
      PetscCalloc1(4,&dDt[i]); // mem allocated in a loop
   }
   /*
      dDt is the partial derivative of the SDD at a stencil point.
      dDt[i][j] represents the derivative wrt
   */
   if (width==1) {
      // Jacobian is independent of position for width=1 case
      dDt[0][0] = -2.0/h2; dDt[0][1] = -2.0/h2; dDt[0][2] = -2.0/h2; dDt[0][3] = 0; // center
      dDt[1][0] = 1/h2;    dDt[1][1] = 0;       dDt[1][2] = 1/h2;    dDt[1][3] = 0; // east
      dDt[2][0] = 1/h2;    dDt[2][1] = 0;       dDt[2][2] = 1/h2;    dDt[2][3] = 0; // west
      dDt[3][0] = 0;       dDt[3][1] = 1/h2;    dDt[3][2] = 0;       dDt[3][3] = 0; // north
      dDt[4][0] = 0;       dDt[4][1] = 1/h2;    dDt[4][2] = 0;       dDt[4][3] = 0; // south
   }
   // loop over each row of J
   for (j = info->ys; j < info->ys+info->ym; j++) {
      row.j = j;
      col[0].j = j;
      y = xymin[1] + (j+1)*hy;
      // loop over each col of J
      for (i = info->xs; i < info->xs+info->xm; i++) {
         row.i = i;
         col[0].i = i;
         x = xymin[0] + (i+1)*hx;
         ncols = 1;
         // Compute the
         if (width==1) {
            // width 1 case is hardcoded because it is simple
            uFwd[0] = (i == info->mx-1) ? user->g_bdry(x+hx,y,0.0,user) : au[j][i+1]; // east
            uBak[0] = (i == 0)          ? user->g_bdry(x-hx,y,0.0,user) : au[j][i-1]; // west
            uFwd[1] = (j == info->my-1) ? user->g_bdry(x,y+hy,0.0,user) : au[j+1][i]; // north
            uBak[1] = (j == 0)          ? user->g_bdry(x,y-hy,0.0,user) : au[j-1][i]; // south
            uFwd[2] = uBak[0];
            uBak[2] = uFwd[0];
            // 2nd dir. deriv
            SDD[0] = (uFwd[0] - 2.0*au[j][i] + uBak[0])/(hx*hx); // horizontal centered-diff
            SDD[1] = (uFwd[1] - 2.0*au[j][i] + uBak[1])/(hy*hy); // vertical centered-diff
            SDD[2] = (uFwd[2] - 2.0*au[j][i] + uBak[2])/(hx*hx); // horizontal centered-diff again
            // stencil step sizes
            hFwd[0] = hx; hBak[0] = hx;
            hFwd[1] = hy; hBak[1] = hy;
            hFwd[2] = hx; hBak[2] = hx;
         } else if (width==2) {
            // get fwd & bak u for north direction
            uFwd[width] = (j+2 > info->my-1) ? user->g_bdry(x,xymax[1],0.0,user) : au[j+2][i]; // north
            uBak[width] = (j-2 < 0)          ? user->g_bdry(x,xymin[1],0.0,user) : au[j-2][i]; // south
            hFwd[width] = (j+2 > info->my-1) ? hy : 2*hy; // north
            hBak[width] = (j-2 < 0)          ? hy : 2*hy; // south
            // get fwd & bak u for east and west
            uFwd[0]   = (i+2 > info->mx-1) ? user->g_bdry(xymax[0],y,0.0,user) : au[j][i+2];  // east
            uBak[0]   = (i-2 < 0)          ? user->g_bdry(xymin[0],y,0.0,user) : au[j][i-2];  // west
            uFwd[Nk-1] = uBak[0];
            uBak[Nk-1] = uFwd[0];
            hFwd[0]   = (i+2 > info->mx-1) ? hx : 2*hx;  // east
            hBak[0]   = (i-2 < 0)          ? hx : 2*hx;  // west
            hFwd[Nk-1] = hBak[0];
            hBak[Nk-1] = hFwd[0];
            // NE direction
            uFwd[1] = (i<info->mx-1 && j<info->my-1)? au[j+1][i+1] : user->g_bdry(x+hx,y+hy,0.0,user);
            hFwd[1] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // NW direction
            uFwd[Nk-2] = (i>0 && j<info->my-1)? au[j+1][i-1] : user->g_bdry(x-hx,y+hy,0.0,user);
            hFwd[Nk-2] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // SW
            uBak[1] = (i>0 && j>0)? au[j-1][i-1] : user->g_bdry(x-hx,y-hy,0.0,user);
            hBak[1] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // SE
            uBak[Nk-2] = (i<info->mx-1 && j>0)? au[j-1][i+1] : user->g_bdry(x+hx,y-hy,0.0,user);
            hBak[Nk-2] = PetscSqrtReal(hx*hx+hy*hy); // doesn't change for width=2
            // Compute SDD
            for (k=0; k<2*width+1; k++) {
               // Use formula for generalized centered difference
               SDD[k] = (hBak[k]*uFwd[k] - (hBak[k]+hFwd[k])*au[j][i] + hFwd[k]*uBak[k])/(hBak[k]*hFwd[k]*(hBak[k]+hFwd[k]));
            }
         }
         // Comparison against SDD and espilon
         for (k=0; k<Nk; k++) {
            SGTE[k] = SDD[k] > user->epsilon;
         }
         // Find the smallest SDD, then check it against eps
         min_k = 0;
         for(k=1; k<Nk; k++) {
            if (SDD[min_k] > SDD[k]) min_k = k;
         }
         if (SDD[min_k]>user->epsilon) {
            min_k = 3;
         }
         // Compute the common factor 2*pi^2*sum(w[k]/max(SDD[k],eps))^(-3)
         common = 0;
         for (k=0; k<Nk; k++) {
            common += user->weights[k]/(SGTE[k]? SDD[k]:user->epsilon);
         }
         common = PetscPowReal(common,-3.0);
         common *= 2*PETSC_PI*PETSC_PI;
         // Compute center value
         v[0] = 0;
         for (k=0; k<Nk; k++) {
            if(SGTE[k]) {
               // v[0] += user->weights[k]/(SDD[k]*SDD[k])*dDt[0][k];
               v[0] += user->weights[k]/(SDD[k]*SDD[k])*(-2.0/(hFwd[k]*hBak[k]));
            }
         }
         v[0] = common*v[0] + dDt[0][min_k]; // jacobian = (common term)*(summation term) + (min-min term)
         if (width==1) {
            // Compute other values for remaining 4*width stencil pts
            count = 1;
            for (di=width; di>-width; di--) {
               dj = -PetscAbsReal(di)+width;
               // forward
               if (i+di<=info->mx-1 && j+dj<=info->my-1) {
                  col[ncols].j = j+dj;
                  col[ncols].i = i+di;
                  v[ncols] = 0;
                  for (k = 0; k<Nk; k++) {
                     if(SGTE[k]) {
                        v[ncols] += user->weights[k]/(SDD[k]*SDD[k])*dDt[count][k];
                     }
                  }
                  v[ncols] = common*v[ncols] + dDt[count][min_k];
                  ncols++;
               }
               // backward
               if (i-di>=0 && j-dj>=0) {
                  col[ncols].j = j-dj;
                  col[ncols].i = i-di;
                  v[ncols] = 0;
                  for (k = 0; k<Nk; k++) {
                     if(SGTE[k]) {
                        v[ncols] += user->weights[k]/(SDD[k]*SDD[k])*dDt[count+1][k];
                     }
                  }
                  v[ncols] = common*v[ncols] + dDt[count+1][min_k];
                  ncols++;
               }
               count += 2;
            }
         }
         // Insert values
         PetscPrintf(PETSC_COMM_WORLD, "Now inserting values for i,j = (%d,%d)\n",i,j);
         for (k=0; k<ncols; k++) {
            PetscPrintf(PETSC_COMM_WORLD, "v[%d][%d] = %f\n",col[k].i,col[k].j,v[k]);
         }
         ierr = MatSetValuesStencil(Jpre,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
      }
   }

   ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   if (J != Jpre) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   }
   // Free memory
   PetscFree2(SDD,SGTE);
   PetscFree2(uFwd,uBak);
   PetscFree2(hFwd,hBak);
   return 0;
}

// Placeholder for 3D Jacobian
PetscErrorCode MA3DJacobianLocal(DMDALocalInfo *info, PetscScalar ***au, Mat J, Mat Jpre, MACtx *user) {
   // PetscErrorCode  ierr;
   // PetscReal   xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag, v[7];
   // PetscInt    i,j,k,ncols;
   // MatStencil  col[7],row;

   // ierr = DMGetBoundingBox(info->da,xyzmin,xyzmax); CHKERRQ(ierr);
   // hx = (xyzmax[0] - xyzmin[0]) / (info->mx - 1);
   // hy = (xyzmax[1] - xyzmin[1]) / (info->my - 1);
   // hz = (xyzmax[2] - xyzmin[2]) / (info->mz - 1);
   // dvol = hx * hy * hz;
   // scx = user->cx * dvol / (hx*hx);
   // scy = user->cy * dvol / (hy*hy);
   // scz = user->cz * dvol / (hz*hz);
   // scdiag = 2.0 * (scx + scy + scz);
   // for (k = info->zs; k < info->zs+info->zm; k++) {
   //    row.k = k;
   //    col[0].k = k;
   //    for (j = info->ys; j < info->ys+info->ym; j++) {
   //       row.j = j;
   //       col[0].j = j;
   //       for (i = info->xs; i < info->xs+info->xm; i++) {
   //          row.i = i;
   //          col[0].i = i;
   //          ncols = 1;
   //          v[0] = scdiag;
   //          if (i>0 && i<info->mx-1 && j>0 && j<info->my-1 && k>0 && k<info->mz-1) {
   //             if (i-1 > 0) {
   //                col[ncols].k = k;    col[ncols].j = j;    col[ncols].i = i-1;
   //                v[ncols++] = - scx;
   //             }
   //             if (i+1 < info->mx-1) {
   //                col[ncols].k = k;    col[ncols].j = j;    col[ncols].i = i+1;
   //                v[ncols++] = - scx;
   //             }
   //             if (j-1 > 0) {
   //                col[ncols].k = k;    col[ncols].j = j-1;  col[ncols].i = i;
   //                v[ncols++] = - scy;
   //             }
   //             if (j+1 < info->my-1) {
   //                col[ncols].k = k;    col[ncols].j = j+1;  col[ncols].i = i;
   //                v[ncols++] = - scy;
   //             }
   //             if (k-1 > 0) {
   //                col[ncols].k = k-1;  col[ncols].j = j;    col[ncols].i = i;
   //                v[ncols++] = - scz;
   //             }
   //             if (k+1 < info->mz-1) {
   //                col[ncols].k = k+1;  col[ncols].j = j;    col[ncols].i = i;
   //                v[ncols++] = - scz;
   //             }
   //          }
   //          ierr = MatSetValuesStencil(Jpre,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
   //       }
   //    }
   // }
   // ierr = MatAssemblyBegin(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   // ierr = MatAssemblyEnd(Jpre,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   // if (J != Jpre) {
   //    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   //    ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
   // }
   return 0;
}

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry, Vec u, MACtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscRandom    rctx;
    switch (it) {
        case ZEROS:
            ierr = VecSet(u,0.0); CHKERRQ(ierr);
            break;
        case RANDOM:
            ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx); CHKERRQ(ierr);
            ierr = VecSetRandom(u,rctx); CHKERRQ(ierr);
            ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
            break;
        default:
            SETERRQ(PETSC_COMM_SELF,4,"invalid InitialType ... how did I get here?\n");
    }
    if (!gbdry) {
        return 0;
    }
    ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
    switch (info.dim) {
        case 1:
        {
            PetscInt  i;
            PetscReal xmax[1], xmin[1], h, x, *au;
            ierr = DMDAVecGetArray(da, u, &au); CHKERRQ(ierr);
            ierr = DMGetBoundingBox(da,xmin,xmax); CHKERRQ(ierr);
            h = (xmax[0] - xmin[0]) / (info.mx - 1);
            for (i = info.xs; i < info.xs + info.xm; i++) {
                if (i==0 || i==info.mx-1) {
                    x = xmin[0] + i * h;
                    au[i] = user->g_bdry(x,0.0,0.0,user);
                }
            }
            ierr = DMDAVecRestoreArray(da, u, &au); CHKERRQ(ierr);
            break;
        }
        case 2:
        {
            PetscInt   i, j;
            PetscReal  xymin[2], xymax[2], hx, hy, x, y, **au;
            ierr = DMDAVecGetArray(da, u, &au); CHKERRQ(ierr);
            ierr = DMGetBoundingBox(da,xymin,xymax); CHKERRQ(ierr);
            hx = (xymax[0] - xymin[0]) / (info.mx - 1);
            hy = (xymax[1] - xymin[1]) / (info.my - 1);
            for (j = info.ys; j < info.ys + info.ym; j++) {
                y = xymin[1] + j * hy;
                for (i = info.xs; i < info.xs + info.xm; i++) {
                    if (i==0 || i==info.mx-1 || j==0 || j==info.my-1) {
                        x = xymin[0] + i * hx;
                        au[j][i] = user->g_bdry(x,y,0.0,user);
                    }
                }
            }
            ierr = DMDAVecRestoreArray(da, u, &au); CHKERRQ(ierr);
            break;
        }
        case 3:
        {
            PetscInt   i, j, k;
            PetscReal  xyzmin[3], xyzmax[3], hx, hy, hz, x, y, z, ***au;
            ierr = DMDAVecGetArray(da, u, &au); CHKERRQ(ierr);
            ierr = DMGetBoundingBox(da,xyzmin,xyzmax); CHKERRQ(ierr);
            hx = (xyzmax[0] - xyzmin[0]) / (info.mx - 1);
            hy = (xyzmax[1] - xyzmin[1]) / (info.my - 1);
            hz = (xyzmax[2] - xyzmin[2]) / (info.mz - 1);
            for (k = info.zs; k < info.zs+info.zm; k++) {
                z = xyzmin[2] + k * hz;
                for (j = info.ys; j < info.ys + info.ym; j++) {
                    y = xyzmin[1] + j * hy;
                    for (i = info.xs; i < info.xs + info.xm; i++) {
                        if (i==0 || i==info.mx-1 || j==0 || j==info.my-1
                                 || k==0 || k==info.mz-1) {
                            x = xyzmin[0] + i * hx;
                            au[k][j][i] = user->g_bdry(x,y,z,user);
                        }
                    }
                }
            }
            ierr = DMDAVecRestoreArray(da, u, &au); CHKERRQ(ierr);
            break;
        }
        default:
            SETERRQ(PETSC_COMM_SELF,5,"invalid dim from DMDALocalInfo\n");
    }
    return 0;
}

PetscErrorCode ComputeWeights(PetscInt width, PetscInt order, MACtx *user) {
   PetscInt   i,M;
   PetscReal  a,h1,h2;
   PetscReal  *theta;

   // Compute angles of L1 stencil points
   M = 2*width+1;
   PetscMalloc1(M,&theta);
   theta[width] = PETSC_PI/2;
   for (i=0; i<width; i++) {
      a = PetscAtanReal((i)/(width-i));
      theta[i]     = a;
      theta[M-i-1] = PETSC_PI - a;
   }
   // Compute quadrature weights
   PetscCalloc1(M,&user->weights);
   if (order==1) {
      // Irregularly spaced trapezoid
      if (width==1) {
         user->weights[0] = PETSC_PI/4.0;
         user->weights[1] = PETSC_PI/2.0;
         user->weights[2] = PETSC_PI/4.0;
      } else { // for width larger than 1
         user->weights[0] = (theta[1]-theta[0])/2.0;
         for (i=1; i<M-1; i++){
            user->weights[i] = (theta[i+1]-theta[i-1])/2.0;
         }
         user->weights[M-1] = (theta[M-1]-theta[M-2])/2.0;
      }
   } else if(order==2) {
      // Irregularly spaced Simpson's
      for(i=0; i<width; i++) {
         h1 = theta[2*i+1] - theta[2*i];
         h2 = theta[2*i+2] - theta[2*i+1];
         a = (h1+h2)/6.0;
         user->weights[2*i]   += a*(2-h2/h1);
         user->weights[2*i+1] += a*(h1+h2)*(h1+h2)/(h1*h2);
         user->weights[2*i+2] += a*(2-h1/h2);
      }
   } else { // for all other orders
      SETERRQ(PETSC_COMM_SELF,5,"Quadarature order > 2 not supported.\n");
   }
   PetscFree(theta);
   return 0;
}