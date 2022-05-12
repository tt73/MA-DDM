#include <petsc.h>
#include "MAfunctions.h"

/*
   This is the residue function for the 1D Monge-Ampere. Which is actually just
      det(D^2u) = u''.
   So we wish to find the root of
      F(u) = {f - u''  in the domain
             {u - g    on the boundary
   Since this is 1D, the ouptut F is represented as an array and so is the input u.
*/
PetscErrorCode MA1DFunctionLocal(DMDALocalInfo *info, PetscReal *u, PetscReal *F, MACtx *user) {
   PetscErrorCode ierr;
   PetscInt   i;
   PetscReal  xmax[1], xmin[1], h, x, ue, uw, f, temp;

   ierr = DMGetBoundingBox(info->da,xmin,xmax); CHKERRQ(ierr);
   h = (xmax[0] - xmin[0]) / (info->mx - 1);
   for (i = info->xs; i < info->xs + info->xm; i++) {
      x = xmin[0] + i * h;
      if (i==0 || i==info->mx-1) {
         // first and last points correspond to boundary
         user->g_bdry(x,0.0,0.0,user,&temp);
         F[i] = u[i] - temp;
      } else {
         // interior points depend on left and right neighbors
         uw = u[i-1];
         ue = u[i+1];
         if(i+1 == info->mx-1) {user->g_bdry(x+h,0.0,0.0,user,&temp); ue = temp;}
         if(i-1 == 0)          {user->g_bdry(x-h,0.0,0.0,user,&temp); uw = temp;}
         user->f_rhs(x,0.0,0.0,user,&f);
         F[i] = (-uw + 2.0*u[i] - ue)/(h*h) + f;
      }
   }
   ierr = PetscLogFlops(9.0*info->xm);CHKERRQ(ierr);
   return 0;
}


/*
   The equation we want to find the root of is F(u). F is the discretized version of the nonlinear system:
      F(u) = f - det(D^2(u)) in the domain and
      F(u) = u - g           on the boundary
   In the implementation, u is a 2D array and it's the 2nd arg. The F is also a 2D array and it's the 3rd arg.

*/
PetscErrorCode MA2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, MACtx *user) {
   PetscErrorCode ierr;
   PetscInt   i, j, k;
   PetscReal  xymin[2], xymax[2], hx, hy, x, y, ue, uw, un, us, temp;
   PetscReal  left, right; // left and right terms in the MA operator approximation
   PetscReal  weights[3], SDD[3]; // quadrature weight and second directional deriv
   ierr = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
   hx = (xymax[0] - xymin[0]) / (info->mx - 1);
   hy = (xymax[1] - xymin[1]) / (info->my - 1);

   // quad weights for order 1 approx
   weights[0] = PETSC_PI/4.0;
   weights[1] = PETSC_PI/2.0;
   weights[2] = PETSC_PI/4.0;

   for (j = info->ys; j < info->ys + info->ym; j++) {
      y = xymin[1] + j * hy;
      for (i = info->xs; i < info->xs + info->xm; i++) {
         x = xymin[0] + i * hx;
         if (i==0 || i==info->mx-1 || j==0 || j==info->my-1) {
            // on the boundary dOmega, F(u) = u - g
            // g_bdry(x,y,z) refers to the exact solution
            user->g_bdry(x,y,0.0,user,&temp);
            aF[j][i] = au[j][i] - temp;
         } else {
            un = au[j+1][i];
            us = au[j-1][i];
            ue = au[j][i+1];
            uw = au[j][i-1];
            if(j+1 == info->my-1) {user->g_bdry(x,y+hy,0.0,user,&temp); un = temp;}
            if(j-1 == 0)          {user->g_bdry(x,y-hy,0.0,user,&temp); us = temp;}
            if(i+1 == info->mx-1) {user->g_bdry(x+hx,y,0.0,user,&temp); ue = temp;}
            if(i-1 == 0)          {user->g_bdry(x-hx,y,0.0,user,&temp); uw = temp;}
            // in the interior Omega, F(u) = f - det+(D^2u)
            // There are only 2 directional derivs when
            SDD[0] = (uw - 2.0*au[j][i] + ue)/(hx*hx); // horizontal centered-diff
            SDD[1] = (un - 2.0*au[j][i] + us)/(hy*hy); // vertical centered-diff
            SDD[2] = (ue - 2.0*au[j][i] + uw)/(hx*hx); // horizontal centered-diff again

            // calculate the left term in the MA operator
            left = 0;
            for (k = 0; k < 3; k++) {
               left = left + weights[k]/PetscMax(SDD[k],user->epsilon);
            }
            left = PetscPowReal(left,-2.0);
            left = left * PETSC_PI*PETSC_PI;

            // calculate the right term in the MA operator
            right = user->epsilon;
            for (k = 0; k < 3; k++) {
               if (right > SDD[k]){ // find minimum
                  right = SDD[k];
               }
            }

            user->f_rhs(x,y,0.0,user,&temp);
            aF[j][i] = (left + right) - temp;
         }
      }
   }
   ierr = PetscLogFlops(11.0*info->xm*info->ym);CHKERRQ(ierr);
   return 0;
}
//ENDFORM2DFUNCTION

PetscErrorCode MA3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au, PetscReal ***aF, MACtx *user) {
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
   h = (xmax[0] - xmin[0]) / (info->mx - 1);
   for (i = info->xs; i < info->xs+info->xm; i++) { // loop over each row of J (mx by mx)
      row.i = i;
      col[0].i = i;
      ncols = 1;
      if (i==0 || i==info->mx-1) {
         v[0] = 1.0;  // if on the boundary, J_{i,j} = 1
      } else {
         v[0] = 2.0/(h*h); // middle J_{i,j} = 2/h^2
         if (i-1 > 0) {
            col[ncols].i = i-1;
            v[ncols++] = -1.0/(h*h); // left J_{i-1,j} = -1/h^2
         }
         if (i+1 < info->mx-1) {
            col[ncols].i = i+1;
            v[ncols++] = -1.0/(h*h); // right J_{i+1,j} = -1/h^2
         }
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
   PetscReal    xymin[2], xymax[2], x, y, hx, hy, h2, v[5];
   PetscReal    SDD[3], weights[3], u, ue, uw, un, us, temp;
   PetscBool    SGTE[3]; // stands for "SDD[k] greather than epsilon"
   PetscReal    common; // common factor for all stencil derivatives
   PetscInt     i,j,k,ncols;
   MatStencil   col[5],row;
   PetscReal    dDt[5][4];
   PetscInt     min_k; // jacobian of the min(min(D^2_theta,eps)) term depends on smallest D^2_theta

   ierr = DMGetBoundingBox(info->da,xymin,xymax); CHKERRQ(ierr);
   hx = (xymax[0] - xymin[0]) / (info->mx - 1);
   hy = (xymax[1] - xymin[1]) / (info->my - 1);
   h2 = hx*hy;

   // quad weights for order 1 approx
   weights[0] = PETSC_PI/4.0;
   weights[1] = PETSC_PI/2.0;
   weights[2] = PETSC_PI/4.0;

   // dDt[c][k] is the partial derivative of the D^2_theta_k at u_ij
   /*
      c = 0,1,2,3,4 corresponding to
      center, north, south, east, and west
   */
   dDt[0][0] = -2.0/h2; dDt[0][1] = -2.0/h2; dDt[0][2] = -2.0/h2; dDt[0][3] = 0;
   dDt[1][0] = 0;       dDt[1][1] = 1/h2;    dDt[1][2] = 0;       dDt[1][3] = 0;       // north
   dDt[2][0] = 0;       dDt[2][1] = 1/h2;    dDt[2][2] = 0;       dDt[2][3] = 0;       // south
   dDt[3][0] = 1/h2;    dDt[3][1] = 0;       dDt[3][2] = 1/h2;    dDt[3][3] = 0;       // east
   dDt[4][0] = 1/h2;    dDt[4][1] = 0;       dDt[4][2] = 1/h2;    dDt[4][3] = 0;       // west

   // loop over each row of J
   for (j = info->ys; j < info->ys+info->ym; j++) {
      row.j = j;
      col[0].j = j;
      y = xymin[1] + j * hy;
      // loop over each col of J
      for (i = info->xs; i < info->xs+info->xm; i++) {
         row.i = i;
         col[0].i = i;
         x = xymin[0] + i * hx;
         ncols = 1;

         if (i==0 || i==info->mx-1 || j==0 || j==info->my-1) {
            // Easy case: the Jacobian is 1 on the boundary
            v[0] = 1.0;
            ierr = MatSetValuesStencil(Jpre,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);

         } else {
            // Interior case: ... see write up
            // in the interior Omega, F(u) = f - det+(D^2u)
            u  = au[j][i];
            uw = au[j][i-1];
            ue = au[j][i+1];
            un = au[j+1][i];
            us = au[j-1][i];
            if(j+1 == info->my-1) {user->g_bdry(x,y+hy,0.0,user,&temp); un = temp;}
            if(j-1 == 0)          {user->g_bdry(x,y-hy,0.0,user,&temp); us = temp;}
            if(i+1 == info->mx-1) {user->g_bdry(x+hx,y,0.0,user,&temp); ue = temp;}
            if(i-1 == 0)          {user->g_bdry(x-hx,y,0.0,user,&temp); uw = temp;}

            // Directional derivs
            SDD[0] = (uw - 2.0*u + ue)/(hx*hx); // horizontal centered-diff
            SDD[1] = (un - 2.0*u + us)/(hy*hy); // vertical centered-diff
            SDD[2] = (ue - 2.0*u + uw)/(hx*hx); // horizontal centered-diff again

            // Comparison against SDD and espilon
            SGTE[0] = SDD[0] > user->epsilon;
            SGTE[1] = SDD[1] > user->epsilon;
            SGTE[2] = SDD[2] > user->epsilon;

            // Find the smallest SDD, then check it against eps
            min_k = 0;
            for(k=1; k<3; k++){
               if (SDD[min_k] > SDD[k]) min_k = k;
            }
            if (SDD[min_k] > user->epsilon){
               // if epsilon is the smallest Jacobian term is 0
               // dDt[.][3] = 0
               min_k = 3;
            }

            // Compute the common factor 2*pi^2*sum(w[k]/max(SDD[k],eps))^(-3)
            common = 0;
            for (k = 0; k < 3; k++) {
               common += weights[k]/(SGTE[k]? SDD[k]:user->epsilon);
            }
            common = PetscPowReal(common,-3.0);
            common *= 2*PETSC_PI*PETSC_PI;

            // Compute center value
            v[0] = 0;
            for (k=0; k<3; k++) {
               if(SGTE[k]) {
                  v[0] += weights[k]/(SDD[k]*SDD[k])*dDt[0][k];
               }
            }
            v[0] = common*v[0] + dDt[0][min_k]; // jacobian = (common term)*(summation term) + (min-min term)

            // Compute north value
            col[ncols].j = j+1;
            col[ncols].i = i;
            v[ncols] = 0;
            for (k = 0; k<3; k++) {
               if(SGTE[k]) {
                  v[ncols] += weights[k]/(SDD[k]*SDD[k])*dDt[1][k];
               }
            }
            v[ncols] = common*v[ncols] + dDt[1][min_k];
            ncols++;

            // Compute south value
            col[ncols].j = j-1;
            col[ncols].i = i;
            v[ncols] = 0;
            for (k = 0; k<3; k++) {
               if(SGTE[k]) {
                  v[ncols] += weights[k]/(SDD[k]*SDD[k])*dDt[2][k];
               }
            }
            v[ncols] = common*v[ncols] + dDt[2][min_k];
            ncols++;

            // east
            col[ncols].j = j;
            col[ncols].i = i+1;
            v[ncols] = 0;
            for (k = 0; k<3; k++) {
               if(SGTE[k]) {
                  v[ncols] += weights[k]/(SDD[k]*SDD[k])*dDt[3][k];
               }
            }
            v[ncols] = common*v[ncols] + dDt[3][min_k];
            ncols++;

            // west
            col[ncols].j = j;
            col[ncols].i = i-1;
            v[ncols] = 0;
            for (k = 0; k<3; k++) {
               if(SGTE[k]) {
                  v[ncols] += weights[k]/(SDD[k]*SDD[k])*dDt[4][k];
               }
            }
            v[ncols] = common*v[ncols] + dDt[4][min_k];
            ncols++;

            ierr = MatSetValuesStencil(Jpre,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
         }
      }
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
   The 3D Jacobian
*/
PetscErrorCode MA3DJacobianLocal(DMDALocalInfo *info, PetscScalar ***au, Mat J, Mat Jpre, MACtx *user) {

   return 0;
}

PetscErrorCode InitialState(DM da, InitialType it, PetscBool gbdry, Vec u, MACtx *user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscRandom    rctx;
    PetscReal temp;
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
                    user->g_bdry(x,0.0,0.0,user,&temp);
                    au[i] = temp;
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
                        user->g_bdry(x,y,0.0,user,&temp);
                        au[j][i] = temp;
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
                            user->g_bdry(x,y,z,user,&temp);
                            au[k][j][i] = temp;
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
