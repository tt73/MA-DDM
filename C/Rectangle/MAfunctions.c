#include <petsc.h>
#include "MAfunctions.h"


/* Initial State - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   This function sets up the intial guess for the outer SNES method.
   The available options are: zeros, random, corner, and pyramid.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
PetscErrorCode InitialState(DM da, InitialType it, Vec u, MACtx *user) {
   PetscErrorCode ierr;
   DMDALocalInfo  info;
   PetscRandom    rctx;
   PetscReal      temp,temp1,temp2;
   PetscFunctionBeginUser;

   switch (it) {
      case ZEROS: // just set u = 0
         ierr = VecSet(u,0.0); CHKERRQ(ierr);
         break;
      case RANDOM:
         ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx); CHKERRQ(ierr);
         ierr = VecSetRandom(u,rctx); CHKERRQ(ierr);
         ierr = PetscRandomDestroy(&rctx); CHKERRQ(ierr);
         break;
      case CORNER:
         /* Corner - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Get largest Dirichlet data on all corners.
            Set the initial guess to that constant.
            We may want to experiment with the smallest values of the 4 corners.
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
         ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
         switch (info.dim) {
            case 1:
            {
               user->g_bdry(user->xmin,0.0,0.0,user,&temp1);
               user->g_bdry(user->xmax,0.0,0.0,user,&temp2);
               temp = PetscMax(temp1,temp2);
               break;
            }
            case 2:
            {
               user->g_bdry(user->xmin,user->ymin,0.0,user,&temp1);
               user->g_bdry(user->xmax,user->ymin,0.0,user,&temp2);
               temp = PetscMax(temp1,temp2);
               user->g_bdry(user->xmin,user->ymax,0.0,user,&temp1);
               user->g_bdry(user->xmax,user->ymax,0.0,user,&temp2);
               temp = PetscMax(temp,PetscMax(temp1,temp2));
               break;
            }
            case 3:
            {
               SETERRQ(PETSC_COMM_SELF,3,"dim = 3 not yet supported\n");
               break;
            }
            default:
               SETERRQ(PETSC_COMM_SELF,5,"invalid dim from DMDALocalInfo\n");

         }
         ierr = VecSet(u,temp); CHKERRQ(ierr);
         break;
      case PYRAMID: // needs fixing
         /*
            Convex pyramid:

            In 1D we have, a < x < b and u(a) = ga and u(b) = gb.
            We want the initial guess to be two line segments connecting
            (a,ga) to ((a+b)/2, gm ) to (b,gb). The midpoint is slightly
            lower than ga and gb. We will let gm = min(ga,gb) - (ga+gb)/2.
            This is equivalent to letting u0(x) = max(ma*(x-a)+ga,mb*(x-b)+gb)
            where the slopes are given by ma = 2(gm-ga)/(b-a), and
            mb = 2(gb-gm)/(b-a).

            In 1D, u0 is shaped like a V. We can generalize this shape to a convex
            pyramid in 2D. We let u0 be the max of 4 lines:
            u0(x,y) = max( mxa*(x-xmin)+gxa,
                           mxb*(x-xmax)+gxb,
                           mya*(y-ymin)+gya,
                           myb*(y-ymax)+gyb)
            In 2D, we let gm = min(gxa,gxb,gya,gyb) - avg(gxa,gxb,gya,gyb).
         */
         ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
         switch (info.dim) {
            case 1: // 1D V-shape
            {
               PetscInt  i;
               PetscReal hx,x,*au,ga,gb,gm,ma,mb;

               ierr = DMDAVecGetArray(da, u, &au); CHKERRQ(ierr);
               hx = (user->xmax - user->xmin)/(PetscReal)(info.mx + 1);
               user->g_bdry(user->xmin,0.0,0.0,user,&ga);
               user->g_bdry(user->xmax,0.0,0.0,user,&gb);
               gm = PetscMin(ga,gb) - (ga+gb)/2.0;
               ma = 2.0*(gm-ga)/(user->xmax-user->xmin);
               mb = 2.0*(gb-gm)/(user->xmax-user->xmin);
               for (i=info.xs; i<info.xs+info.xm; i++) {
                  x = user->xmin + (i+1)*hx;
                  au[i] = PetscMax(ma*(x-user->xmin)+ga,mb*(x-user->xmax)+gb);
               }
               ierr = DMDAVecRestoreArray(da, u, &au); CHKERRQ(ierr);
               break;
            }
            case 2: // 2D pyramid-shape
            {
               PetscInt   i, j;
               PetscReal  hx,hy,x,y,**au,gxa,gxb,gya,gyb,gmin,gm,mxa,mxb,mya,myb,temp;

               ierr = DMDAVecGetArray(da, u, &au); CHKERRQ(ierr);
               hx = (user->xmax - user->xmin)/(PetscReal)(info.mx + 1);
               hy = (user->ymax - user->ymin)/(PetscReal)(info.my + 1);
               for (j = info.ys; j < info.ys + info.ym; j++) {
                  y = user->ymin + (j+1)*hy;
                  user->g_bdry(user->xmin,y,0.0,user,&gya);
                  user->g_bdry(user->xmax,y,0.0,user,&gyb);
                  for (i = info.xs; i < info.xs+info.xm; i++) {
                     x    = user->xmin + (i+1)*hx;
                     temp = PetscMin(gya,gyb);
                     user->g_bdry(x,user->ymin,0.0,user,&gxa);
                     user->g_bdry(x,user->ymax,0.0,user,&gxb);
                     gmin = PetscMin(temp,PetscMin(gxa,gxb));
                     gm   = gmin - (gxa+gxb+gya+gyb)/4.0;
                     // gxm = PetscMin(gxa,gxb) - (gxa+gxb)/2.0;
                     mya = 2.0*(gm-gya)/(user->xmax - user->xmin);
                     myb = 2.0*(gyb-gm)/(user->xmax - user->xmin);
                     mxa = 2.0*(gm-gxa)/(user->ymax - user->ymin);
                     mxb = 2.0*(gxb-gm)/(user->ymax - user->ymin);
                     temp = PetscMax(mya*(y-user->ymin)+gya,myb*(y-user->ymax)+gyb);
                     au[j][i] = PetscMax(temp,PetscMax(mxa*(x-user->xmin)+gxa,mxb*(x-user->xmax)+gxb));
                  }
               }
               ierr = DMDAVecRestoreArray(da, u, &au); CHKERRQ(ierr);
               break;
            }
            case 3:
            {
               SETERRQ(PETSC_COMM_SELF,3,"dim = 3 not yet supported\n");
               break;
            }
            default:
               SETERRQ(PETSC_COMM_SELF,5,"invalid dim from DMDALocalInfo\n");
         }
         break;
      default:
         SETERRQ(PETSC_COMM_SELF,4,"invalid InitialType ... how did I get here?\n");
   }
   PetscFunctionReturn(0);
}


/* 1D Version of the Residue Function

   This is the residue function for the 1D Monge-Ampere. Which is actually just
      det(D^2u) = u''.
   So we wish to find the root of
      F(u) = {f - u''  in the domain
             {u - g    on the boundary
   Since this is 1D, the ouptut F is represented as an array and so is the input u.
*/
PetscErrorCode MA1DFunctionLocal(DMDALocalInfo *info, PetscReal *u, PetscReal *F, MACtx *user) {
   PetscInt  i;
   PetscReal h, x, ue, uw, f;
   PetscFunctionBeginUser;

   h = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   for (i=info->xs; i<info->xs+info->xm; i++) {
      x = user->xmin + (i+1)*h;
      ue = u[i+1];
      uw = u[i-1];
      if(i == info->mx-1) {user->g_bdry(x+h,0.0,0.0,user,&ue);}
      if(i == 0)          {user->g_bdry(x-h,0.0,0.0,user,&uw);}
      user->f_rhs(x,0.0,0.0,user,&f);
      F[i] = (-uw + 2.0*u[i] - ue)/(h*h) + f;
   }
   PetscFunctionReturn(0);
}

/* 2D Version of the Residue Function

   The equation we want to find the root of is F(u). F is the discretized version of the nonlinear system:
      F(u) =  det(D^2(u)) - f
   The 2nd arg **au is a 2D array representing the discretization of u i.e. au[i][j] ~ u(x_i,y_j)
   The 3rd arg **aF is a 2D array representing the discritization of F i.e. aF[i][j] ~ F(u(x_i,y_j))
*/
PetscErrorCode MA2DFunctionLocal(DMDALocalInfo *info, PetscReal **au, PetscReal **aF, MACtx *user) {
   PetscInt   i,j;
   PetscReal  hx,hy,x,y,f;
   PetscReal  DetD2u; // MA operator approximation
   PetscReal *SDD; // second directional deriv
   PetscReal *hFwd, *hBak; // magnitude of step in forward and backward dirs for each direction k
   PetscInt   width, M; // stencil width, and total number of directions
   PetscFunctionBeginUser;

   // get info from DA
   hx = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   hy = (user->ymax-user->ymin)/(PetscReal)(info->my+1);
   // allocate mem for derivative stuff
   width = info->sw;
   M     = 2*width;
   PetscMalloc1(M,&SDD);
   PetscMalloc2(M,&hFwd,M,&hBak);
   // begin loop over all local interior nodes
   for (j=info->ys; j<info->ys+info->ym; j++) {
      y = user->ymin + (j+1)*hy;
      for (i=info->xs; i<info->xs+info->xm; i++) {
         x = user->xmin + (i+1)*hx;
         ComputeSDD(info,au,user,i,j,x,y,SDD,hFwd,hBak);
         ApproxDetD2u(&DetD2u,M,SDD,user);
         user->f_rhs(x,y,0.0,user,&f);
         aF[j][i] = DetD2u - f;
      }
   }
   PetscFree(SDD);
   PetscFree2(hFwd,hBak);
   PetscFunctionReturn(0);
}

// Place holder for 3D residue function
PetscErrorCode MA3DFunctionLocal(DMDALocalInfo *info, PetscReal ***au, PetscReal ***aF, MACtx *user) {
   // PetscErrorCode ierr;
   // PetscInt   i, j, k;
   // PetscReal  xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag,
   //          x, y, z, ue, uw, un, us, uu, ud;
   // PetscFunctionBeginUser;
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
   PetscFunctionReturn(0);
}

/* 1D Jacobian

   The au is the input of size mx.
   The output J is a matrix of size mx by mx.
   The values of J are filled exactly how you would expect...
      For each row, fill in 3 values
*/
PetscErrorCode MA1DJacobianLocal(DMDALocalInfo *info, PetscScalar *au, Mat J, Mat Jpre, MACtx *user) {
   PetscErrorCode  ierr;
   PetscInt        i,ncols;
   PetscReal       h, v[3];
   MatStencil      col[3],row;
   PetscFunctionBeginUser;

   h  = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   for (i=info->xs; i<info->xs+info->xm; i++) { // loop over each row of J (mx by mx)
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
   PetscFunctionReturn(0);
}

/* 2D Jacobian

   The input au is a 2D array of size my by mx.
   The ouptut J is a matrix of size mx*my by mx*my.
   The J is constructed implicitly using stencil values.
*/
PetscErrorCode MA2DJacobianLocal(DMDALocalInfo *info, PetscScalar **au, Mat J, Mat Jpre, MACtx *user) {
   PetscErrorCode ierr;
   PetscInt       i,j,k,Si,Sj,ncols,width,min_k,Nk,Ns,factor; // Nk = # of directions, Ns = # of stencil pts
   PetscReal      x,y,hx,hy,temp,dSDD;
   PetscReal      *v;
   MatStencil     *col;
   MatStencil     row;
   PetscReal      common;
   PetscReal      *hFwd,*hBak; // magnitude of step in forward and backward dirs for each direction k
   PetscReal      *SDD;        // second directional deriv
   PetscInt       *SGTE;       // SDD is greater than epsilon
   PetscBool      regularize;  // true if epsilon is the smallest among SDD
   PetscFunctionBeginUser;

   hx = (user->xmax-user->xmin)/(PetscReal)(info->mx+1);
   hy = (user->ymax-user->ymin)/(PetscReal)(info->my+1);
   width = info->sw;
   Nk    = 2*width; // number of angular forward directions
   Ns    = 4*width+1; // total number of stencil points
   // Allocate memory
   PetscCalloc2(Ns,&v,Ns,&col);
   PetscCalloc2(Nk,&SDD,Nk,&SGTE);
   PetscCalloc2(Nk,&hFwd,Nk,&hBak);
   // loop over each row of J
   for (j = info->ys; j < info->ys+info->ym; j++) {
      row.j = j;
      col[0].j = j;
      y = user->ymin + (j+1)*hy;
      // loop over each col of J
      for (i = info->xs; i < info->xs+info->xm; i++) {
         row.i = i;
         col[0].i = i;
         x = user->xmin + (i+1)*hx;
         ncols = 1;
         ComputeSDD(info,au,user,i,j,x,y,SDD,hFwd,hBak);
         for (k=0; k<Nk; k++) {
            SGTE[k] = SDD[k] > user->epsilon;
         }
         // Loop to find the smallest SDD
         min_k = 0;
         for(k=1; k<Nk; k++) {
            if (SDD[min_k] > SDD[k]) {
               min_k = k;
            }
         }
         // then check it smallest against eps
         if (SDD[min_k]>user->epsilon) {
            min_k = Nk+1;
            regularize = PETSC_TRUE; // if true, derivative of min-min term is 0
         } else {
            regularize = PETSC_FALSE; // if false,
         }
         // Compute the common factor 2*pi^2*sum(w[k]/max(SDD[k],eps))^(-3)
         common = 0;
         for (k=0; k<Nk; k++) {
            factor = (k==0)? 2:1;
            common += factor*(user->weights[k])/(SGTE[k]? SDD[k]:user->epsilon);
         }
         common = 2.0*PETSC_PI*PETSC_PI*PetscPowReal(common,-3.0);
         // Compute center value
         v[0] = 0;
         for (k=0; k<Nk; k++) {
            if(SGTE[k]) {
               factor = (k==0)? 2:1;
               dSDD = -2.0/(hFwd[k]*hBak[k]); // derivative of SDD
               v[0] += factor*(user->weights[k])/(SDD[k]*SDD[k])*dSDD;
            }
         }
         v[0] *= common;
         v[0] += (regularize)? 0 : -2.0/(hFwd[min_k]*hBak[min_k]);
         // Compute partial at remaining 4*width stencil points
         for (k=0; k<Nk; k++) {
            Si = user->Si[k];
            Sj = user->Sj[k];
            factor = (k==0)? 2:1;
            if (i+Si>=0 && i+Si<=info->mx-1 && j+Sj<=info->my-1 && j+Sj>=0) {
               col[ncols].j = j+Sj;
               col[ncols].i = i+Si;
               v[ncols] = 0;
               if(SGTE[k]) {
                  temp = 2.0/(hFwd[k]*(hFwd[k]+hBak[k]));
                  v[ncols] = factor*(user->weights[k])/(SDD[k]*SDD[k])*temp;
               }
               v[ncols] *= common;
               if (!regularize && min_k==k) {
                  v[ncols] += (2.0/(hFwd[min_k]*(hFwd[min_k]+hBak[min_k])));
               }
               ncols++;
            }
            if (i-Si>=0 && i-Si<=info->mx-1 && j-Sj>=0 && j-Sj<=info->my-1) {
               col[ncols].j = j-Sj;
               col[ncols].i = i-Si;
               v[ncols] = 0;
               if(SGTE[k]) {
                  dSDD = 2.0/(hBak[k]*(hFwd[k]+hBak[k]));
                  v[ncols] = factor*(user->weights[k])/(SDD[k]*SDD[k])*dSDD;
               }
               v[ncols] *= common;
               if (!regularize && min_k==k) {
                  v[ncols] += (2.0/(hBak[min_k]*(hFwd[min_k]+hBak[min_k])));
               }
               ncols++;
            }
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
   PetscFree2(hFwd,hBak);
   PetscFunctionReturn(0);
}

// Placeholder for 3D Jacobian
PetscErrorCode MA3DJacobianLocal(DMDALocalInfo *info, PetscScalar ***au, Mat J, Mat Jpre, MACtx *user) {
   // PetscErrorCode  ierr;
   // PetscReal   xyzmin[3], xyzmax[3], hx, hy, hz, dvol, scx, scy, scz, scdiag, v[7];
   // PetscInt    i,j,k,ncols;
   // MatStencil  col[7],row;
   // PetscFunctionBeginUser;
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
   PetscFunctionReturn(0);
}

/* Compute the approximation of the det using the SDD.
   We use the approximation
                  ⎛2d + 1            ⎞-2
                  ⎜ ____             ⎟
                  ⎜ ╲                ⎟
                 2⎜  ╲       w_k     ⎟
     det(D2u) ~ π ⎜  ╱   ──────────  ⎟  + min(min(D_k^2u, ε))
                  ⎜ ╱  max(D_k^2u, ε)⎟     k
                  ⎜ ‾‾‾‾               ⎟
                  ⎝ k = 0            ⎠
   where d is the stencil width, w_k is the quad weight, and D_k^2u is the kth 2nd directional derivative.
   See write-up for more details. The summation term is called `left` and the min-min term is called `right`
   in the code.

   Note: there are 2d+1 directions but the last SDD is the same as the last. Therefore we only sum over 2*d
   terms and the first term is added twice.
*/
PetscErrorCode ApproxDetD2u(PetscReal *DetD2u, PetscInt dim, PetscReal *SDD, MACtx *user) {
   PetscInt  k;
   PetscReal left,right;
   PetscFunctionBeginUser;

   // Compute left term
   left = 0;
   for (k=0; k<dim; k++) {
      // PetscPrintf(PETSC_COMM_WORLD, "SDD[%d]= %f\n",k,SDD[k]);
      left += user->weights[k]/PetscMax(SDD[k],user->epsilon);
   }
   left += user->weights[0]/PetscMax(SDD[0],user->epsilon); // last SDD is same as first SDD
   left =  PETSC_PI*PETSC_PI*PetscPowReal(left,-2.0);
   // calculate the right term
   right = user->epsilon;
   for (k=0; k<dim; k++) {
      if (right > SDD[k]){ // find minimum
         right = SDD[k];
      }
   }
   *DetD2u = left + right;
   PetscFunctionReturn(0);
}


/* This function computes `weights` variable in the MACtx object.
*/
PetscErrorCode ComputeWeights(PetscInt width, PetscInt order, MACtx *user) {
   PetscInt   i,M;
   PetscReal  a,h1,h2;
   PetscReal  *theta;
   PetscFunctionBeginUser;

   // Compute angles of L1 stencil points
   M = 2*width+1;
   PetscMalloc1(M,&theta);
   theta[width] = PETSC_PI/2;
   for (i=0; i<width; i++) {
      a = PetscAtan2Real(i,width-i);
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
   PetscFunctionReturn(0);
}


// Code to print projectinos
PetscErrorCode PrintProjection(DM da, MACtx *user) {
   // PetscErrorCode ierr;
   // DMDALocalInfo  info;
   // PetscInt i,j,k,count,width,N,M,Nk;
   // PetscInt *di, *dj;
   // PetscReal *xFwd, *xBak, *yFwd, *yBak;
   // PetscFunctionBeginUser;
   // // Get dimensions and other info from da
   // ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
   // N = info.my;     // total rows
   // M = info.mx;     // total cols
   // width = info.sw; // stencil width
   // // Compute stencil directions
   // Nk = 4*width+1;
   // PetscCalloc2(Nk,&di,Nk,&dj); // initilize 2 arrays of zeros
   // count = 1;
   // for (k=width; k>-width; k--) {
   //    di[count] = k;
   //    dj[count] = -PetscAbsReal(k)+width;
   //    count++;
   //    di[count] = -k;
   //    dj[count] = PetscAbsReal(k)-width;
   //    count++;
   // }
   // PetscPrintf(PETSC_COMM_WORLD,"Stencil directions for width = %d:\n",width);
   // for (k=0; k<Nk; k++) {
   //    PetscPrintf(PETSC_COMM_WORLD,"(%d) di,dj = (%d,%d)\n",k,di[k],dj[k]);
   // }
   // for (i=0; i<M; i++) {
   //    for (j=0; j<N; j++) {
   //       for (k=0; k<Nk; k++) {
   //          PetscPrintf(PETSC_COMM_WORLD,"i,j,k = (%d,%d,%d)\n",i,j,k);
   //       }
   //    }
   // }
   PetscFunctionReturn(0);
}


/* This function computes Si and Sj in the MACtx object.

   The derivative approximation at a node (i,j) on a square lattice
   depends on several neighboring nodes. Given some integer `width`, the
   set of all nodes (p,q) with L1 distance equal to `width` are a part of
   the stencil. There are always 4*width points. Among them, we only care
   about the directions that form an angle [0, pi) wrt. (i,j) because
   the rest are just mirrored across the center.

   We will call these points the forward stencil point in the kth direction
   for 0 < k < 2*width-1.
*/
PetscErrorCode ComputeFwdStencilDirs(PetscInt width, MACtx *user) {
   PetscInt   k,count,N;
   PetscFunctionBeginUser;

   // Compute angles of L1 stencil points
   N = 2*width;
   PetscCalloc2(N,&user->Si,N,&user->Sj);
   count = 0;
   for (k=width; -width<k; k--) {
      user->Si[count] = k;
      user->Sj[count] = -PetscAbsReal(k)+width;
      count++;
   }
   PetscFunctionReturn(0);
}


/* Compute di and dj
   This is Algorithm 5 in the Appendix.
   This function computes the projected version of Si and Sj when near the boundary.

   Si and Sj are the x and y relative stencil directions (positive or negative).
*/
PetscErrorCode ComputeProjectionIndeces(PetscReal *di, PetscReal *dj, PetscInt i, PetscInt j, PetscInt Si, PetscInt Sj, PetscInt Nx, PetscInt Ny) {
   PetscReal m;
   PetscInt check;
   PetscFunctionBeginUser;

   if (Si==0) {        // purely vertical
      *di = 0;
      *dj = (Sj>0)? Ny-j : -(1+j);
   } else if (Sj==0) { // purely horizontal
      *di = (Si>0)? Nx-i : -(1+i);
      *dj = 0;
   } else {            // mixed
      m = Sj/(PetscReal)Si; // slope
      if (Si>0 && Sj>0) {
         check = PetscAbsReal((Ny-j)/m) < PetscAbsReal(Nx-i);
         // check = Sj > Si;
         *di = (check)? (Ny-j)/m : Nx-i;
         *dj = (check)?     Ny-j : m*(Nx-i);
      } else if (Si>0 && Sj<0) {
         check = PetscAbsReal(Nx-i) < PetscAbsReal((1+j)/m);
         // check = -Sj < Si;
         *di = (check)?     Nx-i : -(1+j)/m;
         *dj = (check)? m*(Nx-i) : -(1+j);
      } else if (Si<0 && Sj<0) {
         check = PetscAbsReal((1+j)/m) < PetscAbsReal(1+i);
         // check = Si > Sj;
         *di = (check)? -(1+j)/m : -(1+i);
         *dj = (check)?   -(1+j) : -m*(1+i);
      } else if (Si<0 && Sj>0) {
         check = PetscAbsReal(1+i) < PetscAbsReal((Ny-j)/m);
         // check = -Si > Sj;
         *di = (check)?   -(1+i) : (Ny-j)/m;
         *dj = (check)? -m*(1+i) : Ny-j;
      } else {
         PetscPrintf(PETSC_COMM_WORLD," -- Unexpected error in projection");
      }
   }
   PetscFunctionReturn(0);
}


/* Computes SDD, hFwd, and hBak
   This gets used in the residue function and the Jacobian.
*/
PetscErrorCode ComputeSDD(DMDALocalInfo *info, PetscReal **au, MACtx *user, PetscInt i, PetscInt j, PetscReal x, PetscReal y, PetscReal *SDD, PetscReal *hFwd, PetscReal *hBak) {
   PetscInt   d,k,M,Ny,Nx,Si,Sj;
   PetscReal  hx,hy,di,dj,temp;
   PetscReal  *uFwd, *uBak; // u in the the forward and backward position for each direction k
   PetscFunctionBeginUser;

   // get info from DA
   Nx = info->mx; Ny = info->my;
   hx = (user->xmax - user->xmin)/(PetscReal)(Nx + 1);
   hy = (user->ymax - user->ymin)/(PetscReal)(Ny + 1);
   d  = info->sw;
   M  = d*2;
   PetscMalloc2(M,&uFwd,M,&uBak);
   if (d==1) {
      // width 1 case is hardcoded because it is simple
      if(j == info->my-1) {user->g_bdry(x,y+hy,0.0,user,&uFwd[1]);} else {uFwd[1] = au[j+1][i];}
      if(j == 0)          {user->g_bdry(x,y-hy,0.0,user,&uBak[1]);} else {uBak[1] = au[j-1][i];}
      if(i == info->mx-1) {user->g_bdry(x+hx,y,0.0,user,&uFwd[0]);} else {uFwd[0] = au[j][i+1];}
      if(i == 0)          {user->g_bdry(x-hx,y,0.0,user,&uBak[0]);} else {uBak[0] = au[j][i-1];}                                  // east
      // 2nd dir. deriv
      SDD[0] = (uFwd[0] - 2.0*au[j][i] + uBak[0])/(hx*hx); // horizontal centered-diff
      SDD[1] = (uFwd[1] - 2.0*au[j][i] + uBak[1])/(hy*hy); // vertical centered-diff
      // step sizes are static
      hFwd[0] = hx; hFwd[1] = hx;
      hBak[0] = hy; hBak[1] = hy;
   } else { // compute SDD for general width
      for (k=0; k<M; k++) { // loop over each stencil directions pairs
         Si = user->Si[k]; // stencil x displacement (positive)
         Sj = user->Sj[k]; // stencil y displacement (positive)
         // Forward point for direction k
         if (i+Si>=0 && i+Si<=Ny-1 && j+Sj<=Ny-1 && j+Sj>=0) {  // in-bound check
            // forward stencil point in the kth direction is in range
            uFwd[k] = au[j+Sj][i+Si];
            hFwd[k] = PetscSqrtReal(hx*hx*Si*Si + hy*hy*Sj*Sj);
         } else {
            // otherwise get the coordinates of the projection
            ComputeProjectionIndeces(&di,&dj,i,j,Si,Sj,Nx,Ny);
            user->g_bdry(x+hx*di,y+hy*dj,0.0,user,&temp);
            uFwd[k] = temp;
            hFwd[k] = PetscSqrtReal(hx*hx*di*di + hy*hy*dj*dj);
         }
         // Backward point for direction k
         if (i-Si>=0 && i-Si<=Ny-1 && j-Sj>=0 && j-Sj <= Ny-1) {
            // backward stencil point in the kth direction is in range
            uBak[k] = au[j-Sj][i-Si];
            hBak[k] = PetscSqrtReal(hx*hx*Si*Si + hy*hy*Sj*Sj);
         } else {
            // otherwise get the coordinates of the projection
            ComputeProjectionIndeces(&di,&dj,i,j,-Si,-Sj,Ny,Ny);
            user->g_bdry(x+hx*di,y+hy*dj,0.0,user,&temp);
            uBak[k] = temp;
            hBak[k] = PetscSqrtReal(hx*hx*di*di + hy*hy*dj*dj);
         }
      }
      // Use formula for generalized centered difference
      for (k=0; k<M; k++) {
         SDD[k] = 2.0*(hBak[k]*uFwd[k] - (hBak[k]+hFwd[k])*au[j][i] + hFwd[k]*uBak[k])/(hBak[k]*hFwd[k]*(hBak[k]+hFwd[k]));
      }
   }
   PetscFree2(uFwd,uBak);
   PetscFunctionReturn(0);
}