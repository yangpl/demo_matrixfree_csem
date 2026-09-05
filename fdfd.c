/* compute matrix vector product for applying Maxwell operator to a vector
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "solver.h"

#define id1(k,j,i) (i + n1*(j + (n2+1)*(k))) //index for Ex
#define id2(k,j,i) (i + (n1+1)*(j + n2*(k))) //index for Ey
#define id3(k,j,i) (i + (n1+1)*(j + (n2+1)*(k))) //index for Ez
#define idx(k,j,i) (i + (n1+1)*(j + n2*(k))) //index for Hx
#define idy(k,j,i) (i + n1*(j + (n2+1)*(k))) //index for Hy
#define idz(k,j,i) (i + n1*(j + n2*(k))) //index for Hz
#define idv(k,j,i) (i + (n1+1)*(j + (n2+1)*(k))) //index for phi
#define idHx(j,i) (i + (n1+1)*(j)) //index for Hxkm1 BC
#define idHy(j,i) (i + n1*(j)) //index for Hy Hykm1 BC

emf_t *emf;
complex *Hxkm1, *Hykm1;//ghost air layer Hx,Hz for boundary condition (BC)

void airbc_init(emf_t *emf);
void airbc_free(emf_t *emf);
void airbc_apply(emf_t *emf, complex *H3k, complex *H1km1, complex *H2km1);

void fdfd_init(emf_t *emf_)
{
  emf = emf_;//pointer for electromagnetic field (EMF)
  if(emf->preco==3){
    emf->lv = (emf->n1+1)*(emf->n2+1)*(emf->n3+1);//length of V
    emf->v = alloc1complex(emf->lv);
    emf->bv = alloc1complex(emf->lv);
    memset(emf->v, 0, emf->lv*sizeof(complex));
    memset(emf->bv, 0, emf->lv*sizeof(complex));
  }
  if(emf->airbc){
    if(emf->preco==1){
      emf->lEx = emf->n1*(emf->n2+1)*(emf->n3+1);//length of Ex
      emf->lEy = (emf->n1+1)*emf->n2*(emf->n3+1);//length of Ey
      emf->lEz = (emf->n1+1)*(emf->n2+1)*emf->n3;//length of Ez
      emf->lE = emf->lEx + emf->lEy + emf->lEz;//length of E=(Ex,Ey,Ez)^T
      emf->xx = alloc1complex(emf->lE);
      emf->yy = alloc1complex(emf->lE);
    }
    emf->n3 = emf->nz + emf->nb;
    airbc_init(emf);//constructor to initialize for airwave BC
    Hxkm1 = alloc1complex((emf->n1+1)*emf->n2);//ghost air layer Hx for BC
    Hykm1 = alloc1complex(emf->n1*(emf->n2+1));//ghost air layer Hy for BC
  }
  emf->lEx = emf->n1*(emf->n2+1)*(emf->n3+1);//length of Ex
  emf->lEy = (emf->n1+1)*emf->n2*(emf->n3+1);//length of Ey
  emf->lEz = (emf->n1+1)*(emf->n2+1)*emf->n3;//length of Ez
  emf->lE = emf->lEx + emf->lEy + emf->lEz;//length of E=(Ex,Ey,Ez)^T
  emf->lHx = (emf->n1+1)*emf->n2*emf->n3;//length of Hx
  emf->lHy = emf->n1*(emf->n2+1)*emf->n2;//length of Hy
  emf->lHz = emf->n1*emf->n2*(emf->n3+1);//length of Hz
  emf->lH = emf->lHx + emf->lHy + emf->lHz;//length of H=(Hx,Hy,Hz)^T
  emf->lv = (emf->n1+1)*(emf->n2+1)*(emf->n3+1);//length of V
}

void fdfd_free(emf_t *emf)
{
  if(emf->preco==3){
    free1complex(emf->v);
    free1complex(emf->bv);
  }
  if(emf->airbc){
    if(emf->preco==1){
      free1complex(emf->xx);
      free1complex(emf->yy);
    }
    airbc_free(emf);//deconstructor to free variables in airwave module
    free1complex(Hxkm1);//free ghost air layer Hx for BC
    free1complex(Hykm1);//free ghost air layer Hz for BC
  }
}


//y=Ax
//d1[i] = x_{i+0.5} - x_{i-0.5}, centered at i
//d1s[i] = x_{i+1} - x_i, half grid staggered, centered at i+0.5
void fdfd_apply(int n, complex *x, complex *y)
{
  int n1, n2, n3;
  int i, j, k, kstart;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  complex *Ex, *Ey, *Ez;
  complex *Hx, *Hy, *Hz;
  complex *bx, *by, *bz;
  double *sigma11, *sigma22, *sigma33;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  complex t1, t2, iwmu;

  iwmu = 1./emf->I_omega_mu0;
  kstart = (emf->airbc)?emf->nb:0;
  n1 = emf->n1;
  n2 = emf->n2;
  n3 = emf->n3;
  d1 = emf->d1;
  d2 = emf->d2;
  d3 = &emf->d3[kstart];
  d1s = emf->d1s;
  d2s = emf->d2s;
  d3s = &emf->d3s[kstart];
  sigma11 = &emf->sigma11[kstart][0][0];
  sigma22 = &emf->sigma22[kstart][0][0];
  sigma33 = &emf->sigma33[kstart][0][0];
  Ex = &x[0];
  Ey = &x[emf->lEx];
  Ez = &x[emf->lEx + emf->lEy];
  Hx = &emf->xh[0];
  Hy = &emf->xh[emf->lHx];
  Hz = &emf->xh[emf->lHx + emf->lHy];
  bx = &y[0];
  by = &y[emf->lEx];
  bz = &y[emf->lEx + emf->lEy];
  memset(y, 0, emf->lE*sizeof(complex));

  for(k=0; k<=n3; k++){
    kp1 = k+1;
    for(j=0; j<=n2; j++){
      jp1 = j+1;
      for(i=0; i<=n1; i++){
	ip1 = i+1;
	
	if(j<n2 && k<n3){
	  t1 = (Ez[id3(k,jp1,i)] - Ez[id3(k,j,i)])/d2s[j];
	  t2 = (Ey[id2(kp1,j,i)] - Ey[id2(k,j,i)])/d3s[k];
	  Hx[idx(k,j,i)] = (t1-t2)*iwmu;//\partial_y Ez - \partial_z Ey
	}
	if(i<n1 && k<n3){
	  t1 = (Ex[id1(kp1,j,i)] - Ex[id1(k,j,i)])/d3s[k];
	  t2 = (Ez[id3(k,j,ip1)] - Ez[id3(k,j,i)])/d1s[i];
	  Hy[idy(k,j,i)] = (t1-t2)*iwmu;//\partial_z Ex - \partial_x Ez
	}
	if(i<n1 && j<n2){
	  t1 = (Ey[id2(k,j,ip1)] - Ey[id2(k,j,i)])/d1s[i];
	  t2 = (Ex[id1(k,jp1,i)] - Ex[id1(k,j,i)])/d2s[j];
	  Hz[idz(k,j,i)] = (t1-t2)*iwmu;//\partial_x Ey - \partial_y Ex
	}
      }
    }
  }

  if(emf->airbc){
    airbc_apply(emf, Hz, Hxkm1, Hykm1);//Hz[z=0]->Hxkm1,Hykm1
    k = 0;
    for(j=0; j<n2; j++){
      jm1 = j-1;
      for(i=0; i<n1; i++){
	im1 = i-1;
	if(j>0){
	  t1 = (Hz[idz(k,j,i)]- Hz[idz(k,jm1,i)])/d2[j];
	  t2 = (Hy[idy(k,j,i)]- Hykm1[idHy(j,i)])/d3[k];
	  bx[id1(k,j,i)] = t1 - t2  - sigma11[id1(k,j,i)]*Ex[id1(k,j,i)];
	  if(emf->fvm) bx[id1(k,j,i)] *= d1s[i]*d2[j]*d3[k]*emf->I_omega_mu0;
	}
	if(i>0){
	  t1 = (Hx[idx(k,j,i)]-Hxkm1[idHx(j,i)])/d3[k];
	  t2 = (Hz[idz(k,j,i)]-Hz[idz(k,j,im1)])/d1[i];
	  by[id2(k,j,i)] = t1 - t2 - sigma22[id2(k,j,i)]*Ey[id2(k,j,i)];
	  if(emf->fvm) by[id2(k,j,i)] *= d1[i]*d2s[j]*d3[k]*emf->I_omega_mu0;
	}
      }//end for i
    }//end for j
  }//end if
  
 
  //\curl again
  for(k=0; k<n3; k++){
    km1 = k-1;
    for(j=0; j<n2; j++){
      jm1 = j-1;
      for(i=0; i<n1; i++){
	im1 = i-1;

	if(j>0 && k>0){
	  t1 = (Hz[idz(k,j,i)]- Hz[idz(k,jm1,i)])/d2[j];
	  t2 = (Hy[idy(k,j,i)]- Hy[idy(km1,j,i)])/d3[k];
	  bx[id1(k,j,i)] = t1 - t2  - sigma11[id1(k,j,i)]*Ex[id1(k,j,i)];
	  if(emf->fvm) bx[id1(k,j,i)] *= d1s[i]*d2[j]*d3[k]*emf->I_omega_mu0;
	}
	if(i>0 && k>0){
	  t1 = (Hx[idx(k,j,i)]-Hx[idx(km1,j,i)])/d3[k];
	  t2 = (Hz[idz(k,j,i)]-Hz[idz(k,j,im1)])/d1[i];
	  by[id2(k,j,i)] = t1 - t2 - sigma22[id2(k,j,i)]*Ey[id2(k,j,i)];
	  if(emf->fvm) by[id2(k,j,i)] *= d1[i]*d2s[j]*d3[k]*emf->I_omega_mu0;
	}
	if(i>0 && j>0){
	  t1 = (Hy[idy(k,j,i)]-Hy[idy(k,j,im1)])/d1[i];
	  t2 = (Hx[idx(k,j,i)]-Hx[idx(k,jm1,i)])/d2[j];
	  bz[id3(k,j,i)] = t1 - t2 - sigma33[id3(k,j,i)]*Ez[id3(k,j,i)];
	  if(emf->fvm) bz[id3(k,j,i)] *= d1[i]*d2[j]*d3s[k]*emf->I_omega_mu0;
	}
      }
    }
  }

}


//y=Ax, less efficient than fdfd_apply() due to repeated differencing of curl E at each index ijk
//d1[i] = x_{i+0.5} - x_{i-0.5}, centered at i
//d1s[i] = x_{i+1} - x_i, half grid staggered, centered at i+0.5
//this routine use Direchlet boundary condition (do not consider air-water boundary boundary condition)
void fdfd_apply_old(int n, complex *x, complex *y)
{
  int n1, n2, n3;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  complex *Ex, *Ey, *Ez;
  complex *bx, *by, *bz;
  complex t1, t2, t3, t4, iwmu;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double *sigma11, *sigma22, *sigma33;

  n1 = emf->n1;
  n2 = emf->n2;
  n3 = emf->n3;
  d1 = emf->d1;
  d2 = emf->d2;
  d3 = emf->d3;
  d1s = emf->d1s;
  d2s = emf->d2s;
  d3s = emf->d3s;
  sigma11 = &emf->sigma11[0][0][0];
  sigma22 = &emf->sigma22[0][0][0];
  sigma33 = &emf->sigma33[0][0][0];
  Ex = &x[0];
  Ey = &x[emf->lEx];
  Ez = &x[emf->lEx + emf->lEy];
  bx = &y[0];
  by = &y[emf->lEx];
  bz = &y[emf->lEx + emf->lEy];

  iwmu = 1./emf->I_omega_mu0;
  memset(y, 0, emf->lE*sizeof(complex));
  for(k=0; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=0; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=0; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	if(j>0 && k>0){
	  t1 = (Ey[id2(k,j,ip1)] - Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)] - Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t2 = (Ey[id2(k,jm1,ip1)] - Ey[id2(k,jm1,i)])/d1s[i] - (Ex[id1(k,j,i)] - Ex[id1(k,jm1,i)])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	  t3 = (Ex[id1(kp1,j,i)] - Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)] - Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t4 = (Ex[id1(k,j,i)] - Ex[id1(km1,j,i)])/d3s[km1] - (Ez[id3(km1,j,ip1)] - Ez[id3(km1,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	  bx[id1(k,j,i)] = ((t1-t2)/d2[j] - (t3-t4)/d3[k])*iwmu - sigma11[id1(k,j,i)]*Ex[id1(k,j,i)];
	  if(emf->fvm) bx[id1(k,j,i)] *= d1s[i]*d2[j]*d3[k]*emf->I_omega_mu0;
	}
	if(i>0 && k>0){
	  t1 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t2 = (Ez[id3(km1,jp1,i)]-Ez[id3(km1,j,i)])/d2s[j] - (Ey[id2(k,j,i)]-Ey[id2(km1,j,i)])/d3s[km1];//\partial_y Ez - \partial_z Ey
	  t3 = (Ey[id2(k,j,ip1)]-Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)]-Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t4 = (Ey[id2(k,j,i)]-Ey[id2(k,j,im1)])/d1s[im1] - (Ex[id1(k,jp1,im1)]-Ex[id1(k,j,im1)])/d2s[j];//\partial_x Ey - \partial_y Ex
	  by[id2(k,j,i)] = ((t1-t2)/d3[k] - (t3-t4)/d1[i])*iwmu - sigma22[id2(k,j,i)]*Ey[id2(k,j,i)];
	  if(emf->fvm) by[id2(k,j,i)] *= d1[i]*d2s[j]*d3[k]*emf->I_omega_mu0;
	}
	if(i>0 && j>0){
	  t1 = (Ex[id1(kp1,j,i)]-Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)]-Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t2 = (Ex[id1(kp1,j,im1)]-Ex[id1(k,j,im1)])/d3s[k] - (Ez[id3(k,j,i)]-Ez[id3(k,j,im1)])/d1s[im1];//\partial_z Ex - \partial_x Ez
	  t3 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t4 = (Ez[id3(k,j,i)]-Ez[id3(k,jm1,i)])/d2s[jm1] - (Ey[id2(kp1,jm1,i)]-Ey[id2(k,jm1,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	  bz[id3(k,j,i)] = ((t1-t2)/d1[i] - (t3-t4)/d2[j])*iwmu - sigma33[id3(k,j,i)]*Ez[id3(k,j,i)];
	  if(emf->fvm) bz[id3(k,j,i)] *= d1[i]*d2[j]*d3s[k]*emf->I_omega_mu0;
	}
      }//end for i
    }//end for j
  }//end for k

}

//direct LU-SGS for double curl equation as the preconditioner
void lusgs_apply(int n, complex *x, complex *y)
{
  int n1, n2, n3;
  int lEx, lEy, lEz, lE;
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  double vol;
  complex *Ex, *Ey, *Ez;
  complex *bx, *by, *bz;
  complex t1, t2, t3, t4, aij_xj, aii, iwmu;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double *sigma11, *sigma22, *sigma33;

  iwmu = 1./emf->I_omega_mu0;
  n1 = emf->n1;
  n2 = emf->n2;
  n3 = emf->n3;
  d1 = emf->d1;
  d2 = emf->d2;
  d3 = emf->d3;
  d1s = emf->d1s;
  d2s = emf->d2s;
  d3s = emf->d3s;
  sigma11 = &emf->sigma11[0][0][0];
  sigma22 = &emf->sigma22[0][0][0];
  sigma33 = &emf->sigma33[0][0][0];
  lEx = n1*(n2+1)*(n3+1);
  lEy = (n1+1)*n2*(n3+1);
  lEz = (n1+1)*(n2+1)*n3;
  lE = lEx + lEy + lEz;
  if(emf->airbc){
    n3 += emf->nb;//extends with air layers
    lEx = n1*(n2+1)*(n3+1);
    lEy = (n1+1)*n2*(n3+1);
    lEz = (n1+1)*(n2+1)*n3;
    lE = lEx + lEy + lEz;
    Ex = &emf->yy[0];
    Ey = &emf->yy[lEx];
    Ez = &emf->yy[lEx + lEy];
    bx = &emf->xx[0];
    by = &emf->xx[lEx];
    bz = &emf->xx[lEx + lEy];
    memset(emf->xx, 0, lE*sizeof(complex));
    memset(emf->yy, 0, lE*sizeof(complex));
    memcpy(&emf->xx[emf->n1*(emf->n2+1)*emf->nb], x, emf->lEx*sizeof(complex));
    memcpy(&emf->xx[(emf->n1+1)*emf->n2*emf->nb + lEx], &x[emf->lEx], emf->lEy*sizeof(complex));
    memcpy(&emf->xx[(emf->n1+1)*(emf->n2+1)*emf->nb + lEx + lEy], &x[emf->lEx + emf->lEy], emf->lEz*sizeof(complex));
  }else{
    Ex = &y[0];
    Ey = &y[emf->lEx];
    Ez = &y[emf->lEx + emf->lEy];
    bx = &x[0];
    by = &x[emf->lEx];
    bz = &x[emf->lEx + emf->lEy];
    memset(y, 0, emf->lE*sizeof(complex));
  }
  
  //forward Gauss-Seidel
  for(k=1; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=0; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ex[id1(k,j,i)] = 0;	  
	t1 = (Ey[id2(k,j,ip1)] - Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)] - Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[id2(k,jm1,ip1)] - Ey[id2(k,jm1,i)])/d1s[i] - (Ex[id1(k,j,i)] - Ex[id1(k,jm1,i)])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[id1(kp1,j,i)] - Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)] - Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[id1(k,j,i)] - Ex[id1(km1,j,i)])/d3s[km1] - (Ez[id3(km1,j,ip1)] - Ez[id3(km1,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	aij_xj = ((t1-t2)/d2[j] - (t3-t4)/d3[k])*iwmu - sigma11[id1(k,j,i)]*Ex[id1(k,j,i)];
	aii = ((1./d2s[j] + 1./d2s[jm1])/d2[j] + (1./d3s[k] + 1./d3s[km1])/d3[k])*iwmu - sigma11[id1(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ex[id1(k,j,i)] = (bx[id1(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=1; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=0; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ey[id2(k,j,i)] = 0;
	t1 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[id3(km1,jp1,i)]-Ez[id3(km1,j,i)])/d2s[j] - (Ey[id2(k,j,i)]-Ey[id2(km1,j,i)])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[id2(k,j,ip1)]-Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)]-Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[id2(k,j,i)]-Ey[id2(k,j,im1)])/d1s[im1] - (Ex[id1(k,jp1,im1)]-Ex[id1(k,j,im1)])/d2s[j];//\partial_x Ey - \partial_y Ex
	aij_xj = ((t1-t2)/d3[k] - (t3-t4)/d1[i])*iwmu - sigma22[id2(k,j,i)]*Ey[id2(k,j,i)];
	aii = ((1./d3s[k] + 1./d3s[km1])/d3[k] + (1./d1s[i] + 1./d1s[im1])/d1[i])*iwmu - sigma22[id2(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ey[id2(k,j,i)] = (by[id2(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=0; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ez[id3(k,j,i)] = 0;
	t1 = (Ex[id1(kp1,j,i)]-Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)]-Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[id1(kp1,j,im1)]-Ex[id1(k,j,im1)])/d3s[k] - (Ez[id3(k,j,i)]-Ez[id3(k,j,im1)])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[id3(k,j,i)]-Ez[id3(k,jm1,i)])/d2s[jm1] - (Ey[id2(kp1,jm1,i)]-Ey[id2(k,jm1,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	aij_xj = ((t1-t2)/d1[i] - (t3-t4)/d2[j])*iwmu - sigma33[id3(k,j,i)]*Ez[id3(k,j,i)];
	aii = ((1./d1s[i] + 1./d1s[im1])/d1[i] + (1./d2s[j] + 1./d2s[jm1])/d2[j])*iwmu - sigma33[id3(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ez[id3(k,j,i)] = (bz[id3(k,j,i)] - aij_xj)/aii;
      }//end for i
    }//end for j
  }//end for k
  //backward Gauss-Seidel
  for(k=n3-1; k>=0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ez[id3(k,j,i)] = 0;
	t1 = (Ex[id1(kp1,j,i)]-Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)]-Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[id1(kp1,j,im1)]-Ex[id1(k,j,im1)])/d3s[k] - (Ez[id3(k,j,i)]-Ez[id3(k,j,im1)])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[id3(k,j,i)]-Ez[id3(k,jm1,i)])/d2s[jm1] - (Ey[id2(kp1,jm1,i)]-Ey[id2(k,jm1,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	aij_xj = ((t1-t2)/d1[i] - (t3-t4)/d2[j])*iwmu - sigma33[id3(k,j,i)]*Ez[id3(k,j,i)];
	aii = ((1./d1s[i] + 1./d1s[im1])/d1[i] + (1./d2s[j] + 1./d2s[jm1])/d2[j])*iwmu - sigma33[id3(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ez[id3(k,j,i)] = (bz[id3(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=n3-1; k>0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>=0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ey[id2(k,j,i)] = 0;
	t1 = (Ez[id3(k,jp1,i)]-Ez[id3(k,j,i)])/d2s[j] - (Ey[id2(kp1,j,i)]-Ey[id2(k,j,i)])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[id3(km1,jp1,i)]-Ez[id3(km1,j,i)])/d2s[j] - (Ey[id2(k,j,i)]-Ey[id2(km1,j,i)])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[id2(k,j,ip1)]-Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)]-Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[id2(k,j,i)]-Ey[id2(k,j,im1)])/d1s[im1] - (Ex[id1(k,jp1,im1)]-Ex[id1(k,j,im1)])/d2s[j];//\partial_x Ey - \partial_y Ex
	aij_xj = ((t1-t2)/d3[k] - (t3-t4)/d1[i])*iwmu - sigma22[id2(k,j,i)]*Ey[id2(k,j,i)];
	aii = ((1./d3s[k] + 1./d3s[km1])/d3[k] + (1./d1s[i] + 1./d1s[im1])/d1[i])*iwmu - sigma22[id2(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ey[id2(k,j,i)] = (by[id2(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=n3-1; k>0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>=0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ex[id1(k,j,i)] = 0;	  
	t1 = (Ey[id2(k,j,ip1)] - Ey[id2(k,j,i)])/d1s[i] - (Ex[id1(k,jp1,i)] - Ex[id1(k,j,i)])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[id2(k,jm1,ip1)] - Ey[id2(k,jm1,i)])/d1s[i] - (Ex[id1(k,j,i)] - Ex[id1(k,jm1,i)])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[id1(kp1,j,i)] - Ex[id1(k,j,i)])/d3s[k] - (Ez[id3(k,j,ip1)] - Ez[id3(k,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[id1(k,j,i)] - Ex[id1(km1,j,i)])/d3s[km1] - (Ez[id3(km1,j,ip1)] - Ez[id3(km1,j,i)])/d1s[i];//\partial_z Ex - \partial_x Ez
	aij_xj = ((t1-t2)/d2[j] - (t3-t4)/d3[k])*iwmu - sigma11[id1(k,j,i)]*Ex[id1(k,j,i)];
	aii = ((1./d2s[j] + 1./d2s[jm1])/d2[j] + (1./d3s[k] + 1./d3s[km1])/d3[k])*iwmu - sigma11[id1(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	Ex[id1(k,j,i)] = (bx[id1(k,j,i)] - aij_xj)/aii;
      }//end for i
    }//end for j
  }//end for k

  if(emf->airbc){
    memcpy(y, &emf->yy[emf->n1*(emf->n2+1)*emf->nb], emf->lEx*sizeof(complex));
    memcpy(&y[emf->lEx], &emf->yy[(emf->n1+1)*emf->n2*emf->nb + lEx], emf->lEy*sizeof(complex));
    memcpy(&y[emf->lEx + emf->lEy], &emf->yy[(emf->n1+1)*(emf->n2+1)*emf->nb + lEx + lEy], emf->lEz*sizeof(complex));
  }
}

//A-phi preconditioning
void aphi_apply(int n, complex *x, complex *y)
{
  int n1, n2, n3;
  int i, j, k, kstart;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  complex t1, t2, t3, t4, t5, t6, aij_xj, aii, iwmu;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double *sigma11, *sigma22, *sigma33;
  complex *ax, *ay, *az;
  complex *bx, *by, *bz;
  complex *Ex, *Ey, *Ez;
  complex *v, *bv;
  double vol;

  kstart = (emf->airbc)?emf->nb:0;
  n1 = emf->n1;
  n2 = emf->n2;
  n3 = emf->n3;
  d1 = emf->d1;
  d2 = emf->d2;
  d3 = &emf->d3[kstart];
  d1s = emf->d1s;
  d2s = emf->d2s;
  d3s = &emf->d3s[kstart];
  sigma11 = &emf->sigma11[kstart][0][0];
  sigma22 = &emf->sigma22[kstart][0][0];
  sigma33 = &emf->sigma33[kstart][0][0];
  v = emf->v;
  bv = emf->bv;
  iwmu = 1./emf->I_omega_mu0;
  
  //1. b= (I GradT)^T*x: [bx,by,bz]=x; bv = -div(x)
  ax = &y[0];
  ay = &y[emf->lEx];
  az = &y[emf->lEx + emf->lEy];
  bx = &x[0];
  by = &x[emf->lEx];
  bz = &x[emf->lEx + emf->lEy];
  memset(bv, 0, emf->lv*sizeof(complex));
  for(k=1; k<n3; k++){
    km1 = k-1;
    for(j=1; j<n2; j++){
      jm1 = j-1;
      for(i=1; i<n1; i++){
	im1 = i-1;
	bv[idv(k,j,i)] = 0;
	bv[idv(k,j,i)] += -(bx[id1(k,j,i)] - bx[id1(k,j,im1)])/d1[i];
	bv[idv(k,j,i)] += -(by[id2(k,j,i)] - by[id2(k,jm1,i)])/d2[j];
	bv[idv(k,j,i)] += -(bz[id3(k,j,i)] - bz[id3(km1,j,i)])/d3[k];
      }
    }
  }
  
  //2. solve Aav*xx = (bx, by, bz, bv)^T by LU-SGS, xx=(ax, ay, az, v)^T
  // LU-SGS is equivalent to apply a forward and backward GS using 0 as initial input vector
  memset(ax, 0, emf->lEx*sizeof(complex));
  memset(ay, 0, emf->lEy*sizeof(complex));
  memset(az, 0, emf->lEz*sizeof(complex));
  memset(v, 0, emf->lv*sizeof(complex));
  //note we never update potential fields a and v at the boundary (assume they=0)
  //forward Gauss-Seidel
  for(k=1; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=0; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	ax[id1(k,j,i)] = 0;	  
	t1 = (ay[id2(k,j,ip1)] - ay[id2(k,j,i)])/d1s[i] - (ax[id1(k,jp1,i)] - ax[id1(k,j,i)])/d2s[j];//\partial_x ay - \partial_y ax
	t2 = (ay[id2(k,jm1,ip1)] - ay[id2(k,jm1,i)])/d1s[i] - (ax[id1(k,j,i)] - ax[id1(k,jm1,i)])/d2s[jm1];//\partial_x ay - \partial_y ax
	t3 = (ax[id1(kp1,j,i)] - ax[id1(k,j,i)])/d3s[k] - (az[id3(k,j,ip1)] - az[id3(k,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	t4 = (ax[id1(k,j,i)] - ax[id1(km1,j,i)])/d3s[km1] - (az[id3(km1,j,ip1)] - az[id3(km1,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	aij_xj = ((t1-t2)/d2[j] - (t3-t4)/d3[k])*iwmu - sigma11[id1(k,j,i)]*(ax[id1(k,j,i)] + (v[idv(k,j,ip1)]-v[idv(k,j,i)])/d1s[i]);
	aii = ((1./d2s[j] + 1./d2s[jm1])/d2[j] + (1./d3s[k] + 1./d3s[km1])/d3[k])*iwmu - sigma11[id1(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	ax[id1(k,j,i)] = (bx[id1(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=1; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=0; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	ay[id2(k,j,i)] = 0;
	t1 = (az[id3(k,jp1,i)]-az[id3(k,j,i)])/d2s[j] - (ay[id2(kp1,j,i)]-ay[id2(k,j,i)])/d3s[k];//\partial_y az - \partial_z ay
	t2 = (az[id3(km1,jp1,i)]-az[id3(km1,j,i)])/d2s[j] - (ay[id2(k,j,i)]-ay[id2(km1,j,i)])/d3s[km1];//\partial_y az - \partial_z ay
	t3 = (ay[id2(k,j,ip1)]-ay[id2(k,j,i)])/d1s[i] - (ax[id1(k,jp1,i)]-ax[id1(k,j,i)])/d2s[j];//\partial_x ay - \partial_y ax
	t4 = (ay[id2(k,j,i)]-ay[id2(k,j,im1)])/d1s[im1] - (ax[id1(k,jp1,im1)]-ax[id1(k,j,im1)])/d2s[j];//\partial_x ay - \partial_y ax
	aij_xj = ((t1-t2)/d3[k] - (t3-t4)/d1[i])*iwmu - sigma22[id2(k,j,i)]*(ay[id2(k,j,i)] + (v[idv(k,jp1,i)]-v[idv(k,j,i)])/d2s[j]);
	aii = ((1./d3s[k] + 1./d3s[km1])/d3[k] + (1./d1s[i] + 1./d1s[im1])/d1[i])*iwmu - sigma22[id2(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	ay[id2(k,j,i)] = (by[id2(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=0; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	az[id3(k,j,i)] = 0;
	t1 = (ax[id1(kp1,j,i)]-ax[id1(k,j,i)])/d3s[k] - (az[id3(k,j,ip1)]-az[id3(k,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	t2 = (ax[id1(kp1,j,im1)]-ax[id1(k,j,im1)])/d3s[k] - (az[id3(k,j,i)]-az[id3(k,j,im1)])/d1s[im1];//\partial_z ax - \partial_x az
	t3 = (az[id3(k,jp1,i)]-az[id3(k,j,i)])/d2s[j] - (ay[id2(kp1,j,i)]-ay[id2(k,j,i)])/d3s[k];//\partial_y az - \partial_z ay
	t4 = (az[id3(k,j,i)]-az[id3(k,jm1,i)])/d2s[jm1] - (ay[id2(kp1,jm1,i)]-ay[id2(k,jm1,i)])/d3s[k];//\partial_y az - \partial_z ay
	aij_xj = ((t1-t2)/d1[i] - (t3-t4)/d2[j])*iwmu - sigma33[id3(k,j,i)]*(az[id3(k,j,i)] + (v[idv(kp1,j,i)] - v[idv(k,j,i)])/d3s[k]);
	aii = ((1./d1s[i] + 1./d1s[im1])/d1[i] + (1./d2s[j] + 1./d2s[jm1])/d2[j])*iwmu - sigma33[id3(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	az[id3(k,j,i)] = (bz[id3(k,j,i)] - aij_xj)/aii;
      }//end for i
    }//end for j
  }//end for k
  for(k=1; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=1; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=1; i<n1; i++){
        ip1 = i+1;
        im1 = i-1;

        v[idv(k,j,i)] = 0;
        t1 = sigma11[id1(k,j,i)] * (ax[id1(k,j,i)] + (v[idv(k,j,ip1)] - v[idv(k,j,i)])/d1s[i]);
        t2 = emf->sigma11[k][j][im1] * (ax[id1(k,j,im1)] + (v[idv(k,j,i)] - v[idv(k,j,im1)])/d1s[im1]);
        t3 = sigma22[id2(k,j,i)] * (ay[id2(k,j,i)] + (v[idv(k,jp1,i)] - v[idv(k,j,i)])/d2s[j]);
        t4 = emf->sigma22[k][jm1][i] * (ay[id2(k,jm1,i)] + (v[idv(k,j,i)] - v[idv(k,jm1,i)])/d2s[jm1]);
        t5 = sigma33[id3(k,j,i)] * (az[id3(k,j,i)] + (v[idv(kp1,j,i)] - v[idv(k,j,i)])/d3s[k]);
        t6 = emf->sigma33[km1][j][i] * (az[id3(km1,j,i)] + (v[idv(k,j,i)] - v[idv(km1,j,i)])/d3s[km1]);
        aij_xj = ((t1-t2)/d1[i] + (t3-t4)/d2[j] + (t5-t6)/d3[k]);
        aii = -((1./d1s[i]+1./d1s[im1])/d1[i] + (1./d2s[j]+1./d2s[jm1])/d2[j] + (1./d3s[k]+1./d3s[km1])/d3[k]);
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
        v[idv(k,j,i)] = (bv[idv(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  //backward Gauss-Seidel
  for(k=n3-1; k>=1; k--){
    kp1 = k+1;
    km1 = k-1;
    for(j=n2-1; j>=1; j--){
      jp1 = j+1;
      jm1 = j-1;
      for(i=n1-1; i>=1; i--){
        ip1 = i+1;
        im1 = i-1;

	v[idv(k,j,i)] = 0;
        t1 = sigma11[id1(k,j,i)] * (ax[id1(k,j,i)] + (v[idv(k,j,ip1)] - v[idv(k,j,i)])/d1s[i]);
        t2 = emf->sigma11[k][j][im1] * (ax[id1(k,j,im1)] + (v[idv(k,j,i)] - v[idv(k,j,im1)])/d1s[im1]);
        t3 = sigma22[id2(k,j,i)] * (ay[id2(k,j,i)] + (v[idv(k,jp1,i)] - v[idv(k,j,i)])/d2s[j]);
        t4 = emf->sigma22[k][jm1][i] * (ay[id2(k,jm1,i)] + (v[idv(k,j,i)] - v[idv(k,jm1,i)])/d2s[jm1]);
        t5 = sigma33[id3(k,j,i)] * (az[id3(k,j,i)] + (v[idv(kp1,j,i)] - v[idv(k,j,i)])/d3s[k]);
        t6 = emf->sigma33[km1][j][i] * (az[id3(km1,j,i)] + (v[idv(k,j,i)] - v[idv(km1,j,i)])/d3s[km1]);
        aij_xj = ((t1-t2)/d1[i] + (t3-t4)/d2[j] + (t5-t6)/d3[k]);
        aii = -((1./d1s[i]+1./d1s[im1])/d1[i] + (1./d2s[j]+1./d2s[jm1])/d2[j] + (1./d3s[k]+1./d3s[km1])/d3[k]);
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
        v[idv(k,j,i)] = (bv[idv(k,j,i)] - aij_xj)/aii;
      }
    }
  }
 
  for(k=n3-1; k>=0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	az[id3(k,j,i)] = 0;
	t1 = (ax[id1(kp1,j,i)]-ax[id1(k,j,i)])/d3s[k] - (az[id3(k,j,ip1)]-az[id3(k,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	t2 = (ax[id1(kp1,j,im1)]-ax[id1(k,j,im1)])/d3s[k] - (az[id3(k,j,i)]-az[id3(k,j,im1)])/d1s[im1];//\partial_z ax - \partial_x az
	t3 = (az[id3(k,jp1,i)]-az[id3(k,j,i)])/d2s[j] - (ay[id2(kp1,j,i)]-ay[id2(k,j,i)])/d3s[k];//\partial_y az - \partial_z ay
	t4 = (az[id3(k,j,i)]-az[id3(k,jm1,i)])/d2s[jm1] - (ay[id2(kp1,jm1,i)]-ay[id2(k,jm1,i)])/d3s[k];//\partial_y az - \partial_z ay
	aij_xj = ((t1-t2)/d1[i] - (t3-t4)/d2[j])*iwmu - sigma33[id3(k,j,i)]*(az[id3(k,j,i)] + (v[idv(kp1,j,i)] - v[idv(k,j,i)])/d3s[k]);
	aii = ((1./d1s[i] + 1./d1s[im1])/d1[i] + (1./d2s[j] + 1./d2s[jm1])/d2[j])*iwmu - sigma33[id3(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	az[id3(k,j,i)] = (bz[id3(k,j,i)] - aij_xj)/aii;
      }//end for i
    }//end for j
  }//end for k
 
  for(k=n3-1; k>0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>=0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	ay[id2(k,j,i)] = 0;
	t1 = (az[id3(k,jp1,i)]-az[id3(k,j,i)])/d2s[j] - (ay[id2(kp1,j,i)]-ay[id2(k,j,i)])/d3s[k];//\partial_y az - \partial_z ay
	t2 = (az[id3(km1,jp1,i)]-az[id3(km1,j,i)])/d2s[j] - (ay[id2(k,j,i)]-ay[id2(km1,j,i)])/d3s[km1];//\partial_y az - \partial_z ay
	t3 = (ay[id2(k,j,ip1)]-ay[id2(k,j,i)])/d1s[i] - (ax[id1(k,jp1,i)]-ax[id1(k,j,i)])/d2s[j];//\partial_x ay - \partial_y ax
	t4 = (ay[id2(k,j,i)]-ay[id2(k,j,im1)])/d1s[im1] - (ax[id1(k,jp1,im1)]-ax[id1(k,j,im1)])/d2s[j];//\partial_x ay - \partial_y ax
	aij_xj = ((t1-t2)/d3[k] - (t3-t4)/d1[i])*iwmu - sigma22[id2(k,j,i)]*(ay[id2(k,j,i)] + (v[idv(k,jp1,i)]-v[idv(k,j,i)])/d2s[j]);
	aii = ((1./d3s[k] + 1./d3s[km1])/d3[k] + (1./d1s[i] + 1./d1s[im1])/d1[i])*iwmu - sigma22[id2(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	ay[id2(k,j,i)] = (by[id2(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  for(k=n3-1; k>0; k--){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=n2-1; j>0; j--){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=n1-1; i>=0; i--){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	ax[id1(k,j,i)] = 0;	  
	t1 = (ay[id2(k,j,ip1)] - ay[id2(k,j,i)])/d1s[i] - (ax[id1(k,jp1,i)] - ax[id1(k,j,i)])/d2s[j];//\partial_x ay - \partial_y ax
	t2 = (ay[id2(k,jm1,ip1)] - ay[id2(k,jm1,i)])/d1s[i] - (ax[id1(k,j,i)] - ax[id1(k,jm1,i)])/d2s[jm1];//\partial_x ay - \partial_y ax
	t3 = (ax[id1(kp1,j,i)] - ax[id1(k,j,i)])/d3s[k] - (az[id3(k,j,ip1)] - az[id3(k,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	t4 = (ax[id1(k,j,i)] - ax[id1(km1,j,i)])/d3s[km1] - (az[id3(km1,j,ip1)] - az[id3(km1,j,i)])/d1s[i];//\partial_z ax - \partial_x az
	aij_xj = ((t1-t2)/d2[j] - (t3-t4)/d3[k])*iwmu - sigma11[id1(k,j,i)]*(ax[id1(k,j,i)] + (v[idv(k,j,ip1)]-v[idv(k,j,i)])/d1s[i]);
	aii = ((1./d2s[j] + 1./d2s[jm1])/d2[j] + (1./d3s[k] + 1./d3s[km1])/d3[k])*iwmu - sigma11[id1(k,j,i)];
	if(emf->fvm){
	  vol = d1s[i]*d2[j]*d3[k];
	  aij_xj *= vol*emf->I_omega_mu0;
	  aii *= vol*emf->I_omega_mu0;
	}
	ax[id1(k,j,i)] = (bx[id1(k,j,i)] - aij_xj)/aii;
      }
    }
  }
  
  //3. y=(I Grad)*xx; E = A + Grad(v)
  Ex = &y[0];
  Ey = &y[emf->lEx];
  Ez = &y[emf->lEx + emf->lEy];
  for(k=1; k<n3; k++){
    for(j=1; j<n2; j++){
      for(i=0; i<n1; i++){
	ip1 = i+1;
	Ex[id1(k,j,i)] = ax[id1(k,j,i)] + (v[idv(k,j,ip1)]-v[idv(k,j,i)])/d1s[i];//Ex = ax + \partial_x v
      }
    }
  }
  for(k=1; k<n3; k++){
    for(j=0; j<n2; j++){
      jp1 = j+1;
      for(i=1; i<n1; i++){
	Ey[id2(k,j,i)] = ay[id2(k,j,i)] + (v[idv(k,jp1,i)]-v[idv(k,j,i)])/d2s[j];//Ey = ay + \partial_y v
      }
    }
  }
  for(k=0; k<n3; k++){
    kp1 = k+1;
    for(j=1; j<n2; j++){
      for(i=1; i<n1; i++){
	Ez[id3(k,j,i)] = az[id3(k,j,i)] + (v[idv(kp1,j,i)]-v[idv(k,j,i)])/d3s[k];//Ez = az + \partial_z v
      }
    }
  }

}

#undef id1
#undef id2
#undef id3
#undef idx
#undef idy
#undef idz
#undef idv


