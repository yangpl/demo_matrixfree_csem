/* inject source, extract EM data at receiver locations
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "emf.h"

#define id1(k,j,i) (i + emf->n1*(j + (emf->n2+1)*(k))) 
#define id2(k,j,i) (i + (emf->n1+1)*(j + emf->n2*(k)))
#define id3(k,j,i) (i + (emf->n1+1)*(j + (emf->n2+1)*(k)))
#define idx(k,j,i) (i + (emf->n1+1)*(j + emf->n2*(k)))
#define idy(k,j,i) (i + emf->n1*(j + (emf->n2+1)*(k)))
#define idz(k,j,i) (i + emf->n1*(j + emf->n2*(k)))

/*< find the index k in x[] such that x[k]<= val <x[k+1] >*/
int find_index(int n, double *x, double val);

/*< inject source >*/
void inject_source(acq_t *acq, emf_t *emf, complex *b)
{
  int ic, isrc, i, j, k;
  double w1, w2, w3, vol;
  complex src_vol;
  complex *bx, *by, *bz;

  int n1 = emf->nx + 2*emf->nb;
  int n2 = emf->ny + 2*emf->nb;
  int n3 = emf->nz + 2*emf->nb;
  bx = &b[0];
  by = &b[emf->lEx];
  bz = &b[emf->lEx + emf->lEy];

  for(ic=0; ic<emf->nchsrc; ic++){
    if(strcmp(emf->chsrc[ic],"Ex")==0){
      //Jx(i+0.5,j,k)
      for(isrc=0; isrc<acq->nsrc; isrc++){
	i = find_index(n1, emf->x1s, acq->src_x1[isrc]);
	j = find_index(n2+1, emf->x2, acq->src_x2[isrc]);
	k = find_index(n3+1, emf->x3, acq->src_x3[isrc]);
	
	w1 = (acq->src_x1[isrc] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
	w2 = (acq->src_x2[isrc] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
	w3 = (acq->src_x3[isrc] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
	vol = (emf->x1s[i+1]- emf->x1s[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3[k+1]- emf->x3[k]);
	src_vol = 1./vol;
	if(emf->fvm) src_vol = emf->I_omega_mu0;
	src_vol /= emf->rho_water;

	//Ex is continuous, convert Jx to Ex, interpolate and back to Ex at local gridpoint
	if(emf->airbc) k -= emf->nb;
	bx[id1(k,j,i)] += src_vol/emf->sigma11[k][j][i]*(1.-w1)*(1.-w2)*(1.-w3);
	bx[id1(k,j,i+1)] += src_vol/emf->sigma11[k][j][i+1]*w1*(1.-w2)*(1.-w3);
	bx[id1(k,j+1,i)] += src_vol/emf->sigma11[k][j+1][i]*(1.-w1)*w2*(1.-w3);
	bx[id1(k,j+1,i+1)] += src_vol/emf->sigma11[k][j+1][i+1]*w1*w2*(1.-w3);
	bx[id1(k+1,j,i)] += src_vol/emf->sigma11[k+1][j][i]*(1.-w1)*(1.-w2)*w3;
	bx[id1(k+1,j,i+1)] += src_vol/emf->sigma11[k+1][j][i+1]*w1*(1.-w2)*w3;
	bx[id1(k+1,j+1,i)] += src_vol/emf->sigma11[k+1][j+1][i]*(1.-w1)*w2*w3;
	bx[id1(k+1,j+1,i+1)] += src_vol/emf->sigma11[k+1][j+1][i+1]*w1*w2*w3;
      }
    }else if(strcmp(emf->chsrc[ic],"Ey")==0){
      //Jy(i,j+0.5,k)
      for(isrc=0; isrc<acq->nsrc; isrc++){
	i = find_index(n1+1, emf->x1, acq->src_x1[isrc]);
	j = find_index(n2, emf->x2s, acq->src_x2[isrc]);
	k = find_index(n3+1, emf->x3, acq->src_x3[isrc]);
	
	w1 = (acq->src_x1[isrc] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
	w2 = (acq->src_x2[isrc] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
	w3 = (acq->src_x3[isrc] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
	vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2s[j+1]- emf->x2s[j])*(emf->x3[k+1]- emf->x3[k]);
	src_vol = 1./vol;
	if(emf->fvm) src_vol = emf->I_omega_mu0;
	src_vol /= emf->rho_water;

	//Ey is continuous, convert Jx to Ey, interpolate and back to Ey at local gridpoint
	if(emf->airbc) k -= emf->nb;
	by[id2(k,j,i)] += src_vol/emf->sigma22[k][j][i]*(1.-w1)*(1.-w2)*(1.-w3);
	by[id2(k,j,i+1)] += src_vol/emf->sigma22[k][j][i+1]*w1*(1.-w2)*(1.-w3);
	by[id2(k,j+1,i)] += src_vol/emf->sigma22[k][j+1][i]*(1.-w1)*w2*(1.-w3);
	by[id2(k,j+1,i+1)] += src_vol/emf->sigma22[k][j+1][i+1]*w1*w2*(1.-w3);
	by[id2(k+1,j,i)] += src_vol/emf->sigma22[k+1][j][i]*(1.-w1)*(1.-w2)*w3;
	by[id2(k+1,j,i+1)] += src_vol/emf->sigma22[k+1][j][i+1]*w1*(1.-w2)*w3;
	by[id2(k+1,j+1,i)] += src_vol/emf->sigma22[k+1][j+1][i]*(1.-w1)*w2*w3;
	by[id2(k+1,j+1,i+1)] += src_vol/emf->sigma22[k+1][j+1][i+1]*w1*w2*w3;
      }
    } else if(strcmp(emf->chsrc[ic],"Ez")==0){
      //Jz(i,j,k+0.5)
      for(isrc=0; isrc<acq->nsrc; isrc++){
	i = find_index(n1+1, emf->x1, acq->src_x1[isrc]);
	j = find_index(n2+1, emf->x2, acq->src_x2[isrc]);
	k = find_index(n3, emf->x3s, acq->src_x3[isrc]);
	
	w1 = (acq->src_x1[isrc] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
	w2 = (acq->src_x2[isrc] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
	w3 = (acq->src_x3[isrc] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);
	vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3s[k+1]- emf->x3s[k]);
	src_vol = 1./vol;
	if(emf->fvm) src_vol = emf->I_omega_mu0;

	//Jz is continous, interpolate directly for Jz
	if(emf->airbc) k -= emf->nb;
	bz[id3(k,j,i)] += src_vol*(1.-w1)*(1.-w2)*(1.-w3);
	bz[id3(k,j,i+1)] += src_vol*w1*(1.-w2)*(1.-w3);
	bz[id3(k,j+1,i)] += src_vol*(1.-w1)*w2*(1.-w3);
	bz[id3(k,j+1,i+1)] += src_vol*w1*w2*(1.-w3);
	bz[id3(k+1,j,i)] += src_vol*(1.-w1)*(1.-w2)*w3;
	bz[id3(k+1,j,i+1)] += src_vol*w1*(1.-w2)*w3;
	bz[id3(k+1,j+1,i)] += src_vol*(1.-w1)*w2*w3;
	bz[id3(k+1,j+1,i+1)] += src_vol*w1*w2*w3;
      }
    }
  
  }

}


void extract_emf(acq_t *acq, emf_t *emf, complex *x, int ifreq)
{
  int i, j, k, kk, ic, irec;
  double w1, w2, w3;
  complex *Ex, *Ey, *Ez;
  complex s;

  int n1 = emf->nx + 2*emf->nb;
  int n2 = emf->ny + 2*emf->nb;
  int n3 = emf->nz + 2*emf->nb;
  Ex = &x[0];
  Ey = &x[emf->lEx];
  Ez = &x[emf->lEx + emf->lEy];
  
  for(ic=0; ic<emf->nchrec; ic++){
    if(strcmp(emf->chrec[ic],"Ex")==0){
      //Ex(i+0.5,j,k)
      for(irec=0; irec<acq->nrec; irec++){
	i = find_index(n1, emf->x1s, acq->rec_x1[irec]);
	j = find_index(n2+1, emf->x2, acq->rec_x2[irec]);
	k = find_index(n3+1, emf->x3, acq->rec_x3[irec]);
	
	w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
	w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
	w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
	if(emf->airbc) k -= emf->nb;

	s = 0;
	s += Ex[id1(k,j,i)]*(1.-w1)*(1.-w2)*(1.-w3);
	s += Ex[id1(k,j,i+1)]*w1*(1.-w2)*(1.-w3);
	s += Ex[id1(k,j+1,i)]*(1.-w1)*w2*(1.-w3);
	s += Ex[id1(k,j+1,i+1)]*w1*w2*(1.-w3);
	s += Ex[id1(k+1,j,i)]*(1.-w1)*(1.-w2)*w3;
	s += Ex[id1(k+1,j,i+1)]*w1*(1.-w2)*w3;
	s += Ex[id1(k+1,j+1,i)]*(1.-w1)*w2*w3;
	s += Ex[id1(k+1,j+1,i+1)]*w1*w2*w3;
	emf->dcal_fd[ic][ifreq][irec] = s;
      }
    }else if(strcmp(emf->chrec[ic],"Ey")==0){
      //Ey(i,j+0.5,k)
      for(irec=0; irec<acq->nrec; irec++){
	i = find_index(n1+1, emf->x1, acq->rec_x1[irec]);
	j = find_index(n2, emf->x2s, acq->rec_x2[irec]);
	k = find_index(n3+1, emf->x3, acq->rec_x3[irec]);

	w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
	w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
	w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
	if(emf->airbc) k -= emf->nb;

	s = 0;
	s += Ey[id2(k,j,i)]*(1.-w1)*(1.-w2)*(1.-w3);
	s += Ey[id2(k,j,i+1)]*w1*(1.-w2)*(1.-w3);
	s += Ey[id2(k,j+1,i)]*(1.-w1)*w2*(1.-w3);
	s += Ey[id2(k,j+1,i+1)]*w1*w2*(1.-w3);
	s += Ey[id2(k+1,j,i)]*(1.-w1)*(1.-w2)*w3;
	s += Ey[id2(k+1,j,i+1)]*w1*(1.-w2)*w3;
	s += Ey[id2(k+1,j+1,i)]*(1.-w1)*w2*w3;
	s += Ey[id2(k+1,j+1,i+1)]*w1*w2*w3;
	emf->dcal_fd[ic][ifreq][irec] = s;
      }
    } else if(strcmp(emf->chrec[ic],"Ez")==0){
      //Ez(i,j,k+0.5)
      for(irec=0; irec<acq->nrec; irec++){
	i = find_index(n1+1, emf->x1, acq->rec_x1[irec]);
	j = find_index(n2+1, emf->x2, acq->rec_x2[irec]);
	k = find_index(n3, emf->x3s, acq->rec_x3[irec]);

	w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
	w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
	w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);
	kk = k;
	if(emf->airbc) k -= emf->nb;

	//interpolation over Jz which is continuous in the presence of interface
	s = 0;
	s += emf->sigma33[kk][j][i]*Ez[id3(k,j,i)]*(1.-w1)*(1.-w2)*(1.-w3);
	s += emf->sigma33[kk][j][i+1]*Ez[id3(k,j,i+1)]*w1*(1.-w2)*(1.-w3);
	s += emf->sigma33[kk][j+1][i]*Ez[id3(k,j+1,i)]*(1.-w1)*w2*(1.-w3);
	s += emf->sigma33[kk][j+1][i+1]*Ez[id3(k,j+1,i+1)]*w1*w2*(1.-w3);
	s += emf->sigma33[kk+1][j][i]*Ez[id3(k+1,j,i)]*(1.-w1)*(1.-w2)*w3;
	s += emf->sigma33[kk+1][j][i+1]*Ez[id3(k+1,j,i+1)]*w1*(1.-w2)*w3;
	s += emf->sigma33[kk+1][j+1][i]*Ez[id3(k+1,j+1,i)]*(1.-w1)*w2*w3;
	s += emf->sigma33[kk+1][j+1][i+1]*Ez[id3(k+1,j+1,i+1)]*w1*w2*w3;
	//convert Jz to Ez: Ez=Jz/sigma_water, receivers sit in water
	emf->dcal_fd[ic][ifreq][irec] = s*emf->rho_water;
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
