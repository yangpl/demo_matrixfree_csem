/* Airwave manipulation using Fourier transform at the top boundary
 *-----------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Anothr: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

fftw_complex *emf_kxky, *emf_kxkyz0;
fftw_plan fft_airwave, ifft_airwave;
int *jx, *jy;
double *wx, *wy;

int find_index(int n, double *x, double val);

/*--------------------------------------------------------------------------*/
void airbc_init(emf_t *emf)
{
  double dkx, dky, kx, ky, kz, x1, x2;
  int i1, i2, j1, j2, istart, ntaper;
  double *win1, *win2;
  
  if(!getparint("n1fft", &emf->n1fft)) emf->n1fft = 1024; //fft_next_fast_size((emf->x1[emf->n1] - emf->x1[0])/emf->d1fft + 1);
  if(!getparint("n2fft", &emf->n2fft)) emf->n2fft = 1024; //fft_next_fast_size((emf->x2[emf->n2] - emf->x2[0])/emf->d2fft + 1);  
  emf->d1fft = (emf->x1[emf->n1] - emf->x1[0])/(emf->n1fft-1);
  emf->d2fft = (emf->x2[emf->n2] - emf->x2[0])/(emf->n2fft-1);
  emf->d3fft = emf->x3node[1] - emf->x3node[0];
  if(iproc==0){
    printf("[n1fft, n2fft]=[%d, %d]\n", emf->n1fft, emf->n2fft);
    printf("[d1fft, d2fft, d3fft]=[%g, %g, %g]\n", emf->d1fft, emf->d2fft, emf->d3fft);
  }

  emf->sH1kxky = alloc2complex(emf->n1fft, emf->n2fft);//scaling factor for H1
  emf->sH2kxky = alloc2complex(emf->n1fft, emf->n2fft);//scaling factor for H2
  
  /* FE3 is not necessary in the air because we do not compute derivates of Hx
   * and Hy in the air: Hx and Hy are derived directly by extrapolation from Hz. */
  emf_kxky = fftw_malloc(sizeof(fftw_complex)*emf->n1fft*emf->n2fft);
  emf_kxkyz0 = fftw_malloc(sizeof(fftw_complex)*emf->n1fft*emf->n2fft);
  /* comapred with FFTW, we have opposite sign convention for time, same sign convetion for space */
  fft_airwave = fftw_plan_dft_2d(emf->n1fft, emf->n2fft, emf_kxky, emf_kxky, FFTW_FORWARD, FFTW_ESTIMATE);
  ifft_airwave = fftw_plan_dft_2d(emf->n1fft, emf->n2fft, emf_kxky, emf_kxky, FFTW_BACKWARD, FFTW_ESTIMATE);  

  /* win1 = alloc1double(emf->n1fft); */
  /* win2 = alloc1double(emf->n2fft); */

  /* ntaper = 0.1*emf->n1fft;// number of points for tapering */
  /* istart = emf->n1fft/2 - ntaper; */
  /* win1[0] = 1.; */
  /* for(i1=1; i1<=emf->n1fft/2; i1++){ */
  /*   win1[i1] = (i1>istart)?cos(0.5*PI*(i1-istart)/(ntaper+1)):1; */
  /*   win1[emf->n1fft-i1] = win1[i1];//mirror window function along Nyquist  */
  /* } */
  /* ntaper = 0.1*emf->n2fft;// number of points for tapering */
  /* istart = emf->n2fft/2 - ntaper; */
  /* win2[0] = 1.; */
  /* for(i2=1; i2<=emf->n2fft/2; i2++){ */
  /*   win2[i2] = (i2>istart)?cos(0.5*PI*(i2-istart)/(ntaper+1)):1; */
  /*   win2[emf->n2fft-i2] = win2[i2];//mirror window function along Nyquist */
  /* } */
  
  dkx = 2.0*PI/(emf->d1fft*emf->n1fft);
  dky = 2.0*PI/(emf->d2fft*emf->n2fft);
  for(i2=0; i2<emf->n2fft; i2++){
    ky = (i2<=emf->n2fft/2)?i2*dky:(i2-emf->n2fft)*dky;
    for(i1=0; i1<emf->n1fft; i1++){
      kx = (i1<=emf->n1fft/2)?i1*dkx:(i1-emf->n1fft)*dkx;

      kz = sqrt(kx*kx + ky*ky);
      emf->sH1kxky[i2][i1] = exp(-kz*0.5*emf->d3fft)*I*kx/(kz+1.e-15);
      emf->sH2kxky[i2][i1] = exp(-kz*0.5*emf->d3fft)*I*ky/(kz+1.e-15);
      emf->sH1kxky[i2][i1] *= cexp(-I*kx*0.5*emf->d1s[0]);//shift half grid due to staggering
      emf->sH2kxky[i2][i1] *= cexp(-I*ky*0.5*emf->d2s[0]);//shift half grid due to staggering
      /* emf->sH1kxky[i2][i1] *= win1[i1]*win2[i2]; */
      /* emf->sH2kxky[i2][i1] *= win1[i1]*win2[i2]; */
    }
  }

  jx = alloc1int(emf->n1fft);
  jy = alloc1int(emf->n2fft);
  wx = alloc1double(emf->n1fft);
  wy = alloc1double(emf->n2fft);
  for(i1=0; i1<emf->n1fft; i1++){
    x1 = emf->x1s[0] + i1*emf->d1fft;
    if(x1<emf->x1s[emf->n1-1]){
      j1 = find_index(emf->n1+1, emf->x1s, x1);
      jx[i1] = j1;
      wx[i1] = (x1-emf->x1s[j1])/(emf->x1s[j1+1]-emf->x1s[j1]);
    }
  }
  for(i2=0; i2<emf->n2fft; i2++){
    x2 = emf->x2s[0] + i2*emf->d2fft;
    if(x2<emf->x2s[emf->n2-1]){
      j2 = find_index(emf->n2+1, emf->x2s, x2);
      jy[i2] = j2;
      wy[i2] = (x2-emf->x2s[j2])/(emf->x2s[j2+1]-emf->x2s[j2]);
    }
  }
  
  /* free1double(win1); */
  /* free1double(win2);    */
}

void airbc_free(emf_t *emf)
{
  free2complex(emf->sH1kxky);
  free2complex(emf->sH2kxky);
  fftw_free(emf_kxky);
  fftw_free(emf_kxkyz0);
  fftw_destroy_plan(fft_airwave);
  fftw_destroy_plan(ifft_airwave);

  free1int(jx);
  free1int(jy);
  free1double(wx);
  free1double(wy);
}

//H1[i,j+0.5,k+0.5], size=(n1+1)*n2*n3;
//H2[i+0.5,j,k+0.5], size=n1*(n2+1)*n3;
//H3[i+0.5,j+0.5,k], size=n1*n2*(n3+1);
//H3k=H3(z=0); H1km1=H1(z=-0.5*dz), H2km1=H2(z=-0.5*dz)
void airbc_apply(emf_t *emf, complex *H3k, complex *H1km1, complex *H2km1)
{
  int i1, i2;
  int j1, j2;
  double x1, x2, w1, w2, s1, s2;

  //interpolate Hz from nonuniform grid to uniform grid
  for(i2=0; i2<emf->n2fft; i2++){
    x2 = emf->x2s[0] + i2*emf->d2fft;
    for(i1=0; i1<emf->n1fft; i1++){
      x1 = emf->x1s[0] + i1*emf->d1fft;

      if(x1<emf->x1s[emf->n1-1] && x2<emf->x2s[emf->n2-1]){//bilinear interpolation
	j1 = jx[i1];
	j2 = jy[i2];
	w1 = wx[i1];
	w2 = wy[i2];
	emf_kxky[i1+emf->n1fft*i2] = (1.-w1)*(1.-w2)*H3k[j1+emf->n1*j2] + w1*(1.-w2)*H3k[j1+1+emf->n1*j2]
	  + (1.-w1)*w2*H3k[j1+emf->n1*(j2+1)] + w1*w2*H3k[j1+1+emf->n1*(j2+1)];
      }else
	emf_kxky[i1+emf->n1fft*i2] = 0.; 
    }
  }
  fftw_execute(fft_airwave);//Hz(x,y,z=0)-->Hz(kx,ky,z=0)
  memcpy(emf_kxkyz0, emf_kxky, emf->n1fft*emf->n2fft*sizeof(fftw_complex));
  for(i2=0; i2<emf->n2fft; i2++){
    for(i1=0; i1<emf->n1fft; i1++){
      //Hz-->Hx, in wavenumber domain
      emf_kxky[i1+emf->n1fft*i2] = emf_kxkyz0[i1+emf->n1fft*i2]*emf->sH1kxky[i2][i1];
    }
  }
  fftw_execute(ifft_airwave);
  for(i2=0; i2<emf->n2fft; i2++){
    for(i1=0; i1<emf->n1fft; i1++){
      emf_kxky[i1+emf->n1fft*i2] /= (emf->n1fft*emf->n2fft);//normalize FFT
    }
  }
  //get Hx[i,j+0.5,k+0.5] on nonuniform grid from uniform grid via bilinear interpolation
  for(i2=0; i2<emf->n2; i2++){
    s2 = (emf->x2s[i2]-emf->x2s[0]);
    j2 = s2/emf->d2fft;
    w2 = s2/emf->d2fft-j2;
    for(i1=0; i1<=emf->n1; i1++){
      s1 = (emf->x1[i1]-emf->x1[0]);
      j1 = s1/emf->d1fft;
      w1 = s1/emf->d1fft-j1;
      //output Hx space domain
      H1km1[i1+(emf->n1+1)*i2] = (1.-w1)*(1.-w2)*emf_kxky[j1+emf->n1fft*j2] + w1*(1.-w2)*emf_kxky[j1+1+emf->n1fft*j2]
	+ (1.-w1)*w2*emf_kxky[j1+emf->n1fft*(j2+1)] + w1*w2*emf_kxky[j1+1+emf->n1fft*(j2+1)];
    }
  }
  
  for(i2=0; i2<emf->n2fft; i2++){
    for(i1=0; i1<emf->n1fft; i1++){
      //Hz-->Hy, in wave number domain
      emf_kxky[i1+emf->n1fft*i2] = emf_kxkyz0[i1+emf->n1fft*i2]*emf->sH2kxky[i2][i1];
    }
  }
  fftw_execute(ifft_airwave);
  for(i2=0; i2<emf->n2fft; i2++){
    for(i1=0; i1<emf->n1fft; i1++){
      emf_kxky[i1+emf->n1fft*i2] /= (emf->n1fft*emf->n2fft);//normalize FFT
    }
  }    
  //get Hy[i+0.5,j,k+0.5] on nonuniform grid from uniform grid by bilinear interpolation
  for(i2=0; i2<=emf->n2; i2++){
    s2 = (emf->x2[i2]-emf->x2[0]);
    j2 = s2/emf->d2fft;
    w2 = s2/emf->d2fft-j2;
    for(i1=0; i1<emf->n1; i1++){
      s1 = (emf->x1s[i1]-emf->x1s[0]);
      j1 = s1/emf->d1fft;
      w1 = s1/emf->d1fft-j1;
      //output Hy space domain
      H2km1[i1+emf->n1*i2] = (1.-w1)*(1.-w2)*emf_kxky[j1+emf->n1fft*j2] + w1*(1.-w2)*emf_kxky[j1+1+emf->n1fft*j2]
	+ (1.-w1)*w2*emf_kxky[j1+emf->n1fft*(j2+1)] + w1*w2*emf_kxky[j1+1+emf->n1fft*(j2+1)];
    }
  }    
  
}

