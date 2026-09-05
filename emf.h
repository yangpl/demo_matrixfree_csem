#ifndef emf_h
#define emf_h

typedef struct {
  int mode; //mode=0, modeling; mode=1, FWI; mode=2, FWI gradient only
  int verb;/* verbose display */
  int reciprocity;
  int fvm;//1=FVM; 0=FDM
  int airbc;//1=air-water/Earth boundary condition; 0=not
  
  int nfreq;//number of frequencies
  float *freqs;//a list of frequencies
  float rho_air;//resistivity in air
  float lextend;//domain extension
  float rho_water;//resistivity in sea water
  float rhomin, rhomax;
  
  int nchsrc, nchrec;//number of active src/rec channels for E and H
  char **chsrc, **chrec;
  int istretch;//1=stretch; 0=no stretching
  float rstretch;//growth rate, i.e. stretching factor
  
  int nx, ny, nz, nb;  //size of the model defined on regular FD grid
  int n1, n2, n3;
  int lv, lE, lH, lEx, lEy, lEz, lHx, lHy, lHz;
  float *x1node, *x2node, *x3node;
  double *x1, *x2, *x3;
  double *x1s, *x2s, *x3s;
  double *d1, *d2, *d3;
  double *d1s, *d2s, *d3s;
  double d1min, d2min, d3min;
  double x1min, x2min, x3min;
  double x1max, x2max, x3max;

  int n1fft, n2fft;
  double d1fft, d2fft, d3fft;
  complex **sH1kxky, **sH2kxky;
  
  float ***rho11, ***rho22, ***rho33;
  double ***sigma11, ***sigma22, ***sigma33;
  float _Complex ****E1, ****E2, ****E3;
  float _Complex ***dcal_fd;
  
  int ne, nh;//length of vector E=(Ex,Ey,Ez)^T
  complex *xe, *xh, *be;
  complex *v, *bv;
  complex *xx, *yy;
  complex I_omega_mu0;
  int preco;//1=precondiiton; 0=not
  int algopt;//1=A-V using LU-SGS; 2=line Gauss-Seidel in x
  int niter;//total number of iterations
  int nrestart;//number of iterations for restart
  int mp;//GMRES(mp) preconditioner
  float tol;//convergence tolerance

  int param_extended; //1=output parameters after domain extension
} emf_t; /* type of electromagnetic field (emf)  */

#endif
