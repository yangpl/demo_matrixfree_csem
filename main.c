/* Preconditioned BiCGStab/GMRES for solving Maxwell equation
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "emf.h"
#include "acq.h"
#include "solver.h"

int iproc, nproc, ierr;

void emf_init(emf_t *emf);
void emf_free(emf_t *emf);

void acq_init(acq_t *acq, emf_t * emf);
void acq_free(acq_t *acq);

void inject_source(acq_t *acq, emf_t *emf, complex *b);
void extract_emf(acq_t *acq, emf_t *emf, complex *x, int ifreq);
void write_data(acq_t *acq, emf_t *emf, char *fname, float _Complex ***dcal_fd);

void extend_model_init(acq_t *acq, emf_t *emf);
void extend_model_free(emf_t *emf);

void fdfd_init(emf_t *emf);
void fdfd_free(emf_t *emf);
void fdfd_apply(int n, complex *x, complex *y);
void lusgs_apply(int n, complex *x, complex *y);
void aphi_apply(int n, complex *x, complex *y);

void gmres_init(int n, int m, op_t Aop);
void gmres_free();
void gmres_apply(int n, complex *x, complex *y);

void bicgstab(int n, complex *x, complex *b, op_t Aop, int niter, double tol, int verb);
void gmres(int n, complex *x, complex *b, op_t Aop, int niter, double tol, int m, int verb);
void fbicgstab(int n, complex *x, complex *b, op_t Aop, op_t invMop, int niter, double tol, int verb);
void fgmres(int n, complex *x, complex *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb);

//matrix-free implementation
int main(int argc, char **argv)
{
  emf_t *emf;
  acq_t *acq;
  int ifreq;
  float tol, *time;
  char fname[sizeof("emf_0000.txt")];
  
  MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  initargs(argc, argv);

  emf = malloc(sizeof(emf_t));
  acq = malloc(sizeof(acq_t));
  if(!getparint("mode", &emf->mode)) emf->mode = 0;
  if(!getparint("verb", &emf->verb)) emf->verb = (iproc==0)?1:0;
  if(!getparint("niter", &emf->niter)) emf->niter = 1000;
  if(!getparint("preco", &emf->preco)) emf->preco = 1;//0=no precondition; 1=A-phi
  if(!getparint("algopt", &emf->algopt)) emf->algopt = 1;//1=BICSTAB, 2=GMRES
  if(!getparfloat("tol", &emf->tol)) emf->tol = 1e-8;
  if(!getparint("nrestart", &emf->nrestart)) emf->nrestart = 10;
  if(!getparint("mp", &emf->mp)) emf->mp = 30;//m for preconditioner
  if(emf->verb){
    printf("preco=%d (1=precondition; 0=unpreconditioned)\n", emf->preco);
    printf("algopt=%d (1=BICGStab; 2=GMRES)\n", emf->algopt);
    printf("niter=%d\n", emf->niter);
    printf("nrestart=%d\n", emf->nrestart);
    printf("mp=%d (GMRES(mp) preconditioner)\n", emf->mp);
    printf("tol=%.2e\n", emf->tol);
  }
  
  emf_init(emf);
  acq_init(acq, emf);
  extend_model_init(acq, emf);
  fdfd_init(emf);
  if(emf->verb) printf("dof=%d\n", emf->lE);
  emf->ne = emf->lE;
  emf->nh = emf->lH;
  emf->xe = alloc1complex(emf->ne);//vector xe=(Ex,Ey,Ez)^T
  emf->xh = alloc1complex(emf->nh);//vector xh=(Hx,Hy,Hz)^T
  emf->be = alloc1complex(emf->ne);//vector be=(bx,by,bz)^T
  emf->dcal_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  memset(&emf->dcal_fd[0][0][0], 0, acq->nrec*emf->nfreq*emf->nchrec*sizeof(float _Complex));
  time = malloc(emf->nfreq*sizeof(float));
  memset(time, 0, emf->nfreq*sizeof(float));
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    tol = emf->tol/emf->freqs[ifreq];//frequency
    emf->I_omega_mu0 = I*2.*PI*emf->freqs[ifreq]*mu0;
    if(iproc==0) printf("------freq=%g Hz, tol=%e--------\n", emf->freqs[ifreq], tol);

    memset(emf->xe, 0, emf->ne*sizeof(complex));//initialization
    memset(emf->xh, 0, emf->nh*sizeof(complex));//initialization
    memset(emf->be, 0, emf->ne*sizeof(complex));//initialization        
    inject_source(acq, emf, emf->be);

    if(iproc==0) time[ifreq] = MPI_Wtime();//record start time
    if(emf->preco==0){
      if(emf->algopt==1) bicgstab(emf->ne, emf->xe, emf->be, fdfd_apply, emf->niter, tol, emf->verb);
      if(emf->algopt==2) gmres(emf->ne, emf->xe, emf->be, fdfd_apply, emf->niter/emf->nrestart, tol, emf->nrestart, emf->verb);
    }else if(emf->preco==1){
      if(emf->algopt==1) fbicgstab(emf->ne, emf->xe, emf->be, fdfd_apply, lusgs_apply, emf->niter, tol, emf->verb);
      if(emf->algopt==2) fgmres(emf->ne, emf->xe, emf->be, fdfd_apply, lusgs_apply, emf->niter/emf->nrestart, tol, emf->nrestart, emf->verb);
    }else if(emf->preco==2){
      gmres_init(emf->ne, emf->mp, fdfd_apply);
      if(emf->algopt==1) fbicgstab(emf->ne, emf->xe, emf->be, fdfd_apply, gmres_apply, emf->niter, tol, emf->verb);
      if(emf->algopt==2) fgmres(emf->ne, emf->xe, emf->be, fdfd_apply, gmres_apply, emf->niter/emf->nrestart, tol, emf->nrestart, emf->verb);
      gmres_free();      
    }else if(emf->preco==3){
      if(emf->algopt==1) fbicgstab(emf->ne, emf->xe, emf->be, fdfd_apply, aphi_apply, emf->niter, tol, emf->verb);
      if(emf->algopt==2) fgmres(emf->ne, emf->xe, emf->be, fdfd_apply, aphi_apply, emf->niter/emf->nrestart, tol, emf->nrestart, emf->verb);
    }    
    if(iproc==0) time[ifreq] = MPI_Wtime()-time[ifreq];//compute elapsed time for ifreq-th freq
    if(iproc==0) printf("elapsed time: %g sec\n", time[ifreq]);
    
    extract_emf(acq, emf, emf->xe, ifreq);
  }
  if(iproc==0){
    sprintf(fname, "emf_%04d.txt", 1);//only 1st source
    write_data(acq, emf, fname, emf->dcal_fd);
  }

  free1complex(emf->xe);
  free1complex(emf->xh);
  free1complex(emf->be);
  free(time);
  free3complexf(emf->dcal_fd);
  fdfd_free(emf);
  extend_model_free(emf);
  emf_free(emf);
  acq_free(acq);
  
  free(emf);
  free(acq);

  MPI_Finalize();

  return 0;
}

//=====================================================
void Acsr_init(emf_t *emf, int ifreq);
void Acsr_free();
void Acsr_apply(int n, complex *x, complex *y);
void Acsr_lusgs_apply(int n, complex *x, complex *y);

//matrix-storage implementation
int main2(int argc, char **argv)
{
  emf_t *emf;
  acq_t *acq;
  int ifreq;
  char fname[sizeof("emf_0000.txt")];
  
  MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  initargs(argc, argv);

  emf = malloc(sizeof(emf_t));
  acq = malloc(sizeof(acq_t));
  if(!getparint("mode", &emf->mode)) emf->mode = 0;
  if(!getparint("verb", &emf->verb)) emf->verb = (iproc==0)?1:0;
  
  emf_init(emf);
  acq_init(acq, emf);
  extend_model_init(acq, emf);
  emf->ne = emf->lE;
  emf->xe = alloc1complex(emf->ne);//vector x=(Ex,Ey,Ez)^T
  emf->be = alloc1complex(emf->ne);//vector b=(bx,by,bz)^T
  emf->dcal_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  memset(&emf->dcal_fd[0][0][0], 0, acq->nrec*emf->nfreq*emf->nchrec*sizeof(float _Complex));
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    if(emf->verb) printf("----freq=%g Hz--------\n", emf->freqs[ifreq]);

    Acsr_init(emf, ifreq);
    memset(emf->xe, 0, emf->ne*sizeof(complex));//initialization
    memset(emf->be, 0, emf->ne*sizeof(complex));//initialization        
    inject_source(acq, emf, emf->be);

    if(emf->preco==0){
      if(emf->algopt==1) bicgstab(emf->ne, emf->xe, emf->be, Acsr_apply, emf->niter, emf->tol, emf->verb);
      if(emf->algopt==2) gmres(emf->ne, emf->xe, emf->be, Acsr_apply, emf->niter/emf->nrestart, emf->tol, emf->nrestart, emf->verb);
    }else if(emf->preco==1){
      if(emf->algopt==1) fbicgstab(emf->ne, emf->xe, emf->be, Acsr_apply, Acsr_lusgs_apply, emf->niter, emf->tol, emf->verb);
      if(emf->algopt==2) fgmres(emf->ne, emf->xe, emf->be, Acsr_apply, Acsr_lusgs_apply, emf->niter/emf->nrestart, emf->tol, emf->nrestart, emf->verb);
    }
    
    Acsr_free();
    extract_emf(acq, emf, emf->xe, ifreq);

  }
  if(iproc==0){
    sprintf(fname, "emf_%04d.txt", 1);//only 1st source
    write_data(acq, emf, fname, emf->dcal_fd);
  }
  
  free1complex(emf->xe);
  free1complex(emf->be);
  free3complexf(emf->dcal_fd);
  extend_model_free(emf);
  emf_free(emf);
  acq_free(acq);
  
  free(emf);
  free(acq);

  MPI_Finalize();

  return 0;
}

//=====================================================
void Acoo_init(emf_t *emf, int ifreq);
void Acoo_free();
void Acoo_solve(emf_t *emf);

//matrix-storage implementation
int main3(int argc, char **argv)
{
  emf_t *emf;
  acq_t *acq;
  int ifreq;
  char fname[sizeof("emf_0000.txt")];
  
  MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  initargs(argc, argv);

  emf = malloc(sizeof(emf_t));
  acq = malloc(sizeof(acq_t));
  if(!getparint("mode", &emf->mode)) emf->mode = 0;
  if(!getparint("verb", &emf->verb)) emf->verb = (iproc==0)?1:0;
  
  emf_init(emf);
  acq_init(acq, emf);
  extend_model_init(acq, emf);
  emf->ne = emf->lE;
  emf->be = alloc1complex(emf->ne);//vector b=(bx,by,bz)^T
  emf->dcal_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  memset(&emf->dcal_fd[0][0][0], 0, acq->nrec*emf->nfreq*emf->nchrec*sizeof(float _Complex));
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    if(emf->verb) printf("----freq=%g Hz--------\n", emf->freqs[ifreq]);

    Acoo_init(emf, ifreq);
    memset(emf->be, 0, emf->ne*sizeof(complex));//initialization        
    inject_source(acq, emf, emf->be);
    for(int k=0; k<emf->ne; k++) emf->be[k] *= emf->I_omega_mu0;//Js-->i*omega*mu0*Js

    Acoo_solve(emf);
    Acoo_free();
    extract_emf(acq, emf, emf->be, ifreq);
  }
  if(iproc==0){
    sprintf(fname, "emf_%04d.txt", 1);
    write_data(acq, emf, fname, emf->dcal_fd);
  }

  free1complex(emf->be);
  free3complexf(emf->dcal_fd);
  extend_model_free(emf);
  emf_free(emf);
  acq_free(acq);
  
  free(emf);
  free(acq);

  MPI_Finalize();

  return 0;
}

