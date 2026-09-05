/* linear solvers - BICGSTAB and GMRES with right preconditioning
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "solver.h"

//compute dot product s=<x,y>=x^T*conj(y), this definiton is required according to Saad book
complex cdotprod(int n, complex *a, complex *b);

//linear solver using BiCGStab, algorithm 7.7 in Saad book
void bicgstab(int n, complex *x, complex *b, op_t Aop, int niter, double tol, int verb)
{
  int i, iter;
  double rs0, rs, relres;
  complex rho_old, rho_new, alpha, beta, omega;
  FILE *fp;
  complex *r = malloc(n*sizeof(complex));
  complex *r0 = malloc(n*sizeof(complex));//rprime0
  complex *p = malloc(n*sizeof(complex));
  complex *v = malloc(n*sizeof(complex));
  complex *s = malloc(n*sizeof(complex));
  complex *t = malloc(n*sizeof(complex));
  
  Aop(n, x, v);//v=Ax
  for(i=0; i<n; i++) {
    r[i] = b[i]-v[i];//r=b-Ax
    p[i] = r[i];
    r0[i] = r[i];
  }
  rho_old = cdotprod(n, r, r0);
  rs = creal(cdotprod(n, r, r));
  rs0 = rs;
  if(rs==0.) return;
  relres = sqrt(rs/rs0);//relative residual
  if(verb){
    fp = fopen("iterate_bicgstab.txt", "w");
    fprintf(fp, "iteration \t relres\n");
    fprintf(fp, "%d \t %e\n", 0, relres);
    fclose(fp);
    printf("bicgstab iter=%d relres=%e\n", 0, relres);
  }
  
  for(iter=0; iter<niter; ++iter){
    if(verb && iter>0){
      fp = fopen("iterate_bicgstab.txt", "a");
      fprintf(fp, "%d \t %e\n", iter, relres);
      fclose(fp);
      printf("bicgstab iter=%d relres=%e\n", iter, relres);
    }

    Aop(n, p, v);//v=Ap
    alpha = rho_old/cdotprod(n, v, r0);
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];

    Aop(n, s, t);//t=As
    omega = cdotprod(n, t, s)/creal(cdotprod(n, t, t));

    for(i=0; i<n; i++){
      x[i] += alpha*p[i] + omega*s[i];
      r[i] = s[i] - omega*t[i];
    }
    rs = creal(cdotprod(n, r, r));
    relres = sqrt(rs/rs0);
    if(relres<tol) {
      if(verb) printf("converged at iter=%d\n", iter);
      break;
    }

    rho_new = cdotprod(n, r, r0);
    beta = (rho_new/rho_old)*alpha/omega;
    for(i=0; i<n; i++) p[i] = r[i] + beta*(p[i]-omega*v[i]);

    rho_old = rho_new;
  }
    
  free(r);
  free(r0);
  free(p);
  free(v);
  free(s);
  free(t);
}

//Reference:
//[1] Jie Chen, 2016 JSC, Right preconditioned/Flexible BiCGStab paper
//[2] Vogel, Judith A. "Flexible BiCG and flexible Bi-CGSTAB for nonsymmetric linear systems."
//    Applied Mathematics and Computation 188.1 (2007): 226-233.
void fbicgstab(int n, complex *x, complex *b, op_t Aop, op_t invMop, int niter, double tol, int verb)
{
  int i, iter;
  double rs0, rs, relres;
  complex rho_old, rho_new, alpha, beta, omega;
  FILE *fp;
  complex *r = malloc(n*sizeof(complex));
  complex *r0 = malloc(n*sizeof(complex));
  complex *p = malloc(n*sizeof(complex));
  complex *pp = malloc(n*sizeof(complex));
  complex *v = malloc(n*sizeof(complex));
  complex *s = malloc(n*sizeof(complex));
  complex *t = malloc(n*sizeof(complex));
  
  Aop(n, x, v);//v=Ax
  for(i=0; i<n; i++) {
    r[i] = b[i]-v[i];//r=b-Ax
    p[i] = r[i];
    r0[i] = r[i];
  }
  rho_old = cdotprod(n, r, r0);
  rs = creal(cdotprod(n, r, r));
  rs0 = rs;
  if(rs==0.){
    free(r);
    free(r0);
    free(p);
    free(pp);
    free(v);
    free(s);
    free(t);
    return;
  }
  relres = sqrt(rs/rs0);
  if(verb) {
    fp = fopen("iterate_fbicgstab.txt", "w");
    fprintf(fp, "iteration \t relres\n");
    fprintf(fp, "%d \t %e\n", iter, relres);
    fclose(fp);
    printf("f-bicgstab iter=%d relres=%e\n", 0, relres);
  }
  
  for(iter=0; iter<niter; ++iter){
    if(iter>0 && verb) {
      fp = fopen("iterate_fbicgstab.txt", "a");
      fprintf(fp, "%d \t %e\n", iter, relres);
      fclose(fp);
      printf("f-bicgstab iter=%d relres=%e\n", iter, relres);
    }

    invMop(n, p, pp);//pp=invM*p
    Aop(n, pp, v);//v=A*pp
    alpha = cdotprod(n, v, r0);
    alpha = rho_old/alpha;
    for(i=0; i<n; i++) s[i] = r[i] - alpha*v[i];
    
    invMop(n, s, r);//ss=invM*s
    Aop(n, r, t);//t=A*ss
    omega = cdotprod(n, s, t)/creal(cdotprod(n, t, t));

    for(i=0; i<n; i++){
      x[i] += alpha*pp[i] + omega*r[i];
      r[i] = s[i] - omega*t[i];
    }
    rs = creal(cdotprod(n, r, r));
    relres = sqrt(rs/rs0);
    if(relres<tol) {
      if(verb) printf("converged at iter=%d\n", iter);
      break;
    }
    
    rho_new = cdotprod(n, r, r0);
    beta = (rho_new/rho_old)*alpha/omega;
    for(i=0; i<n; i++) p[i] = r[i] + beta*(p[i]-omega*v[i]);

    rho_old = rho_new;
  }
  
  free(r);
  free(r0);
  free(p);
  free(pp);
  free(v);
  free(s);
  free(t);
}

//GMRES without preconditioning
void gmres(int n, complex *x, complex *b, op_t Aop, int niter, double tol, int m, int verb)
{
  int i, j, k, iter;
  double beta, tmp, r0;
  complex ss;
  FILE *fp;
  complex *w = alloc1complex(n);
  complex *r = alloc1complex(n);
  complex **v = alloc2complex(n, m+1);
  complex **h = alloc2complex(m, m+1);
  complex *g = alloc1complex(m+1);
  complex *y = alloc1complex(m+1);
  complex *c = alloc1complex(m+1);
  double *s = alloc1double(m+1);
  
  for(iter=0; iter<niter; iter++){
    Aop(n, x, w);//w=A*x
    for(i=0; i<n; i++) r[i] = b[i] -w[i];
    beta = sqrt(cdotprod(n, r, r));
    if(beta==0.0) return;
    for(i=0; i<n; i++) v[0][i] = r[i]/beta;
    memset(g, 0, (m+1)*sizeof(complex));
    g[0] = beta;
    memset(&h[0][0], 0, (m+1)*m*sizeof(complex));
    if(iter==0){
      r0 = beta;
      if(verb){
	fp = fopen("iterate_gmres.txt", "w");
	fprintf(fp, "iteration \t error\n");
  	fprintf(fp, "%d \t %e\n", iter, beta/r0);
	fclose(fp);
	printf("gmres iter=0 relres=%e\n", beta/r0);
      }
    }
    
    /*---------------------------------------------------*/
    for(j=0; j<m; j++){
      Aop(n, v[j], w);//r=Av;
      for(i=0; i<=j; i++) {
	h[i][j] = cdotprod(n, w, v[i]);
	for(k=0; k<n; k++) w[k] -= h[i][j]*v[i][k];
      }
      h[j+1][j] = sqrt(creal(cdotprod(n, w, w)));

      if(cabs(h[j+1][j])==0.0) { m=j+1; break; }
      for(i=0; i<n; i++) v[j+1][i] = w[i]/h[j+1][j];

      //solve least-squares problem by QR factorization using Givens rotations
      //min \|g - H(1:m,1:m) y\|^2, g=(beta,0,...,0)
      //min \|G*g - G*Hy\|^2, G=G(i-1,i,theta)*...*G(2,3,theta)*G(1,2,theta)
      if(j>0){
	for(i=0; i<j; i++){
	  //apply G12, G23,..., G_{j-1,j} to the last column of H_{j,*}
	  ss = conj(c[i])*h[i][j] + s[i]*h[i+1][j];
	  h[i+1][j] = -s[i]*h[i][j] + c[i]*h[i+1][j];
	  h[i][j] = ss;
	}
      }

      //compute c=cos(theta) and s=sin(theta)
      tmp = sqrt(h[j][j]*conj(h[j][j]) + h[j+1][j]*conj(h[j+1][j]));
      s[j] = creal(h[j+1][j])/tmp;
      c[j] = h[j][j]/tmp;
      h[j][j] = conj(c[j])*h[j][j] + s[j]*h[j+1][j];
      h[j+1][j] = 0.;
      //g=G(j,j+1,theta)g with g[j+1]=0
      g[j+1] = -s[j]*g[j];
      g[j] = conj(c[j])*g[j];

      tmp = cabs(g[j+1])/r0;
      if(tmp<tol) { m=j+1; return; }
      if(verb){
	fp = fopen("iterate_gmres.txt", "a");
	fprintf(fp, "%d \t %e\n", iter*m+j+1, tmp);
	fclose(fp);
	printf("gmres iter=%d relres=%e\n", iter*m+j+1, tmp);
      }
    }

    //now, H becomes an upper triangule matrix, problem min\|g-Hy\|^2 is g=Hy
    //solve it by backward substitution, y = H(1:m,1:m)\g(1:m)
    y[m-1] = g[m-1]/h[m-1][m-1];
    for(i=m-2; i>=0; i--){
      y[i] = g[i];
      for(j=i+1; j<m; j++) y[i] -= h[i][j]*y[j];
      y[i] /= h[i][i];
    }

    //x=x0+Vm*y
    for(j=0; j<m; j++){
      for(i=0; i<n; i++){
	x[i] += y[j]*v[j][i];
      }
    }      
  }//end for iter
  
  free1complex(r);
  free1complex(w);
  free2complex(v);
  free2complex(h);
  free1complex(g);
  free1complex(y);
  free1complex(c);
  free1double(s);
}


//GMRES with right preconditioning
void gmres_rp(int n, complex *x, complex *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb)
{
  int i, j, k, iter, mm, conv;
  double beta, tmp, r0;
  complex ss;
  FILE *fp;
  complex *w = alloc1complex(n);
  complex *r = alloc1complex(n);
  complex **v = alloc2complex(n, m+1);
  complex **h = alloc2complex(m+1, m);
  complex *g = alloc1complex(m+1);
  complex *y = alloc1complex(m+1);
  complex *c = alloc1complex(m+1);
  double *s = alloc1double(m+1);

  conv = 1;
  for(iter=0; iter<niter; iter++){
    /*-------------------------------------------------------------*/
    Aop(n, x, w);//w=A*x
    for(i=0; i<n; i++) r[i] = b[i]-w[i];
    beta = sqrt(creal(cdotprod(n, r, r)));
    if(beta==0.) break;
    for(i=0; i<n; i++) v[0][i] = r[i]/beta;
    memset(g, 0, (m+1)*sizeof(complex));
    g[0] = beta;
    memset(&h[0][0], 0, (m+1)*m*sizeof(complex));
    if(iter==0){
      r0 = beta;
      if(verb){
	fp = fopen("iterate_fgmres.txt", "w");
	fprintf(fp, "iteration \t error\n");
	fprintf(fp, "%d \t %e\n", iter, beta/r0);
	fclose(fp);
	printf("f-gmres iter=0 relres=%e\n", beta/r0);
      }
    }
    
    for(j=0; j<m; j++){
      invMop(n, v[j], r);//r=invM*v;
      Aop(n, r, w);//w=Ar;
      for(i=0; i<=j; i++) {
	h[j][i] = cdotprod(n, w, v[i]);
	for(k=0; k<n; k++) w[k] -= h[j][i]*v[i][k];
      }
      h[j][j+1] = sqrt(creal(cdotprod(n, w, w)));

      if(cabs(h[j][j+1])==0.) { conv = 1; break; }
      for(i=0; i<n; i++) v[j+1][i] = w[i]/h[j][j+1];

      //solve least-squares problem by QR factorization using Given rotations
      //min \|g - H(1:m,1:m) y\|^2, g=(beta,0,...,0)
      //min \|G*g - G*Hy\|^2, G=G(i-1,i,theta)*...*G(2,3,theta)*G(1,2,theta)
      if(j>0){
	for(i=0; i<j; i++){
	  //apply G12, G23,..., G_{j-1,j} to the last column of H_{j,*}
	  ss = conj(c[i])*h[j][i] + s[i]*h[j][i+1];
	  h[j][i+1] = -s[i]*h[j][i] + c[i]*h[j][i+1];
	  h[j][i] = ss;
	}
      }

      //compute c=cos(theta) and s=sin(theta)
      tmp = sqrt(creal(h[j][j]*conj(h[j][j]) + h[j][j+1]*conj(h[j][j+1])));
      s[j] = creal(h[j][j+1])/tmp;
      c[j] = h[j][j]/tmp;
      //g=G(j,j+1,theta)g with g[j+1]=0
      g[j+1] = -s[j]*g[j];
      g[j] = conj(c[j])*g[j];
      h[j][j] = conj(c[j])*h[j][j] + s[j]*h[j][j+1];

      tmp = cabs(g[j+1])/r0;
      if(verb){
	fp = fopen("iterate_fgmres.txt", "a");
	fprintf(fp, "%d \t %e\n", iter*m+j+1, tmp);
	fclose(fp);
	printf("f-gmres iter=%d relres=%e\n", iter*m+j+1, tmp);
      }
      if(tmp<tol) { conv = 1; break; }
    }

    //now, H becomes an upper triangule matrix, problem min\|g-Hy\|^2 is g=Hy
    //solve it by backward substitution, y = H(1:m,1:m)\g(1:m)
    mm = MIN(m, j+1);
    y[mm-1] = g[mm-1]/h[mm-1][mm-1];
    for(i=mm-2; i>=0; i--){
      y[i] = g[i];
      for(j=i+1; j<mm; j++) y[i] -= h[j][i]*y[j];
      y[i] /= h[i][i];
    }

    //x=x0+invM*Vm*y
    memset(w, 0, n*sizeof(complex));
    for(j=0; j<mm; j++){
      for(i=0; i<n; i++){
	w[i] += y[j]*v[j][i];
      }
    }
    invMop(n, w, r);
    for(i=0; i<n; i++) x[i] += r[i];
    if(conv) break;
  }//end for iter
  
  free1complex(r);
  free1complex(w);
  free2complex(v);
  free2complex(h);
  free1complex(g);
  free1complex(y);
  free1complex(c);
  free1double(s);
}


//Flexible GMRES
void fgmres(int n, complex *x, complex *b, op_t Aop, op_t invMop, int niter, double tol, int m, int verb)
{
  int i, j, k, mm, iter, conv;
  double beta, tmp, r0, relres;
  complex ss;
  FILE *fp;
  complex *w = alloc1complex(n);
  complex *r = alloc1complex(n);
  complex **v = alloc2complex(n, m+1);
  complex **h = alloc2complex(m+1, m);
  complex **z = alloc2complex(n, m);
  complex *g = alloc1complex(m+1);
  complex *y = alloc1complex(m+1);
  complex *c = alloc1complex(m+1);
  double *s = alloc1double(m+1);

  conv = 0;
  for(iter=0; iter<niter; iter++){
    /*-------------------------------------------------------------*/
    Aop(n, x, w);//w=A*x
    for(i=0; i<n; i++) r[i] = b[i] - w[i];
    beta = sqrt(creal(cdotprod(n, r, r)));
    if(beta==0.) break;
    for(i=0; i<n; i++) v[0][i] = r[i]/beta;
    memset(g, 0, (m+1)*sizeof(complex));
    g[0] = beta;
    memset(&h[0][0], 0, (m+1)*m*sizeof(complex));
    relres = (iter>0)?(beta/r0):1;
    if(iter==0){
      r0 = beta;
      if(verb){
	fp = fopen("iterate_fgmres.txt", "w");
	fprintf(fp, "iteration \t error\n");
	fprintf(fp, "%d \t %e\n", iter, relres);
	fclose(fp);
	printf("f-gmres iter=0 relres=%e\n", relres);
      }
    }
    
    for(j=0; j<m; j++){
      invMop(n, v[j], z[j]);//z=invM*v;
      Aop(n, z[j], w);//w=Az;
      for(i=0; i<=j; i++){
	h[j][i] = cdotprod(n, w, v[i]);
	for(k=0; k<n; k++) w[k] -= h[j][i]*v[i][k];
      }
      h[j][j+1] = sqrt(creal(cdotprod(n, w, w)));
      if(cabs(h[j][j+1])==0.) { conv = 1; break; }
      for(i=0; i<n; i++) v[j+1][i] = w[i]/h[j][j+1];

      //solve least-squares problem by QR factorization using Given rotations
      //min \|g - H(1:m,1:m) y\|^2, g=(beta,0,...,0)
      //min \|G*g - G*Hy\|^2, G=G(i-1,i,theta)*...*G(2,3,theta)*G(1,2,theta)
      if(j>0){
	for(i=0; i<j; i++){
	  //apply G12, G23,..., G_{j-1,j} to the last column of H_{j,*}
	  ss = conj(c[i])*h[j][i] + s[i]*h[j][i+1];
	  h[j][i+1] = -s[i]*h[j][i] + c[i]*h[j][i+1];
	  h[j][i] = ss;
	}
      }

      //compute c=cos(theta) and s=sin(theta)
      tmp = sqrt(creal(h[j][j]*conj(h[j][j]) + h[j][j+1]*conj(h[j][j+1])));
      s[j] = creal(h[j][j+1])/tmp;
      c[j] = h[j][j]/tmp;
      //g=G(j,j+1,theta)g with g[j+1]=0
      g[j+1] = -s[j]*g[j];
      g[j] = conj(c[j])*g[j];
      h[j][j] = conj(c[j])*h[j][j] + s[j]*h[j][j+1];

      relres = cabs(g[j+1])/r0;
      if(verb){
	fp = fopen("iterate_fgmres.txt", "a");
	fprintf(fp, "%d \t %e\n", iter*m+j+1, relres);
	fclose(fp);
	printf("f-gmres iter=%d relres=%e\n", iter*m+j+1, relres);
      }
      if(relres<tol) { conv = 1; break;}
    }
    
    //now, H becomes an upper triangule matrix, problem min\|g-Hy\|^2 is g=Hy
    //solve it by backward substitution, y = H(1:m,1:m)\g(1:m)
    mm = MIN(m, j+1);
    y[mm-1] = g[mm-1]/h[mm-1][mm-1];
    for(i=mm-2; i>=0; i--){
      y[i] = g[i];
      for(j=i+1; j<mm; j++) y[i] -= h[j][i]*y[j];
      y[i] /= h[i][i];
    }

    //x=x0+invMM*Vmm*y
    for(j=0; j<mm; j++){
      for(i=0; i<n; i++){
	x[i] += y[j]*z[j][i];
      }
    }
    if(conv) break;
  }//end for iter
  
  free1complex(r);
  free1complex(w);
  free2complex(v);
  free2complex(z);
  free2complex(h);
  free1complex(g);
  free1complex(y);
  free1complex(c);
  free1double(s);
}
