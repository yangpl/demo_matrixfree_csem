/* GMRES preconditioner
 *------------------------------------------------------------------------
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "solver.h"

//compute dot product s=<x,y>=x^T*conj(y), this definiton is required according to Saad book
complex cdotprod(int n, complex *a, complex *b);

typedef struct{
  int n;
  int m;
  complex *w;
  complex *r;
  complex **v;
  complex **h;
  complex *g;
  complex *y;
  complex *c;
  double *s;
  op_t Aop;
} gmres_t;
gmres_t *gmres_pre;

void gmres_init(int n, int m, op_t Aop)
{
  gmres_pre = (gmres_t*)malloc(sizeof(gmres_t));
  gmres_pre->n = n;
  gmres_pre->m = m;
  gmres_pre->w = alloc1complex(n);
  gmres_pre->r = alloc1complex(n);
  gmres_pre->v = alloc2complex(n, m+1);
  gmres_pre->h = alloc2complex(m, m+1);
  gmres_pre->g = alloc1complex(m+1);
  gmres_pre->y = alloc1complex(m+1);
  gmres_pre->c = alloc1complex(m+1);
  gmres_pre->s = alloc1double(m+1);
  gmres_pre->Aop = Aop;
}

void gmres_free()
{
  free1complex(gmres_pre->r);
  free1complex(gmres_pre->w);
  free2complex(gmres_pre->v);
  free2complex(gmres_pre->h);
  free1complex(gmres_pre->g);
  free1complex(gmres_pre->y);
  free1complex(gmres_pre->c);
  free1double(gmres_pre->s);
  free(gmres_pre);
}

void gmres_apply(int n, complex *b, complex *x)
{
  int i, j, k;
  double beta, tmp;
  complex ss;
  op_t Aop = gmres_pre->Aop;
  int m = gmres_pre->m;
  complex *w = gmres_pre->w;
  complex *r = gmres_pre->r;
  complex **v = gmres_pre->v;
  complex **h = gmres_pre->h;
  complex *g = gmres_pre->g;
  complex *y = gmres_pre->y;
  complex *c = gmres_pre->c;
  double *s = gmres_pre->s;
  
  //Aop(n, x, w);//w=A*x
  memset(x, 0, n*sizeof(complex));
  for(i=0; i<n; i++) r[i] = b[i];// -w[i];
  beta = sqrt(cdotprod(n, r, r));
  if(beta==0.0) return;
  for(i=0; i<n; i++) v[0][i] = r[i]/beta;
  memset(g, 0, (m+1)*sizeof(complex));
  g[0] = beta;
  memset(&h[0][0], 0, (m+1)*m*sizeof(complex));
    
  for(j=0; j<m; j++){
    Aop(n, v[j], w);//r=Av;
    for(i=0; i<=j; i++) {
      h[i][j] = cdotprod(n, w, v[i]);
      for(k=0; k<n; k++) w[k] -= h[i][j]*v[i][k];
    }
    h[j+1][j] = sqrt(creal(cdotprod(n, w, w)));

    //if(cabs(h[j+1][j])==0.0) { m=j+1; break; }
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

    //tmp = cabs(g[j+1])/r0;
    //if(tmp<tol) { m=j+1; return; }
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
  
}






