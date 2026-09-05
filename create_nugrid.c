/* generate nonuniform grid using geometric progression 
 * determine the optimal ratio q by root finding using fixed point iteration 
 * Reference: 
 *  Controlled-source electromagnetic modeling using a high-order finite-difference 
 *  time-domain method on a nonuniform grid. P Yang, R Mittet Geophysics 88 (2), E53-E67
 *--------------------------------------------------------------------
 *
 *   Copyright (c) 2020-2025, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *--------------------------------------------------------------------*/
#include <math.h>

/*< create nonuniform grid using fixed point iteration >*/
double create_nugrid(int n, double len, double dx, double *x, int istretch)
{
  int i, iter;
  double q, qq;
  double eps = 1e-15;
  int niter = 100;

  if(n*dx>=len||!istretch) {
    for(i=0; i<=n; i++) x[i] = i*dx;
    return 1;
  }
  
  q = 1.1;
  qq = 1;
  for(iter=0; iter<niter; ++iter){
    qq = pow(len*(q-1.)/dx + 1., 1./n);
    if(fabs(qq-q)<eps) break;
    q = qq;
  }

  for(x[0]=0,i=1; i<=n; i++)
    x[i] = (pow(q,i) - 1.)*dx/(q-1.);

  return q;
}
