#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* ========================================================
   make a 3D conductivity model including multiple interfaces
   (layer cake: air, water, sediments, basement, etc.)
   --------------------------------------------------------
   Inputs:
   nx, ny, nz  : grid sizes
   xs[nx], ys[ny], zs[nz] : fine grid cell centers
   nlayer      : number of materials
   zsurf[nlayer-1][ny*nx] : depth surfaces between layers
   (layer ll is between zsurf[ll-1] and zsurf[ll])
   Must be monotone increasing in depth.
   sigma_layer[nlayer][3]   : conductivity for each layer
   if anisotropic=0: only use [ll][0]
   if anisotropic=1: use [ll][0]=σxx, [ll][1]=σyy, [ll][2]=σzz
   anisotropic : 0 or 1

   Outputs:
   isotropic: sigma[nz*ny*nx] (NULL if anisotropic=1)
   anisotropic: sigma_xx, sigma_yy, sigma_zz (NULL if anisotropic=0)
   ======================================================== */
void make_model(int nx, int ny, int nz, const float *xs, const float *ys, const float *zs,
		int nlayer,
		const float **zsurf,   /* array of (nlayer-1) pointers, each [ny*nx] */
		const float rho_layer[][3],
		int anisotropic,/* 0 or 1 */
		float *rho,/* NULL if anisotropic=1 */
		float *rho_xx, float *rho_yy, float *rho_zz /*NULL if anisotropic=0*/)
{
  int i,j,k,ll;
  float *dzs   = (float*) malloc(nz*sizeof(float));
  float *z_top = (float*) malloc(nz*sizeof(float));
  float *z_bot = (float*) malloc(nz*sizeof(float));

  /* voxel thickness */
  for (k=0;k<nz;k++) {
    if (k==0) dzs[k] = zs[1]-zs[0];
    else if (k==nz-1) dzs[k] = zs[nz-1]-zs[nz-2];
    else dzs[k] = 0.5*((zs[k+1]-zs[k])+(zs[k]-zs[k-1]));
    z_top[k] = zs[k]-0.5*dzs[k];
    z_bot[k] = zs[k]+0.5*dzs[k];
  }

  for (j=0;j<ny;j++) {
    for (i=0;i<nx;i++) {
      /* collect vertical interface depths at this (x,y) */
      float zcut[32]; /* support up to 32 layers */
      if (nlayer > 32) { fprintf(stderr,"Too many layers!\n"); exit(1); }
      zcut[0] = -1e30; /* top open */
      for (ll=1;ll<nlayer;ll++) {
	zcut[ll] = zsurf[ll-1][j*nx+i];
      }
      zcut[nlayer] = 1e30; /* bottom open */

      for (k=0;k<nz;k++) {
	float sum_frac = 0.0;
	float effx=0, effy=0, effz=0;
	int idx=(k*ny+j)*nx+i;

	/* voxel bounds */
	float zt=z_top[k], zb=z_bot[k], dz=dzs[k];

	for (ll=0;ll<nlayer;ll++) {
	  /* overlap of voxel with this layer */
	  float zlow=fmax(zt,zcut[ll]);
	  float zhigh=fmin(zb,zcut[ll+1]);
	  float thick=zhigh-zlow;
	  if (thick<=0) continue;
	  float frac=thick/dz;
	  sum_frac+=frac;

	  if (anisotropic) {
	    /*arithmetic average over sigma_xx and sigma_yy*/
	    effx += frac/rho_layer[ll][0];
	    effy += frac/rho_layer[ll][1];
	    /* harmonic average over sigma_zz */
	    effz += frac*rho_layer[ll][2];
	  } else {
	    effx += frac/rho_layer[ll][0];
	  }
	}

	if (anisotropic) {
	  rho_xx[idx]=1./effx;
	  rho_yy[idx]=1./effy;
	  rho_zz[idx]=effz;
	} else {
	  rho[idx]=1./effx;
	}
	
	if (fabs(sum_frac-1.0)>1e-4) fprintf(stderr,"Warning: voxel fractions sum=%.3f at (i=%d,j=%d,k=%d)\n", sum_frac,i,j,k);
      }
    }
  }

  free(dzs);
  free(z_top);
  free(z_bot);
}


/* generate nonuniform grid using geometric progression */
/* determine the optimal ratio r by root finding using fixed point iteration */
float create_nugrid(int n, float len, float dx)
{
  int i;
  float r, rr;
  float eps = 1e-15;

  if(fabs(n*dx-len)<eps) {
    return 1;
  }
  
  //fixed point iterations to find the root-"r" of the equation len = dx *(r^n-1)/(r-1)
  // assume n intervals (n+1 points/nodes)
  r = 1.1;
  rr = 1;
  while(1){
    rr = pow(len*(r-1.)/dx + 1., 1./n);
    if(fabs(rr-r)<eps) break;
    r = rr;
  }
  return r;
}


void create_dense_model(int nx, int ny, int nz,
			float xmin, float ymin, float zmin,
			float xmax, float ymax, float zmax)
{
  FILE *fp;

  float dx = (xmax-xmin)/nx;
  float dy = (ymax-ymin)/ny;
  float dz = (zmax-zmin)/nz;

  float *x1 = malloc((nx+1)*sizeof(float));
  float *x2 = malloc((ny+1)*sizeof(float));
  float *x3 = malloc((nz+1)*sizeof(float));
  for(int i=0; i<=nx; i++) x1[i] = xmin + i*dx;
  for(int j=0; j<=ny; j++) x2[j] = ymin + j*dy;
  for(int k=0; k<=nz; k++) x3[k] = zmin + k*dz;

  //generate FD cell centers
  float *xs = malloc(nx*sizeof(float));
  float *ys = malloc(ny*sizeof(float));
  float *zs = malloc(nz*sizeof(float));
  for(int i=0; i<nx; i++) xs[i] = 0.5*(x1[i] + x1[i+1]);
  for(int j=0; j<ny; j++) ys[j] = 0.5*(x2[j] + x2[j+1]);
  for(int k=0; k<nz; k++) zs[k] = 0.5*(x3[k] + x3[k+1]);
  
  //=======================================================
  /* 4 layers: air, water, formation-1, formation-2 + resistor */
  int nlayer=4;
  float PI = 3.14159265;
  float Lx = 9000 - (-9000);
  float eta = PI*4/Lx;
  /* resistivities: air, seawater, rock */
  float rho_layer[4][3] = {
    {1e8, 1e8, 1e8},/* air */
    {0.3, 0.3, 0.3},/* seawater */
    {1.0, 1.0, 1.0},/* rock/formation-1 */
    {2.0, 2.0, 2.0} /* rock/formation-2 */
  };

  /* two interfaces: air-water, water-rock */
  float *surf1 = malloc(nx*ny*sizeof(float));
  float *surf2 = malloc(nx*ny*sizeof(float));
  float *surf3 = malloc(nx*ny*sizeof(float));
  for (int j=0;j<ny;j++){
    for (int i=0;i<nx;i++){
      surf1[j*nx+i]=0;/* sea surface z1(x,y) */
      surf2[j*nx+i]=1000. + 100.*sin(eta*xs[i]);/* seabed z2(x,y) */
      surf3[j*nx+i]=2500.;      
    }
  }
  const float *zsurf[3]={surf1, surf2, surf3};
  
  fp = fopen("ftopo.txt", "w");
  for (int i=0;i<nx;i++) fprintf(fp, "%e \t %e\n", xs[i], surf2[i]);
  fclose(fp);

  float *rho=(float*)malloc(nx*ny*nz*sizeof(float));
  make_model(nx, ny, nz, xs, ys, zs, nlayer, zsurf, rho_layer, 0, rho, NULL, NULL, NULL);  /* 1=anisotropic */
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
	if(xs[i]>=-5000 && xs[i]<5000 && zs[k]>=2000 && zs[k]<2100){
	  int idx = i + nx*(j + ny*k);
	  rho[idx] = 100.;//add a resistor
	}
      }
    }
  }
  //output resistivity files
  fp = fopen("frho", "wb");
  fwrite(rho, nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

  //=====================================================
  int nr = 81;
  float dr = Lx/(nr - 1);//receiver spacing

  fp = fopen("receivers.txt", "w");
  fprintf(fp, "x \t y \t z \t azimuth \t dip \t irec\n");
  for(int i=0; i<nr; i++){
    float xi = -9000 + i*dr;
    float zi = 1000. + 100.*sin(eta*xi);//3 periods of sine shape
    fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %d\n", xi, 0., zi, 0., 0., i+1);
  }
  fclose(fp);

  fp = fopen("sources.txt", "w");
  fprintf(fp, "x \t y \t z \t azimuth \t dip \t isrc\n");
  fprintf(fp, "%f \t %f \t %f \t %f \t %f \t %d\n", 0., 0., 830.0, 0., 0., 1);
  fclose(fp);

  fp = fopen("src_rec_table.txt", "w");
  fprintf(fp, "isrc \t irec\n");
  for(int i=0; i<nr; i++) fprintf(fp, "%d \t %d\n", 1, i+1);
  fclose(fp);
  
  free(rho);
  free(surf1);
  free(surf2);
  free(surf3);
  free(x1);
  free(x2);
  free(x3);
  free(xs);
  free(ys);
  free(zs);
}

int main(int argc, char *argv[])
{
  FILE *fp;
  int i, j, k;
  int ix, iy, iz;
  int nx = 2000;
  int ny = 1000;
  int nz = 1000;
  float xmin = -10e3;
  float ymin = -10e3;
  float zmin = 0;
  float xmax = 10e3;
  float ymax = 10e3;
  float zmax = 4e3;

  float *xx = malloc((nx+1)*sizeof(float));
  float *yy = malloc((ny+1)*sizeof(float));
  float *zz = malloc((nz+1)*sizeof(float));
  float dx = (xmax-xmin)/nx;
  float dy = (ymax-ymin)/ny;
  float dz = (zmax-zmin)/nz;
  for(i=0; i<=nx; i++) xx[i] = xmin + i*dx;
  for(j=0; j<=ny; j++) yy[j] = ymin + j*dy;
  for(k=0; k<=nz; k++) zz[k] = zmin + k*dz;

  create_dense_model(nx, ny, nz, xmin, ymin, zmin, xmax, ymax, zmax);
  printf("dense model created!\n");
  
  float *rho = malloc(nx*ny*nz*sizeof(float));
  fp = fopen("frho", "rb");
  fread(rho, nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

  int n1 = 125;
  int n2 = 125;
  int n3 = 100;
  float *x1 = malloc((n1+1)*sizeof(float));
  float *x2 = malloc((n2+1)*sizeof(float));
  float *x3 = malloc((n3+1)*sizeof(float));

  float d1 = (xmax-xmin)/(float)(n1);
  for(i=0; i<=n1; i++) x1[i] = xmin + i*d1;
  float d2 = (ymax-ymin)/(float)(n2);
  for(j=0; j<=n2; j++) x2[j] = ymin + j*d2;
  /*
  float d3 = 25;
  float r = create_nugrid(n3, zmax-zmin, d3);
  x3[0] = 0;
  for(k=1; k<=n3; k++) x3[k] = x3[k-1] + d3*pow(r, k-1);
  */
  /* float d3 = (zmax-zmin)/(float)(n3); */
  /* for(k=0; k<=n3; k++) x3[k] = zmin + d3*k; */
  
  float d3 = 20;
  for(k=0; k<=n3; k++) x3[k] = zmin + k*d3;
  
  iz = 65;
  float zstart = iz*d3;//nonuniform grid at depth>=zstart
  int nh = n3-iz;
  float r = create_nugrid(nh, zmax-zstart, d3);
  for(k=iz+1; k<=n3; k++){
    i = k-iz;
    x3[k] = x3[k-1] + d3*pow(r, i-1);
  }

  fp = fopen("fx1", "wb");
  fwrite(x1, (n1+1)*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("fx2", "wb");
  fwrite(x2, (n2+1)*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("fx3", "wb");
  fwrite(x3, (n3+1)*sizeof(float), 1, fp);
  fclose(fp);

  float a, b, s;
  float left, right, overlap;
  float *rho1, *rho2, *rho3;

  printf("average sigma=1./rho\n");
  //======================================================
  rho1 = malloc(n1*ny*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(iy=0; iy<ny; iy++){

      i = 0;
      for(j=0; j<n1; j++) {
	a = x1[j];
	b = x1[j+1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<nx && xx[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nx && xx[k]<b) {
	  overlap = 0.0;
	  left = (xx[k]>a)? xx[k]:a;
	  right = (xx[k+1]<b)? xx[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += 1./rho[k + nx*(iy + ny*iz)]*overlap;
	  k++;
	}
	rho1[j + n1*(iy + ny*iz)] = s/(b-a);
      }
      
    }//end for iy
  }//end for iz
  printf("x done!\n");
  
  rho2 = malloc(n1*n2*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n2; j++) {
	a = x2[j];
	b = x2[j + 1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<ny && yy[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<ny && yy[k]<b) {
	  overlap = 0.0;
	  left = (yy[k]>a)? yy[k]:a;
	  right = (yy[k+1]<b)? yy[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho1[ix + n1*(k + ny*iz)]*overlap;
	  k++;
	}
	rho2[ix + n1*(j + n2*iz)] = s/(b-a);
      }
      
    }//end for ix
  }//end for iz
  free(rho1);
  printf("y done!\n");

  rho3 = malloc(n1*n2*n3*sizeof(float));
  for(iy=0; iy<n2; iy++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n3; j++) {
	a = x3[j];
	b = x3[j + 1];
	s = 0.0;
	/* advance i to skip source cells entirely left of this target cell */
	while(i<nz && zz[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nz && zz[k]<b) {
	  overlap = 0.0;
	  left = (zz[k]>a)? zz[k]:a;
	  right = (zz[k+1]<b)? zz[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho2[ix + n1*(iy + n2*k)]*overlap;
	  k++;
	}
	rho3[ix + n1*(iy + n2*j)] = s/(b-a);
	rho3[ix + n1*(iy + n2*j)] = 1./rho3[ix + n1*(iy + n2*j)];//convert averaged 1./rho to rho
      }
      
    }//end for ix
  }//end for iy
  free(rho2);
  printf("z done!\n");
  
  fp = fopen("frho11", "wb");
  fwrite(rho3, n1*n2*n3*sizeof(float), 1, fp);
  fclose(fp);

  fp = fopen("frho22", "wb");
  fwrite(rho3, n1*n2*n3*sizeof(float), 1, fp);
  fclose(fp);

  free(rho3);

  printf("average rho\n");
  //======================================================
  rho1 = malloc(n1*ny*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(iy=0; iy<ny; iy++){

      i = 0;
      for(j=0; j<n1; j++) {
	a = x1[j];
	b = x1[j + 1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<nx && xx[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nx && xx[k]<b) {
	  overlap = 0.0;
	  left = (xx[k]>a)? xx[k]:a;
	  right = (xx[k+1]<b)? xx[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho[k + nx*(iy + ny*iz)]*overlap;
	  k++;
	}
	rho1[j + n1*(iy + ny*iz)] = s/(b-a);
      }
      
    }//end for iy
  }//end for iz
  printf("x done!\n");
  
  rho2 = malloc(n1*n2*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n2; j++) {
	a = x2[j];
	b = x2[j + 1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<ny && yy[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<ny && yy[k]<b) {
	  overlap = 0.0;
	  left = (yy[k]>a)? yy[k]:a;
	  right = (yy[k+1]<b)? yy[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho1[ix + n1*(k + ny*iz)]*overlap;
	  k++;
	}
	rho2[ix + n1*(j + n2*iz)] = s/(b-a);
      }
      
    }//end for ix
  }//end for iz
  free(rho1);
  printf("y done!\n");

  rho3 = malloc(n1*n2*n3*sizeof(float));
  for(iy=0; iy<n2; iy++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n3; j++) {
	a = x3[j];
	b = x3[j + 1];
	s = 0.0;
	/* advance i to skip source cells entirely left of this target cell */
	while(i<nz && zz[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nz && zz[k]<b) {
	  overlap = 0.0;
	  left = (zz[k]>a)? zz[k]:a;
	  right = (zz[k+1]<b)? zz[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho2[ix + n1*(iy + n2*k)]*overlap;
	  k++;
	}
	rho3[ix + n1*(iy + n2*j)] = s/(b-a);
      }
      
    }//end for ix
  }//end for iy
  free(rho2);
  printf("z done!\n");
  
  fp = fopen("frho33", "wb");
  fwrite(rho3, n1*n2*n3*sizeof(float), 1, fp);
  fclose(fp);

  free(rho3);
  
  free(rho);
  free(xx);
  free(yy);
  free(zz);
}


