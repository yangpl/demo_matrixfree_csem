#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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


int main(int argc, char *argv[])
{
  long filesize;
  FILE *fp;
  int nx, ny, nz;
  float *xx, *yy, *zz;
  int i, j, k;
  int ix, iy, iz;
  int istat;
  int na, nb;
  float r;
  
  //------- read coordinates x1, x2, x3 ------------
  fp = fopen("x1", "rb");
  if(fp==NULL) { printf("cannot open file x1\n"); exit(1); }
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  nx = filesize/sizeof(float)-1;
  xx = malloc((nx+1)*sizeof(float));
  rewind(fp);
  istat = fread(xx, (nx+1)*sizeof(float), 1, fp);
  fclose(fp);
  
  fp = fopen("x2", "rb");
  if(fp==NULL) { printf("cannot open file x2\n"); exit(1); }
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  ny = filesize/sizeof(float)-1;
  yy = malloc((ny+1)*sizeof(float));
  rewind(fp);
  istat = fread(yy, (ny+1)*sizeof(float), 1, fp);
  fclose(fp);
  
  fp = fopen("x3", "rb");
  if(fp==NULL) { printf("cannot open file x3\n"); exit(1); }
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  nz = filesize/sizeof(float)-1;
  zz = malloc((nz+1)*sizeof(float));
  rewind(fp);
  istat = fread(zz, (nz+1)*sizeof(float), 1, fp);
  fclose(fp);

  float xmin = xx[0];
  float xmax = xx[nx];
  float ymin = yy[0];
  float ymax = yy[ny];
  float zmin = zz[0];
  float zmax = zz[nz];
  printf("[xmin,xmax]=[%g,%g], Lx=%g\n", xmin, xmax, xmax-xmin);
  printf("[ymin,ymax]=[%g,%g], Ly=%g\n", ymin, ymax, ymax-ymin);
  printf("[zmin,zmax]=[%g,%g], Lz=%g\n", zmin, zmax, zmax-zmin);

  int n1 = 180;
  int n2 = 120;
  int n3 = 120;
  float *x1 = malloc((n1+1)*sizeof(float));
  float *x2 = malloc((n2+1)*sizeof(float));
  float *x3 = malloc((n3+1)*sizeof(float));
  
  float xs = 390275.0;
  float ys = 7517812.0;
  float zs = 849.783813;

  float d1 = 100;
  na = n1/2;
  nb = n1 - na;
  x1[na] = xs;
  r = create_nugrid(na, xs-xmin, d1);
  printf("left strech r=%g\n", r);
  for(i=0; i<na; i++) x1[na-i-1] = x1[na-i] - d1*pow(r,i);
  r = create_nugrid(nb, xmax-xs, d1);
  printf("right strech r=%g\n", r);
  for(i=0; i<nb; i++) x1[na+i+1] = x1[na+i] + d1*pow(r,i);
  
  float d2 = 100;
  na = n2/2;
  nb = n2 - na;
  x2[na] = ys;
  r = create_nugrid(na, ys-ymin, d2);
  printf("front strech r=%g\n", r);
  for(i=0; i<na; i++) x2[na-i-1] = x2[na-i] - d2*pow(r,i);
  r = create_nugrid(nb, ymax-ys, d2);
  printf("rear strech r=%g\n", r);
  for(i=0; i<nb; i++) x2[na+i+1] = x2[na+i] + d2*pow(r,i);
  /*
  float d1 = (xmax-xmin)/(float)(n1);
  for(i=0; i<=n1; i++) x1[i] = xmin + i*d1;
  float d2 = (ymax-ymin)/(float)(n2);
  for(j=0; j<=n2; j++) x2[j] = ymin + j*d2;
  printf("d1=%g\n", d1);
  printf("d2=%g\n", d2);
  */

  float d3 = 30;
  for(k=0; k<=n3; k++) x3[k] = zmin + k*d3;
  iz = 60;
  float zstart = iz*d3;//nonuniform grid at depth>=zstart
  int nh = n3-iz;
  r = create_nugrid(nh, zmax-zstart, d3);
  printf("bottom strech r=%g\n", r);
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
  float *rho;
  rho = malloc(nx*ny*nz*sizeof(float));

  //======================================================
  printf("average sigma_h=1./rho_h\n");
  fp = fopen("Rh", "rb");
  fread(rho, nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);

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

  //======================================================
  printf("average rho_v\n");
  fp = fopen("Rv", "rb");
  fread(rho, nx*ny*nz*sizeof(float), 1, fp);
  fclose(fp);
  
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


