/* initialize and free arrays allocated by electromagnetic field pointer
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

int cmpfunc(const void *a, const void *b) { return ( *(int*)a - *(int*)b ); }

void emf_init(emf_t *emf)
{
  int i1, i2, i3, ifreq, ic;
  int istat, n;
  char *fx1, *fx2, *fx3;
  char *frho11, *frho22, *frho33;
  long filesize;
  FILE *fp;

  //------- read frequencies ------------------
  if(!(emf->nfreq=countparval("freqs"))) err("Need freqs= vector");
  emf->freqs=alloc1float(emf->nfreq);
  getparfloat("freqs", emf->freqs);/* a list of frequencies separated by comma */
  qsort(emf->freqs, emf->nfreq, sizeof(float), cmpfunc);/*sort frequencies in ascending order*/
  //-------- read active source channels ------
  if((emf->nchsrc=countparval("chsrc"))!=0) {
    emf->chsrc=(char**)alloc1(emf->nchsrc, sizeof(void*));
    getparstringarray("chsrc", emf->chsrc);
    /* active source channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchsrc=1;
    emf->chsrc=(char**)alloc1(emf->nchsrc, sizeof(void*));
    emf->chsrc[0]="Ex";
  }
  //-------- read active receiver channels ------
  if((emf->nchrec=countparval("chrec"))!=0) {
    emf->chrec=(char**)alloc1(emf->nchrec, sizeof(void*));
    getparstringarray("chrec", emf->chrec);
    /* active receiver channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchrec=1;
    emf->chrec=(char**)alloc1(emf->nchrec, sizeof(void*));
    emf->chrec[0] = "Ex";
  }
  if(emf->verb){
    printf("----------- modelling parameters -------\n");
    printf("freqs=");
    for(ifreq=0; ifreq<emf->nfreq; ++ifreq) printf("%g,", emf->freqs[ifreq]);
    printf("\n");
    printf("Active source channels=");
    for(ic=0; ic<emf->nchsrc; ++ic) printf("%s,", emf->chsrc[ic]);
    printf("\n");
    printf("Active recever channels=");
    for(ic=0; ic<emf->nchrec; ++ic) printf("%s,", emf->chrec[ic]);
    printf("\n");
  }

  if(!getparint("reciprocity", &emf->reciprocity)) emf->reciprocity = 0;/*1=switch source receiver locations */
  if(!getparint("fvm", &emf->fvm)) emf->fvm = 0;//1=FVM; 0=FDM
  if(!getparint("nb", &emf->nb)) emf->nb = 15; /* number of layers padded as absorbing boundaries */
  if(!getparfloat("rho_air", &emf->rho_air)) emf->rho_air = 1e6;//air resistivity
  if(!getparint("istretch", &emf->istretch)) emf->istretch = 1; /* 1=grid stretch; 0=no stretching */
  if(!getparint("airbc", &emf->airbc)) emf->airbc = 0;//1=air-water BC without air; 0=air included
  if(emf->airbc && emf->nb==0) err("airbc=1 assumes input model without air, you must ensure nb>0");
  if(emf->verb){
    printf("reciprocity=%d\n", emf->reciprocity);
    printf("fvm=%d\n", emf->fvm);
    printf("airbc=%d (1=use air BC; 0=including air)\n", emf->airbc);
    printf("istretch=%d\n", emf->istretch);
    printf("rho_air=%g \n", emf->rho_air);
    printf("nb=%d \n", emf->nb);
  }
  //------- read coordinates x1, x2, x3 ------------
  if(!(getparstring("fx1", &fx1))) err("Need fx1= ");
  if(!(getparstring("fx2", &fx2))) err("Need fx2= ");
  if(!(getparstring("fx3", &fx3))) err("Need fx3= ");

  fp = fopen(fx1, "rb");
  if(fp==NULL) err("cannot open file fx1=%s\n", fx1);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  emf->nx = filesize/sizeof(float)-1;
  emf->x1node = alloc1float(emf->nx+1);
  rewind(fp);
  istat = fread(emf->x1node, (emf->nx+1)*sizeof(float), 1, fp);
  fclose(fp);
  
  fp = fopen(fx2, "rb");
  if(fp==NULL) err("cannot open file fx2=%s\n", fx2);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  emf->ny = filesize/sizeof(float)-1;
  emf->x2node = alloc1float(emf->ny+1);
  rewind(fp);
  istat = fread(emf->x2node, (emf->ny+1)*sizeof(float), 1, fp);
  fclose(fp);
  
  fp = fopen(fx3, "rb");
  if(fp==NULL) err("cannot open file fx3=%s\n", fx3);
  fseek(fp, 0, SEEK_END);
  filesize = ftell(fp);
  emf->nz = filesize/sizeof(float)-1;
  emf->x3node = alloc1float(emf->nz+1);
  rewind(fp);
  istat = fread(emf->x3node, (emf->nz+1)*sizeof(float), 1, fp);
  fclose(fp);
  
  if(emf->verb){
    printf("----------- input model -----------\n");
    printf("nx=%d\n", emf->nx);
    printf("ny=%d\n", emf->ny);
    printf("nz=%d\n", emf->nz);
    printf("model domain [x1min, x1max]=[%g, %g]\n", emf->x1node[0], emf->x1node[emf->nx]);
    printf("model domain [x2min, x2max]=[%g, %g]\n", emf->x2node[0], emf->x2node[emf->ny]);
    printf("model domain [x3min, x3max]=[%g, %g]\n", emf->x3node[0], emf->x3node[emf->nz]);
    printf("input Lx=%g\n", emf->x1node[emf->nx]-emf->x1node[0]);
    printf("input Ly=%g\n", emf->x2node[emf->ny]-emf->x2node[0]);
    printf("input Lz=%g\n", emf->x3node[emf->nz]-emf->x3node[0]);
  }

  //---------- read model rho11, rho22, rho33 ------------------
  emf->rho11 = alloc3float(emf->nx, emf->ny, emf->nz);
  emf->rho22 = alloc3float(emf->nx, emf->ny, emf->nz);
  emf->rho33 = alloc3float(emf->nx, emf->ny, emf->nz);
  n = emf->nx*emf->ny*emf->nz;//compute total number of FD cells

  if(!(getparstring("frho11", &frho11))) err("Need frho11= ");
  if(!(getparstring("frho22", &frho22))) err("Need frho22= ");
  if(!(getparstring("frho33", &frho33))) err("Need frho33= ");

  fp = fopen(frho11, "rb");
  if(fp==NULL) err("cannot open file frho11=%s\n", frho11);
  istat = fread(&emf->rho11[0][0][0], sizeof(float), n, fp);
  if(istat != n) err("size not match frho11: file=%d (nx+1)*(ny+1)*(nz+1)=%d\n", istat, n);
  fclose(fp);

  fp = fopen(frho22, "rb");
  if(fp==NULL) err("cannot open file frho22=%s\n", frho22);
  istat = fread(&emf->rho22[0][0][0], sizeof(float), n, fp);
  if(istat != n) err("size not match frho22: file=%d (nx+1)*(ny+1)*(nz+1)=%d\n", istat, n);
  fclose(fp);

  fp = fopen(frho33, "rb");
  if(fp==NULL) err("cannot open file frho33=%s\n", frho33);
  istat = fread(&emf->rho33[0][0][0], sizeof(float), n, fp);
  if(istat != n) err("size not match frho33: file=%d (nx+1)*(ny+1)*(nz+1)=%d\n", istat, n);
  fclose(fp);

  emf->rhomax = emf->rho11[0][0][0];
  emf->rhomin = emf->rho11[0][0][0];
  for(i3=0; i3<emf->nz; i3++){
    for(i2=0; i2<emf->ny; i2++){
      for(i1=0; i1<emf->nx; i1++){
	emf->rhomax = MAX(emf->rhomax, emf->rho11[i3][i2][i1]);
	emf->rhomin = MIN(emf->rhomin, emf->rho11[i3][i2][i1]);
	emf->rhomax = MAX(emf->rhomax, emf->rho22[i3][i2][i1]);
	emf->rhomin = MIN(emf->rhomin, emf->rho22[i3][i2][i1]);
	emf->rhomax = MAX(emf->rhomax, emf->rho33[i3][i2][i1]);
	emf->rhomin = MIN(emf->rhomin, emf->rho33[i3][i2][i1]);
      }
    }
  }
  emf->rho_water = emf->rhomin;//water resistivity, assume input model without air
  if(emf->verb){
    printf("[rhomin,rhomax]=[%g,%g]\n", emf->rhomin, emf->rhomax);
    printf("rho_water=%g\n", emf->rho_water);
  }

}


void emf_free(emf_t *emf)
{
  free1float(emf->freqs);
  free1float(emf->x1node);
  free1float(emf->x2node);
  free1float(emf->x3node);
  free3float(emf->rho11);
  free3float(emf->rho22);
  free3float(emf->rho33);
}

