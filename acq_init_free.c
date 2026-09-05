/* read acquisition file for source and receiver geometry
 *--------------------------------------------------------------------------
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include <mpi.h>
#include "cstd.h"
#include "acq.h"
#include "emf.h"
 
void acq_init(acq_t *acq, emf_t * emf)
/*< read acquisition file to initialize survey geometry >*/
{
  static int nd = 5000;//maximum dimensions for the number of source and receiver
  float src_x1[nd], src_x2[nd], src_x3[nd], src_hd[nd], src_pit[nd];/* source receiver coordinates */
  float rec_x1[nd], rec_x2[nd], rec_x3[nd], rec_hd[nd], rec_pit[nd];/* source receiver coordinates */
  float x, y, z, hd, pit;
  int isrc, irec, iseof, idx;
  char *fsrc, *frec;
  FILE *fp;

  if(iproc==0) printf("---------- read acquisition ------------\n");
  if(!getparfloat("lensrc", &acq->lensrc)) acq->lensrc = 278;
  if(!getparfloat("lenrec", &acq->lenrec)) acq->lenrec = 8;
  if(!getparint("nsubsrc", &acq->nsubsrc)) acq->nsubsrc = 1;
  if(!getparint("nsubrec", &acq->nsubrec)) acq->nsubrec = 1;
  if(!(getparstring("fsrc", &fsrc))) err("Need fsrc= ");
  /* file to specify all possible source locations */
  if(!(getparstring("frec", &frec))) err("Need frec= ");
  /* file to specify all possible receiver locations */  

  /*============================================*/
  /* step 2: read all possible source locations */
  /*============================================*/
  fp = fopen(fsrc,"r");
  if(fp==NULL) err("file does not exist!"); 
  iseof = fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  isrc = 0;
  while(1){
    /* (northing,easting,depth)=(y,x,z);   azimuth = heading;  dip=pitch */
    iseof=fscanf(fp,"%f %f %f %f %f %d", &x, &y, &z, &hd, &pit, &idx);
    if(iseof==EOF)
      break;
    else{
      src_x1[isrc] = x;
      src_x2[isrc] = y;
      src_x3[isrc] = z;
      src_hd[isrc] = hd;
      src_pit[isrc] = pit;
      isrc++;
    }
  }
  acq->nsrc = 1; //isrc;
  fclose(fp);
    
  /*==============================================*/
  /* step 1: read all possible receiver locations */
  /*==============================================*/
  fp = fopen(frec,"r");
  if(fp==NULL) err("file does not exist!"); 
  iseof = fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  irec = 0;
  while(1){
    /* (northing,easting,depth)=(y,x,z);   azimuth = heading;  dip=pitch */
    iseof=fscanf(fp,"%f %f %f %f %f %d", &x, &y, &z, &hd, &pit, &idx);
    if(iseof==EOF)
      break;
    else{
      rec_x1[irec] = x;
      rec_x2[irec] = y;
      rec_x3[irec] = z;
      rec_hd[irec] = hd;
      rec_pit[irec] = pit;
      irec++;
    }
  }
  acq->nrec = irec;
  fclose(fp);

  //-----------------------------------------------------------
  acq->x1min = src_x1[0];
  acq->x1max = src_x1[0];
  acq->x2min = src_x2[0];
  acq->x2max = src_x2[0];
  acq->x3min = src_x3[0];
  acq->x3max = src_x3[0];
  acq->nsrc = 1; /* assume 1 source per process by default */
  acq->src_x1 = alloc1float(acq->nsrc);
  acq->src_x2 = alloc1float(acq->nsrc);
  acq->src_x3 = alloc1float(acq->nsrc);
  acq->src_azimuth = alloc1float(acq->nsrc);
  acq->src_dip = alloc1float(acq->nsrc);
  for(isrc=0; isrc<acq->nsrc; ++isrc){
    acq->src_x1[isrc] = src_x1[isrc];
    acq->src_x2[isrc] = src_x2[isrc];
    acq->src_x3[isrc] = src_x3[isrc];
    acq->src_azimuth[isrc] = src_hd[isrc];
    acq->src_dip[isrc] = src_pit[isrc];
          
    acq->x1min = MIN(acq->x1min, acq->src_x1[isrc]);
    acq->x1max = MAX(acq->x1max, acq->src_x1[isrc]);
    acq->x2min = MIN(acq->x2min, acq->src_x2[isrc]);
    acq->x2max = MAX(acq->x2max, acq->src_x2[isrc]);
    acq->x3min = MIN(acq->x3min, acq->src_x3[isrc]);
    acq->x3max = MAX(acq->x3max, acq->src_x3[isrc]);
  }/* end for isrc */

  //-----------------------------------------------------------
  acq->rec_x1 = alloc1float(acq->nrec);
  acq->rec_x2 = alloc1float(acq->nrec);
  acq->rec_x3 = alloc1float(acq->nrec);
  acq->rec_azimuth = alloc1float(acq->nrec);
  acq->rec_dip = alloc1float(acq->nrec);
  for(irec=0; irec<acq->nrec; ++irec){//we always have: acq->nrec <= acq->nrec_total
    //nrec < nrec_total if only inline data are used
    //idx=index of the receivers associated with current source or common receiver gather
    acq->rec_x1[irec] = rec_x1[irec];
    acq->rec_x2[irec] = rec_x2[irec];
    acq->rec_x3[irec] = rec_x3[irec];
    acq->rec_azimuth[irec] = rec_hd[irec];
    acq->rec_dip[irec] = rec_pit[irec];
          
    acq->x1min = MIN(acq->x1min, acq->rec_x1[irec]);
    acq->x1max = MAX(acq->x1max, acq->rec_x1[irec]);
    acq->x2min = MIN(acq->x2min, acq->rec_x2[irec]);
    acq->x2max = MAX(acq->x2max, acq->rec_x2[irec]);
    acq->x3min = MIN(acq->x3min, acq->rec_x3[irec]);
    acq->x3max = MAX(acq->x3max, acq->rec_x3[irec]);
  }/* end for irec */
  if(emf->x1node[0]>acq->x1min) err("x1 origin > acq->x1min");
  if(emf->x1node[emf->nx]<acq->x1max) err("x1 end < acq->x1min");
  if(emf->x2node[0]>acq->x2min) err("x2 origin > acq->x2min");
  if(emf->x2node[emf->ny]<acq->x2max) err("x2 end < acq->x2min");
  if(emf->x3node[0]>acq->x3min) err("x3 origin > acq->x3min");
  if(emf->x3node[emf->nz]<acq->x3max) err("x3 end < acq->x3min");
  
  if(emf->verb){
    printf("lensrc=%g, nsubsrc=%d\n", acq->lensrc, acq->nsubsrc);
    printf("lenrec=%g, nsubrec=%d\n", acq->lenrec, acq->nsubrec);
    printf("src-rec domain [x1min,x1max]=[%g, %g]\n", acq->x1min, acq->x1max);
    printf("src-rec domain [x2min,x2max]=[%g, %g]\n", acq->x2min, acq->x2max);
    printf("src-rec domain [x3min,x3max]=[%g, %g]\n", acq->x3min, acq->x3max);
    printf("isrc=%d (x,y,z)=(%.2f, %.2f, %.2f)\n", 1, acq->src_x1[0], acq->src_x2[0], acq->src_x3[0]);
  }
}

/*< free the allocated variables for acquisition >*/
void acq_free(acq_t *acq)
{
  free(acq->src_x1);
  free(acq->src_x2);
  free(acq->src_x3);
  free(acq->src_azimuth);
  free(acq->src_dip);

  free(acq->rec_x1);
  free(acq->rec_x2);
  free(acq->rec_x3);
  free(acq->rec_azimuth);
  free(acq->rec_dip);

}

