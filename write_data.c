/* read and write CSEM data in binary
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 */
#include <unistd.h>
#include "cstd.h"
 

#include "acq.h"
#include "emf.h"

void write_data(acq_t *acq, emf_t *emf, char *fname, float _Complex ***dcal_fd)
/*< write synthetic data according to shot/process index >*/
{
  FILE *fp;
  int isrc, irec, ichrec, ifreq;
  float dp_re, dp_im;

  fp=fopen(fname,"w");
  if(fp==NULL) err("error opening file for writing");
  fprintf(fp, "iTx 	 iRx    ichrec  ifreq 	 emf_real 	 emf_imag\n");
  isrc = 1;//index starts from 1
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acq->nrec; irec++){
	dp_re = creal(dcal_fd[ichrec][ifreq][irec]);
	dp_im = cimag(dcal_fd[ichrec][ifreq][irec]);
	fprintf(fp, "%d \t %d \t %d \t %d \t %e \t %e\n",
		isrc, irec+1, ichrec+1, ifreq+1, dp_re, dp_im);
      }
    }
  }
  fclose(fp);
}


