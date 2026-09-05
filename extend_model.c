#include "cstd.h"
#include "acq.h"
#include "emf.h"

double create_nugrid(int n, double len, double dx, double *x, int istretch);

/*< generate staggered grid coordinate >*/
void generate_staggered_xs_dx(int n1, double *x1, double *x1s, double *d1, double *d1s)
{
  int i;
  
  for(i=0; i<n1; i++) {
    d1s[i] = x1[i+1] - x1[i];
    x1s[i] = 0.5*(x1[i] + x1[i+1]);
  }
  d1s[n1] = d1s[n1-1];
  x1s[n1] = x1[n1] + 0.5*d1s[n1-1];

  for(i=1; i<=n1; i++) d1[i] = x1s[i] - x1s[i-1];
  d1[0] = d1[1];
}


void extend_model_init(acq_t *acq, emf_t *emf)
{
  int i, j, k;
  int im1, jm1, km1;
  double dx, *x;
  double ***sig11, ***sig22, ***rho33;
  FILE *fp;

  if(emf->verb)  printf("------------ extend model -------------\n");
  emf->n1 = emf->nx + 2*emf->nb;//total number of FD cells in x direction
  emf->n2 = emf->ny + 2*emf->nb;//total number of FD cells in y direction 
  emf->n3 = emf->nz + 2*emf->nb;//total number of FD cells in z direction 
  emf->d1min = emf->x1node[1]-emf->x1node[0];
  for(i=0; i<emf->nx; i++) emf->d1min = MIN(emf->d1min, emf->x1node[i+1]-emf->x1node[i]);
  emf->d2min = emf->x2node[1]-emf->x2node[0];
  for(j=0; j<emf->ny; j++) emf->d2min = MIN(emf->d2min, emf->x2node[j+1]-emf->x2node[j]);
  emf->d3min = emf->x3node[1]-emf->x3node[0];
  for(k=0; k<emf->nz; k++) emf->d3min = MIN(emf->d3min, emf->x3node[k+1]-emf->x3node[k]);

  if(!getparfloat("lextend", &emf->lextend)) emf->lextend = (emf->nb>0)?40e3:0;//domain extension
  if(emf->nb>0){
    dx = emf->x1node[emf->nx]-emf->x1node[0];
    emf->lextend = MAX(emf->lextend, 2.2*dx);
    dx = emf->x2node[emf->ny]-emf->x2node[0];
    emf->lextend = MAX(emf->lextend, 2.2*dx);
    dx = emf->x3node[emf->nz]-emf->x3node[0];
    emf->lextend = MAX(emf->lextend, 2.2*dx);
  }  
  if(emf->verb){
    printf("[d1min, d2min, d3min]=[%g, %g, %g]\n", emf->d1min, emf->d2min, emf->d3min);
    printf("lextend=%g \n", emf->lextend);
  }


  emf->x1 = alloc1double(emf->n1+1);
  emf->x2 = alloc1double(emf->n2+1);
  emf->x3 = alloc1double(emf->n3+1);
  emf->x1s = alloc1double(emf->n1+1);
  emf->x2s = alloc1double(emf->n2+1);
  emf->x3s = alloc1double(emf->n3+1);

  emf->d1 = alloc1double(emf->n1+1);
  emf->d2 = alloc1double(emf->n2+1);
  emf->d3 = alloc1double(emf->n3+1);
  emf->d1s = alloc1double(emf->n1+1);
  emf->d2s = alloc1double(emf->n2+1);
  emf->d3s = alloc1double(emf->n3+1);

  for(i=0; i<=emf->nx; i++) emf->x1[i+emf->nb] = emf->x1node[i];
  for(j=0; j<=emf->ny; j++) emf->x2[j+emf->nb] = emf->x2node[j];
  for(k=0; k<=emf->nz; k++) emf->x3[k+emf->nb] = emf->x3node[k];
  if(emf->nb>0){//if we extend domain nb cells on each side
    x = alloc1double(emf->nb+1);
    //left of x
    dx = emf->x1node[1] - emf->x1node[0];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);
    for(i=0; i<emf->nb; i++) emf->x1[emf->nb-1-i] = emf->x1[emf->nb-i] - dx*pow(emf->rstretch,i);
    if(emf->verb) printf("left stretching factor: r=%g\n", emf->rstretch);
    //right of x
    dx = emf->x1node[emf->nx] - emf->x1node[emf->nx-1];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);
    for(i=0; i<emf->nb; i++) emf->x1[emf->nb+emf->nx+1+i] = emf->x1[emf->nb+emf->nx+i] + dx*pow(emf->rstretch,i);
    if(emf->verb) printf("right stretching factor: r=%g\n", emf->rstretch);
    //front of y
    dx = emf->x2node[1] - emf->x2node[0];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);
    for(j=0; j<emf->nb; j++) emf->x2[emf->nb-1-j] = emf->x2[emf->nb-j] - dx*pow(emf->rstretch,j);
    if(emf->verb) printf("front stretching factor: r=%g\n", emf->rstretch);
    //rear of y
    dx = emf->x2node[emf->ny] - emf->x2node[emf->ny-1];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);
    for(j=0; j<emf->nb; j++) emf->x2[emf->nb+emf->ny+1+j] = emf->x2[emf->nb+emf->ny+j] + dx*pow(emf->rstretch,j);
    if(emf->verb) printf("rear stretching factor: r=%g\n", emf->rstretch);
    //top of z
    dx = emf->x3node[1] - emf->x3node[0];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);//determine stretching factor using thickness of air layer 
    for(k=0; k<emf->nb; k++) emf->x3[emf->nb-1-k] = emf->x3[emf->nb-k] - dx*pow(emf->rstretch,k);
    if(emf->verb) printf("top stretching factor: r=%g\n", emf->rstretch);
    //bottom of z
    dx = emf->x3node[emf->nz] - emf->x3node[emf->nz-1];
    emf->rstretch = create_nugrid(emf->nb, emf->lextend, dx, x, emf->istretch);
    for(k=0; k<emf->nb; k++) emf->x3[emf->nb+emf->nz+1+k] = emf->x3[emf->nb+emf->nz+k] + dx*pow(emf->rstretch,k);
    if(emf->verb) printf("bottom stretching factor: r=%g\n", emf->rstretch);
    free1double(x);
  }
  if(emf->verb){
    printf("n1=nx+2*nb=%d\n", emf->n1);
    printf("n2=ny+2*nb=%d\n", emf->n2);
    printf("n3=nz+2*nb=%d\n", emf->n3);
    printf("extended domain [x1min,x1max]=[%g, %g]\n", emf->x1[0], emf->x1[emf->n1]);
    printf("extended domain [x2min,x2max]=[%g, %g]\n", emf->x2[0], emf->x2[emf->n2]);
    printf("extended domain [x3min,x3max]=[%g, %g]\n", emf->x3[0], emf->x3[emf->n3]);
    printf("extended Lx=%g\n", emf->x1[emf->n1]-emf->x1[0]);
    printf("extended Ly=%g\n", emf->x2[emf->n2]-emf->x2[0]);
    printf("extended Lz=%g\n", emf->x3[emf->n3]-emf->x3[0]);
  }

  //this will be used for extraction of EM data at receiver locations
  generate_staggered_xs_dx(emf->n1, emf->x1, emf->x1s, emf->d1, emf->d1s);
  generate_staggered_xs_dx(emf->n2, emf->x2, emf->x2s, emf->d2, emf->d2s);
  generate_staggered_xs_dx(emf->n3, emf->x3, emf->x3s, emf->d3, emf->d3s);

  //extend input model parameters, still sitting in the center of FD volume
  sig11 = alloc3double(emf->n1, emf->n2, emf->n3);
  sig22 = alloc3double(emf->n1, emf->n2, emf->n3);
  rho33 = alloc3double(emf->n1, emf->n2, emf->n3);
  //fill the interior part
  for(k=0; k<emf->nz; k++){
    for(j=0; j<emf->ny; j++){
      for(i=0; i<emf->nx; i++){
	sig11[k+emf->nb][j+emf->nb][i+emf->nb] = 1./emf->rho11[k][j][i];
	sig22[k+emf->nb][j+emf->nb][i+emf->nb] = 1./emf->rho22[k][j][i];
	rho33[k+emf->nb][j+emf->nb][i+emf->nb] = emf->rho33[k][j][i];
      }
    }
  }
  if(emf->nb>0){//if we extend the domain nb cells on each side
    for(k=0; k<emf->n3; k++){
      for(j=0; j<emf->n2; j++){
	for(i=0; i<emf->nb; i++){
	  sig11[k][j][i] = sig11[k][j][emf->nb];
	  sig11[k][j][emf->n1-1-i] = sig11[k][j][emf->n1-1-emf->nb];
	  sig22[k][j][i] = sig22[k][j][emf->nb];
	  sig22[k][j][emf->n1-1-i] = sig22[k][j][emf->n1-1-emf->nb];
	  rho33[k][j][i] = rho33[k][j][emf->nb];
	  rho33[k][j][emf->n1-1-i] = rho33[k][j][emf->n1-1-emf->nb];
	}
      }
    }
    for(k=0; k<emf->n3; k++){
      for(j=0; j<emf->nb; j++){
	for(i=0; i<emf->n1; i++){
	  sig11[k][j][i] = sig11[k][emf->nb][i];
	  sig11[k][emf->n2-1-j][i] = sig11[k][emf->n2-1-emf->nb][i];
	  sig22[k][j][i] = sig22[k][emf->nb][i];
	  sig22[k][emf->n2-1-j][i] = sig22[k][emf->n2-1-emf->nb][i];
	  rho33[k][j][i] = rho33[k][emf->nb][i];
	  rho33[k][emf->n2-1-j][i] = rho33[k][emf->n2-1-emf->nb][i];
	}
      }
    }
    for(k=0; k<emf->nb; k++){
      for(j=0; j<emf->n2; j++){
	for(i=0; i<emf->n1; i++){
	  sig11[k][j][i] = 1./emf->rho_air;//fill in with air
	  sig11[emf->n3-1-k][j][i] = sig11[emf->n3-1-emf->nb][k][i];
	  sig22[k][j][i] = 1./emf->rho_air;//fill in with air
	  sig22[emf->n3-1-k][j][i] = sig22[emf->n3-1-emf->nb][k][i];
	  rho33[k][j][i] = emf->rho_air;//fill in with air
	  rho33[emf->n3-1-k][j][i] = rho33[emf->n3-1-emf->nb][k][i];
	}
      }
    }  
  }//end if
  
  /*------------ read resistivity models ------------------*/
  emf->sigma11 = alloc3double(emf->n1, emf->n2+1, emf->n3+1);
  emf->sigma22 = alloc3double(emf->n1+1, emf->n2, emf->n3+1);
  emf->sigma33 = alloc3double(emf->n1+1, emf->n2+1, emf->n3);
  //average sigma11
  for(k=1; k<emf->n3; k++){
    km1 = k-1;
    for(j=1; j<emf->n2; j++){
      jm1 = j-1;
      for(i=0; i<emf->n1; i++){
	emf->sigma11[k][j][i] = 0.25*((emf->d2s[j]*sig11[k][j][i] + emf->d2s[jm1]*sig11[k][jm1][i])*emf->d3s[k] +
				      (emf->d2s[j]*sig11[km1][j][i] + emf->d2s[jm1]*sig11[km1][jm1][i])*emf->d3s[km1])/(emf->d2[j]*emf->d3[k]);
      }
    }
  }
  for(j=0; j<=emf->n2; j++){
    for(i=0; i<emf->n1; i++){
      emf->sigma11[0][j][i] = emf->sigma11[1][j][i];
      emf->sigma11[emf->n3][j][i] = emf->sigma11[emf->n3-1][j][i];
    }
  }
  for(k=0; k<=emf->n3; k++){
    for(i=0; i<emf->n1; i++){
      emf->sigma11[k][0][i] = emf->sigma11[k][1][i];
      emf->sigma11[k][emf->n2][i] = emf->sigma11[k][emf->n2-1][i];
    }
  }
  //average sigma22
  for(k=1; k<emf->n3; k++){
    km1 = k-1;
    for(j=0; j<emf->n2; j++){
      for(i=1; i<emf->n1; i++){
	im1 = i-1;
	emf->sigma22[k][j][i] = 0.25*((emf->d1s[i]*sig22[k][j][i] + emf->d1s[im1]*sig22[k][j][im1])*emf->d3s[k] +
				      (emf->d1s[i]*sig22[km1][j][i] + emf->d1s[im1]*sig22[km1][j][im1])*emf->d3s[km1])/(emf->d1[i]*emf->d3[k]);
      }
    }
  }
  for(j=0; j<emf->n2; j++){
    for(i=0; i<=emf->n1; i++){
      emf->sigma22[0][j][i] = emf->sigma22[1][j][i];
      emf->sigma22[emf->n3][j][i] = emf->sigma22[emf->n3-1][j][i];
    }
  }
  for(k=0; k<=emf->n3; k++){
    for(j=0; j<emf->n2; j++){
      emf->sigma22[k][j][0] = emf->sigma22[k][j][1];
      emf->sigma22[k][j][emf->n1] = emf->sigma22[k][j][emf->n1-1];
    }
  }
  //average rho33
  for(k=0; k<emf->n3; k++){
    for(j=1; j<emf->n2; j++){
      jm1 = j-1;
      for(i=1; i<emf->n1; i++){
	im1 = i-1;
	emf->sigma33[k][j][i] = 0.25*((emf->d1s[i]*rho33[k][j][i] + emf->d1s[im1]*rho33[k][j][im1])*emf->d2s[j] +
				      (emf->d1s[i]*rho33[k][jm1][i] + emf->d1s[im1]*rho33[k][jm1][im1])*emf->d2s[jm1])/(emf->d1[i]*emf->d2[j]);
	emf->sigma33[k][j][i] = 1./emf->sigma33[k][j][i];//convert rho33 to sigma33
      }
    }
  }
  for(k=0; k<emf->n3; k++){
    for(j=0; j<=emf->n2; j++){
      emf->sigma33[k][j][0] = emf->sigma33[k][j][1];
      emf->sigma33[k][j][emf->n1] = emf->sigma33[k][j][emf->n1-1];
    }
  }
  for(k=0; k<emf->n3; k++){
    for(i=0; i<=emf->n1; i++){
      emf->sigma33[k][0][i] = emf->sigma33[k][1][i];
      emf->sigma33[k][emf->n2][i] = emf->sigma33[k][emf->n2-1][i];
    }
  }

  if(!getparint("param_extended", &emf->param_extended)) emf->param_extended = 0;
  if(emf->param_extended && emf->verb){
    for(k=0; k<emf->n3; k++){
      for(j=0; j<emf->n2; j++){
	for(i=0; i<emf->n1; i++){
	  //convert sig11,sig22 to rho11,rho22
	  sig11[k][j][i] = 1./sig11[k][j][i];
	  sig22[k][j][i] = 1./sig22[k][j][i];
	}
      }
    }

    fp = fopen("frho11_extended", "wb");
    fwrite(&sig11[0][0][0], emf->n1*emf->n2*emf->n3*sizeof(double), 1, fp);
    fclose(fp);

    fp = fopen("frho22_extended", "wb");
    fwrite(&sig22[0][0][0], emf->n1*emf->n2*emf->n3*sizeof(double), 1, fp);
    fclose(fp);
    
    fp = fopen("frho33_extended", "wb");
    fwrite(&rho33[0][0][0], emf->n1*emf->n2*emf->n3*sizeof(double), 1, fp);
    fclose(fp);
    
    fp = fopen("fx1_extended", "wb");
    fwrite(emf->x1, (emf->n1+1)*sizeof(double), 1, fp);
    fclose(fp);

    fp = fopen("fx2_extended", "wb");
    fwrite(emf->x2, (emf->n2+1)*sizeof(double), 1, fp);
    fclose(fp);

    fp = fopen("fx3_extended", "wb");
    fwrite(emf->x3, (emf->n3+1)*sizeof(double), 1, fp);
    fclose(fp);
  }
  
  free3double(sig11);
  free3double(sig22);
  free3double(rho33);
}

void extend_model_free(emf_t *emf)
{
  free1double(emf->x1);
  free1double(emf->x2);
  free1double(emf->x3);

  free1double(emf->x1s);
  free1double(emf->x2s);
  free1double(emf->x3s);

  free1double(emf->d1);
  free1double(emf->d2);
  free1double(emf->d3);

  free1double(emf->d1s);
  free1double(emf->d2s);
  free1double(emf->d3s);

  free3double(emf->sigma11);
  free3double(emf->sigma22);
  free3double(emf->sigma33);
  
}

