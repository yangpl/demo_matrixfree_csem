/* build matrix in COO format for discretized Maxwell eqn in 2nd order
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2025 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

typedef struct{
  int nrow;
  int ncol;
  int nnz;//number of non-zeros
  int *row_ptr;//row pointers
  int *col_ind;//column indices
  complex *val;//values of the matrix A
} ccsr_t;//compressed sparse row format (compressed row storage, CRS)
ccsr_t *Acsr;

void Acsr_init(emf_t *emf, int ifreq)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int nnz, row_ind, nnzold;
  double tmp;
  complex aii;
  
  int n1 = emf->n1;
  int n2 = emf->n2;
  int n3 = emf->n3;
  emf->I_omega_mu0 = I*2.*PI*emf->freqs[ifreq]*mu0;
  emf->lEx = emf->n1*(emf->n2+1)*(emf->n3+1);//length of Ex
  emf->lEy = (emf->n1+1)*emf->n2*(emf->n3+1);//length of Ey
  emf->lEz = (emf->n1+1)*(emf->n2+1)*emf->n3;//length of Ez
  emf->lE = emf->lEx + emf->lEy + emf->lEz;//length of E=(Ex,Ey,Ez)^T

#define id1(k,j,i) (i + n1*(j + (n2+1)*(k))) 
#define id2(k,j,i) (i + (n1+1)*(j + n2*(k)))
#define id3(k,j,i) (i + (n1+1)*(j + (n2+1)*(k)))

  nnz = 0;
  //========= Eqn for Ex component ===========
  for(k=0; k<=n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<=n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<n1; i++){
	ip1 = i+1;
	im1 = i-1;

	//--------------------------------1
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1 && j<n2 && k<=n3){
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  nnz++;
	}

	//---------------------------------2
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  nnz++;
	}
	//Ex(i+0.5,j,k)
	
	//----------------------------------3
	//Ey(i+1, j-0.5, k)
	if(ip1<=n1 && jm1>=0 && k<=n3){
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  nnz++;
	}
	
	//----------------------------------4
	//Ex(i+0.5, j-1, k)
	if(jm1>=0){
	  nnz++;
	}
	//Ex(i+0.5,j,k)

	//----------------------------------5
	//Ex(i+0.5, j, k+1)
	if(kp1<=n3){
	  nnz++;
	}	
	//Ex(i+0.5,j,k)
	
	//-----------------------------------6
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1 && j<=n2 && k<n3){
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  nnz++;
	}

	//-----------------------------------7
	//Ex(i+0.5, j, k)
	//Ex(i+0.5, j, k-1)
	if(km1>=0){
	  nnz++;
	}
	//---------------------------------8
	//Ez(i+1, j, k-0.5)
	if(ip1<=n1 && j<=n2 && km1>=0){
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  nnz++;
	}
	//---------------------------------9
	//Ex(i+0.5, j, k)
	nnz++;
      }
    }
  }
  //=========== Eqn for Ey =================
  for(k=0; k<=n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<=n1; i++){
	ip1 = i+1;
	im1 = i-1;

	//-------------------------------1
	//Ez(i, j+1, k+0.5)
	if(i<=n1 && jp1<=n2 && k<n3){
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  nnz++;
	}
	
	//-------------------------------2
	//Ey(i, j+0.5, k+1)
	if(kp1<=n3){
	  nnz++;
	}	  
	//Ey(i, j+0.5, k)

	//-------------------------------3
	//Ez(i, j+1, k-0.5)
	if(i<=n1 && jp1<=n2 && km1>=0){
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  nnz++;
	}	

	//--------------------------------4
	//Ey(i, j+0.5, k)
	//Ey(i, j+0.5, k-1)
	if(km1>=0){
	  nnz++;
	}

	//---------------------------------5
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1){
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	
	//----------------------------------6
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  nnz++;
	}

	//--------------------------------7
	//Ey(i, j+0.5, k)
	//Ey(i-1, j+0.5, k)
	if(im1>=0){
	  nnz++;
	}

	//-------------------------------8
	//Ex(i-0.5, j+1, k)
	if(im1>=0 && jp1<=n2 && k<=n3){
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  nnz++;
	}

	//--------------------------------9
	//Ey(i, j+0.5, k)
	nnz++;
      }
    }
  }
  //========= Eqn for Ez =====================
  for(k=0; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<=n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<=n1; i++){
	ip1 = i+1;
	im1 = i-1;

	//----------------------------1
	//Ex(i+0.5, j, k+1)
	if(i<n1 && j<=n2 && kp1<=n3){
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  nnz++;
	}
	//-----------------------------2
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1){
	  nnz++;
	}
	//Ez(i, j, k+0.5)

	//----------------------------3
	//Ex(i-0.5, j, k+1)
	if(im1>=0 && j<=n2 && kp1<=n3){
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  nnz++;
	}

	//----------------------------4
	//Ez(i, j, k+0.5)
	//Ez(i-1, j, k+0.5)
	if(im1>=0){
	  nnz++;
	}

	//----------------------------5
	//Ez(i, j+1, k+0.5)
	if(jp1<=n2){
	  nnz++;
	}
	//Ez(i, j, k+0.5)

	//-----------------------------6
	//Ey(i, j+0.5, k+1)
	if(i<=n1 && j<n2 && kp1<=n3){
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  nnz++;
	}

	//----------------------------7	
	//Ez(i, j, k+0.5)
	//Ez(i, j-1, k+0.5)
	if(jm1>=0){
	  nnz++;
	}

	//----------------------------8
	//Ey(i, j-0.5, k+1)
	if(i<=n1 && jm1>=0 && kp1<=n3){
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  nnz++;
	}

	//----------------------------9
	//Ez(i, j, k+0.5)
	nnz++;
      }//end for i
    }//end for j
  }//end for k

  Acsr = malloc(sizeof(ccsr_t));
  Acsr->nrow = emf->lEx + emf->lEy + emf->lEz;
  Acsr->ncol = emf->lEx + emf->lEy + emf->lEz;
  Acsr->nnz = nnz;
  Acsr->row_ptr = malloc((Acsr->nrow+1)*sizeof(int));
  Acsr->col_ind = malloc(Acsr->nnz*sizeof(int));
  Acsr->val = malloc(Acsr->nnz*sizeof(complex));

  nnz = 0;
  Acsr->row_ptr[0] = 0;
  //========= Eqn for Ex component ===========
  for(k=0; k<=n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<=n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<n1; i++){
	ip1 = i+1;
	im1 = i-1;

	aii = 0;//diagonal coefficient
	//--------------------------------1
	tmp = 1./(emf->d2[j]*emf->d1s[i]);
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1 && j<n2 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,j,ip1) + emf->lEx;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,j,i) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//---------------------------------2
	tmp = 1./(emf->d2[j]*emf->d2s[j]);
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,jp1,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i+0.5,j,k)
	aii += tmp;
	
	//----------------------------------3
	tmp = 1./(emf->d2[j]*emf->d1s[i]);
	//Ey(i+1, j-0.5, k)
	if(ip1<=n1 && jm1>=0 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,jm1,ip1) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,jm1,i) + emf->lEx;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	
	//----------------------------------4
	tmp = 1./(emf->d2[j]*emf->d2s[MAX(jm1,0)]);
	//Ex(i+0.5, j-1, k)
	if(jm1>=0){
	  Acsr->col_ind[nnz] = id1(k,jm1,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i+0.5,j,k)
	aii += tmp;

	//----------------------------------5
	tmp = 1./(emf->d3[k]*emf->d3s[k]);
	//Ex(i+0.5, j, k+1)
	if(kp1<=n3){
	  Acsr->col_ind[nnz] = id1(kp1,j,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}	
	//Ex(i+0.5,j,k)
	aii += tmp;
	
	//-----------------------------------6
	tmp = 1./(emf->d3[k]*emf->d1s[i]);
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1 && j<=n2 && k<n3){
	  Acsr->col_ind[nnz] = id3(k,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  Acsr->col_ind[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//-----------------------------------7
	tmp = 1./(emf->d3[k]*emf->d3s[MAX(km1,0)]);
	//Ex(i+0.5, j, k)
	aii += tmp;
	//Ex(i+0.5, j, k-1)
	if(km1>=0){
	  Acsr->col_ind[nnz] = id1(km1,j,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//---------------------------------8
	tmp = 1./(emf->d3[k]*emf->d1s[i]);
	//Ez(i+1, j, k-0.5)
	if(ip1<=n1 && j<=n2 && km1>=0){
	  Acsr->col_ind[nnz] = id3(km1,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  Acsr->col_ind[nnz] = id3(km1,j,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//---------------------------------9
	//Ex(i+0.5, j, k)
	row_ind = id1(k,j,i);//row index
	Acsr->row_ptr[row_ind+1] = nnz;//row pointer
	Acsr->col_ind[nnz] = id1(k,j,i);//column index
	Acsr->val[nnz] = aii - emf->I_omega_mu0*emf->sigma11[k][j][i];//matrix element at ii-row, Acsr->col[nnz]-column
	nnz++;
      }
    }
  }
  printf("#nnz for Ex=%d\n", nnz);
  nnzold = nnz;
  //=========== Eqn for Ey =================
  for(k=0; k<=n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<=n1; i++){
	ip1 = i+1;
	im1 = i-1;

	aii = 0;
	//-------------------------------1
	tmp = 1./(emf->d3[k]*emf->d2s[j]);
	//Ez(i, j+1, k+0.5)
	if(i<=n1 && jp1<=n2 && k<n3){
	  Acsr->col_ind[nnz] = id3(k,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  Acsr->col_ind[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	
	//-------------------------------2
	tmp = 1./(emf->d3[k]*emf->d3s[k]);	
	//Ey(i, j+0.5, k+1)
	if(kp1<=n3){
	  Acsr->col_ind[nnz] = id2(kp1,j,i) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}	  
	//Ey(i, j+0.5, k)
	aii += tmp;

	//-------------------------------3
	tmp = 1./(emf->d3[k]*emf->d2s[j]);
	//Ez(i, j+1, k-0.5)
	if(i<=n1 && jp1<=n2 && km1>=0){
	  Acsr->col_ind[nnz] = id3(km1,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  Acsr->col_ind[nnz] = id3(km1,j,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}	

	//--------------------------------4
	tmp = 1./(emf->d3[k]*emf->d3s[MAX(km1,0)]);
	//Ey(i, j+0.5, k)
	aii += tmp;
	//Ey(i, j+0.5, k-1)
	if(km1>=0){
	  Acsr->col_ind[nnz] = id2(km1,j,i) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//---------------------------------5
	tmp = 1./(emf->d1[i]*emf->d1s[i]);
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1){
	  Acsr->col_ind[nnz] = id2(k,j,ip1) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	aii += tmp;
	
	//----------------------------------6
	tmp = 1./(emf->d1[i]*emf->d2s[j]);
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,jp1,i);//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,j,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//--------------------------------7
	tmp = 1./(emf->d1[i]*emf->d1s[MAX(im1,0)]);
	//Ey(i, j+0.5, k)
	aii += tmp;
	//Ey(i-1, j+0.5, k)
	if(im1>=0){
	  Acsr->col_ind[nnz] = id2(k,j,im1) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//-------------------------------8
	tmp = 1./(emf->d1[i]*emf->d2s[j]);
	//Ex(i-0.5, j+1, k)
	if(im1>=0 && jp1<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,jp1,im1);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,j,im1);//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}

	//--------------------------------9
	//Ey(i, j+0.5, k)
	row_ind = id2(k,j,i) + emf->lEx;//row index
	Acsr->row_ptr[row_ind + 1] = nnz;//row pointer
	Acsr->col_ind[nnz] = id2(k,j,i) + emf->lEx;//column index
	Acsr->val[nnz] = aii - emf->I_omega_mu0*emf->sigma22[k][j][i];
	nnz++;
      }
    }
  }
  printf("#nnz for Ey=%d\n", nnz-nnzold);
  nnzold = nnz;
  //========= Eqn for Ez =====================
  for(k=0; k<n3; k++){
    kp1 = k+1;
    km1 = k-1;
    for(j=0; j<=n2; j++){
      jp1 = j+1;
      jm1 = j-1;
      for(i=0; i<=n1; i++){
	ip1 = i+1;
	im1 = i-1;

	aii = 0;
	//----------------------------1
	tmp = 1./(emf->d1[i]*emf->d3s[k]);
	//Ex(i+0.5, j, k+1)
	if(i<n1 && j<=n2 && kp1<=n3){
	  Acsr->col_ind[nnz] = id1(kp1,j,i);//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,j,i);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//-----------------------------2
	tmp = 1./(emf->d1[i]*emf->d1s[i]);
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1){
	  Acsr->col_ind[nnz] = id3(k,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	aii += tmp;

	//----------------------------3
	tmp = 1./(emf->d1[i]*emf->d3s[k]);
	//Ex(i-0.5, j, k+1)
	if(im1>=0 && j<=n2 && kp1<=n3){
	  Acsr->col_ind[nnz] = id1(kp1,j,im1);//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  Acsr->col_ind[nnz] = id1(k,j,im1);//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}

	//----------------------------4
	tmp = 1./(emf->d1[i]*emf->d1s[MAX(im1,0)]);
	//Ez(i, j, k+0.5)
	aii += tmp;
	//Ez(i-1, j, k+0.5)
	if(im1>=0){
	  Acsr->col_ind[nnz] = id3(k,j,im1) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//----------------------------5
	tmp = 1./(emf->d2[j]*emf->d2s[j]);
	//Ez(i, j+1, k+0.5)
	if(jp1<=n2){
	  Acsr->col_ind[nnz] = id3(k,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	aii += tmp;

	//-----------------------------6
	tmp = 1./(emf->d2[j]*emf->d3s[k]);
	//Ey(i, j+0.5, k+1)
	if(i<=n1 && j<n2 && kp1<=n3){
	  Acsr->col_ind[nnz] = id2(kp1,j,i) + emf->lEx;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,j,i) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//----------------------------7	
	tmp = 1./(emf->d2[j]*emf->d2s[MAX(jm1,0)]);
	//Ez(i, j, k+0.5)
	aii += tmp;
	//Ez(i, j-1, k+0.5)
	if(jm1>=0){
	  Acsr->col_ind[nnz] = id3(k,jm1,i) + emf->lEx + emf->lEy;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}

	//----------------------------8
	tmp = 1./(emf->d2[j]*emf->d3s[k]);
	//Ey(i, j-0.5, k+1)
	if(i<=n1 && jm1>=0 && kp1<=n3){
	  Acsr->col_ind[nnz] = id2(kp1,jm1,i) + emf->lEx;//column index
	  Acsr->val[nnz] = -tmp;
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  Acsr->col_ind[nnz] = id2(k,jm1,i) + emf->lEx;//column index
	  Acsr->val[nnz] = tmp;
	  nnz++;
	}

	//----------------------------9
	//Ez(i, j, k+0.5)
	row_ind = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	Acsr->row_ptr[row_ind + 1] = nnz;//row pointer
	Acsr->col_ind[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	Acsr->val[nnz] = aii - emf->I_omega_mu0*emf->sigma33[k][j][i];
	nnz++;
      }//end for i
    }//end for j
  }//end for k
  printf("#nnz for Ez=%d\n", nnz-nnzold);

#undef id1
#undef id2
#undef id3

}

void Acsr_free()
{
  free1int(Acsr->row_ptr);
  free1int(Acsr->col_ind);
  free1complex(Acsr->val);
  free(Acsr);
}

//y=Ax
void Acsr_apply(int n, complex *x, complex *y)
{
  int i, j, k;

  for(i=0; i<Acsr->nrow; i++){
    y[i] = 0;
    for(k=Acsr->row_ptr[i]; k<Acsr->row_ptr[i+1]; k++){
      j = Acsr->col_ind[k];//a_ij=Acsr->val[k]
      y[i] += Acsr->val[k]*x[j];//y_i += a_ij*x_j
    }
  }  
}


//LU-SGS=SGS initialized by 0, y=Minv*x (or My=x)
void Acsr_lusgs_apply(int n, complex *x, complex *y)
{
  int i, j, k;
  complex s, aii;
  
  //-------------------------------------------------------
  memset(y, 0, Acsr->nrow*sizeof(complex));
  //forward Gauss-Seidel
  for(i=0; i<Acsr->nrow; i++){
    s = 0;
    aii = 1;
    for(k=Acsr->row_ptr[i]; k<Acsr->row_ptr[i+1]; k++){
      j = Acsr->col_ind[k];
      if(j==i) aii = Acsr->val[k];
      else s += Acsr->val[k]*y[j];
    }
    y[i] = (x[i] - s)/aii;
  }
  //backward Gauss-Seidel
  for(i=Acsr->nrow-1; i>=0; i--){
    s = 0;
    aii = 1;
    for(k=Acsr->row_ptr[i]; k<Acsr->row_ptr[i+1]; k++){
      j = Acsr->col_ind[k];
      if(j==i) aii = Acsr->val[k];
      else s += Acsr->val[k]*y[j];
    }
    y[i] = (x[i] - s)/aii;
  }
}
