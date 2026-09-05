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
#include "mpi.h"
#include "zmumps_c.h"

typedef struct{
  int nrow;
  int ncol;
  int nnz;//number of non-zeros
  int *row;//row indices
  int *col;//column indices
  complex *val;//values of the matrix A
} ccoo_t;//coordinate format (COO)
ccoo_t *Acoo;

void Acoo_init(emf_t *emf, int ifreq)
{
  int i, j, k;
  int nnz, nnzold;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
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
	//Ez(i, j-1, k+0.5)
	if(jm1>=0){
	  nnz++;
	}
	//Ez(i, j, k+0.5)

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

  Acoo = malloc(sizeof(ccoo_t));
  Acoo->nrow = emf->lEx + emf->lEy + emf->lEz;
  Acoo->ncol = emf->lEx + emf->lEy + emf->lEz;
  Acoo->nnz = nnz;
  Acoo->row = malloc(nnz*sizeof(int));
  Acoo->col = malloc(nnz*sizeof(int));
  Acoo->val = malloc(nnz*sizeof(complex));
  
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

	aii = 0;//diagonal coefficient
	//--------------------------------1
	tmp = 1./(emf->d2[j]*emf->d1s[i]);
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1 && j<n2 && k<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id2(k,j,ip1) + emf->lEx;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id2(k,j,i) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//---------------------------------2
	tmp = 1./(emf->d2[j]*emf->d2s[j]);
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id1(k,jp1,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ex(i+0.5,j,k)
	  aii += tmp;
	}
	
	//----------------------------------3
	tmp = 1./(emf->d2[j]*emf->d1s[i]);
	//Ey(i+1, j-0.5, k)
	if(ip1<=n1 && jm1>=0 && k<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id2(k,jm1,ip1) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id2(k,jm1,i) + emf->lEx;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	
	//----------------------------------4
	tmp = 1./(emf->d2[j]*emf->d2s[MAX(jm1,0)]);
	//Ex(i+0.5, j-1, k)
	if(jm1>=0){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id1(k,jm1,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ex(i+0.5,j,k)
	  aii += tmp;
	}

	//----------------------------------5
	tmp = 1./(emf->d3[k]*emf->d3s[k]);
	//Ex(i+0.5, j, k+1)
	if(kp1<=n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id1(kp1,j,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ex(i+0.5,j,k)
	  aii += tmp;
	}	
	
	//-----------------------------------6
	tmp = 1./(emf->d3[k]*emf->d1s[i]);
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1 && j<=n2 && k<n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id3(k,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//-----------------------------------7
	tmp = 1./(emf->d3[k]*emf->d3s[MAX(km1,0)]);
	//Ex(i+0.5, j, k-1)
	if(km1>=0){
	  //Ex(i+0.5, j, k)
	  aii += tmp;
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id1(km1,j,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//---------------------------------8
	tmp = 1./(emf->d3[k]*emf->d1s[i]);
	//Ez(i+1, j, k-0.5)
	if(ip1<=n1 && j<=n2 && km1>=0){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id3(km1,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  Acoo->row[nnz] = id1(k,j,i);//row index
	  Acoo->col[nnz] = id3(km1,j,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//---------------------------------9
	//Ex(i+0.5, j, k)
	Acoo->row[nnz] = id1(k,j,i);//row index
	Acoo->col[nnz] = id1(k,j,i);//column index
	Acoo->val[nnz] = aii - emf->I_omega_mu0*emf->sigma11[k][j][i];//matrix element at ii-row, Acoo->col[nnz]-column
	nnz++;
      }
    }
  }
  if(emf->verb) printf("#nnz-Ex=%d\n", nnz);
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

	aii = 0;//diagonal coefficient
	//-------------------------------1
	tmp = 1./(emf->d3[k]*emf->d2s[j]);
	//Ez(i, j+1, k+0.5)
	if(i<=n1 && jp1<=n2 && k<n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id3(k,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ez(i, j, k+0.5)
	if(i<=n1 && j<=n2 && k<n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	
	//-------------------------------2
	tmp = 1./(emf->d3[k]*emf->d3s[k]);	
	//Ey(i, j+0.5, k+1)
	if(kp1<=n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id2(kp1,j,i) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ey(i, j+0.5, k)
	  aii += tmp;
	}	  

	//-------------------------------3
	tmp = 1./(emf->d3[k]*emf->d2s[j]);
	//Ez(i, j+1, k-0.5)
	if(i<=n1 && jp1<=n2 && km1>=0){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id3(km1,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ez(i, j, k-0.5)
	if(i<=n1 && j<=n2 && km1>=0){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id3(km1,j,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}	

	//--------------------------------4
	tmp = 1./(emf->d3[k]*emf->d3s[MAX(km1,0)]);
	//Ey(i, j+0.5, k-1)
	if(km1>=0){
	  //Ey(i, j+0.5, k)
	  aii += tmp;
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id2(km1,j,i) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//---------------------------------5
	tmp = 1./(emf->d1[i]*emf->d1s[i]);
	//Ey(i+1, j+0.5, k)
	if(ip1<=n1){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id2(k,j,ip1) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ey(i, j+0.5, k)
	  aii += tmp;
	}
	
	//----------------------------------6
	tmp = 1./(emf->d1[i]*emf->d2s[j]);
	//Ex(i+0.5, j+1, k)
	if(i<n1 && jp1<=n2 && k<=n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id1(k,jp1,i);//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id1(k,j,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//--------------------------------7
	tmp = 1./(emf->d1[i]*emf->d1s[MAX(im1,0)]);
	//Ey(i-1, j+0.5, k)
	if(im1>=0){
	  //Ey(i, j+0.5, k)
	  aii += tmp;
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id2(k,j,im1) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//-------------------------------8
	tmp = 1./(emf->d1[i]*emf->d2s[j]);
	//Ex(i-0.5, j+1, k)
	if(im1>=0 && jp1<=n2 && k<=n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id1(k,jp1,im1);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	  Acoo->col[nnz] = id1(k,j,im1);//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}

	//--------------------------------9
	//Ey(i, j+0.5, k)
	Acoo->row[nnz] = id2(k,j,i) + emf->lEx;//row index
	Acoo->col[nnz] = id2(k,j,i) + emf->lEx;//column index
	Acoo->val[nnz] = aii - emf->I_omega_mu0*emf->sigma22[k][j][i];
	nnz++;
      }
    }
  }
  if(emf->verb) printf("#nnz-Ey=%d\n", nnz-nnzold);
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

	aii = 0;//diagonal coefficient
	//----------------------------1
	tmp = 1./(emf->d1[i]*emf->d3s[k]);
	//Ex(i+0.5, j, k+1)
	if(i<n1 && j<=n2 && kp1<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id1(kp1,j,i);//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ex(i+0.5, j, k)
	if(i<n1 && j<=n2 && k<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id1(k,j,i);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//-----------------------------2
	tmp = 1./(emf->d1[i]*emf->d1s[i]);
	//Ez(i+1, j, k+0.5)
	if(ip1<=n1){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id3(k,j,ip1) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ez(i, j, k+0.5)
	  aii += tmp;
	}

	//----------------------------3
	tmp = 1./(emf->d1[i]*emf->d3s[k]);
	//Ex(i-0.5, j, k+1)
	if(im1>=0 && j<=n2 && kp1<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id1(kp1,j,im1);//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ex(i-0.5, j, k)
	if(im1>=0 && j<=n2 && k<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id1(k,j,im1);//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}

	//----------------------------4
	tmp = 1./(emf->d1[i]*emf->d1s[MAX(im1,0)]);
	//Ez(i-1, j, k+0.5)
	if(im1>=0){
	  //Ez(i, j, k+0.5)
	  aii += tmp;
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id3(k,j,im1) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//----------------------------5
	tmp = 1./(emf->d2[j]*emf->d2s[j]);
	//Ez(i, j+1, k+0.5)
	if(jp1<=n2){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id3(k,jp1,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ez(i, j, k+0.5)
	  aii += tmp;
	}

	//-----------------------------6
	tmp = 1./(emf->d2[j]*emf->d3s[k]);
	//Ey(i, j+0.5, k+1)
	if(i<=n1 && j<n2 && kp1<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id2(kp1,j,i) + emf->lEx;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}
	//Ey(i, j+0.5, k)
	if(i<=n1 && j<n2 && k<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id2(k,j,i) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}

	//----------------------------7	
	tmp = 1./(emf->d2[j]*emf->d2s[MAX(jm1,0)]);
	//Ez(i, j-1, k+0.5)
	if(jm1>=0){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id3(k,jm1,i) + emf->lEx + emf->lEy;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	  //Ez(i, j, k+0.5)
	  aii += tmp;
	}

	//----------------------------8
	tmp = 1./(emf->d2[j]*emf->d3s[k]);
	//Ey(i, j-0.5, k+1)
	if(i<=n1 && jm1>=0 && kp1<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id2(kp1,jm1,i) + emf->lEx;//column index
	  Acoo->val[nnz] = -tmp;
	  nnz++;
	}
	//Ey(i, j-0.5, k)
	if(i<=n1 && jm1>=0 && k<=n3){
	  Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	  Acoo->col[nnz] = id2(k,jm1,i) + emf->lEx;//column index
	  Acoo->val[nnz] = tmp;
	  nnz++;
	}

	//----------------------------9
	//Ez(i, j, k+0.5)
	Acoo->row[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//row index
	Acoo->col[nnz] = id3(k,j,i) + emf->lEx + emf->lEy;//column index
	Acoo->val[nnz] = aii - emf->I_omega_mu0*emf->sigma33[k][j][i];
	nnz++;
      }//end for i
    }//end for j
  }//end for k
  if(emf->verb) printf("#nnz-Ez=%d\n", nnz-nnzold);
  
#undef id1
#undef id2
#undef id3

}

void Acoo_free()
{
  free1int(Acoo->row);
  free1int(Acoo->col);
  free1complex(Acoo->val);
  free(Acoo);
}

#define USE_COMM_WORLD -987654


void Acoo_solve(emf_t *emf)
{
  int error = 0;  
  ZMUMPS_STRUC_C id;

  if(iproc == 0) {
    printf("====================================================\n");
    printf("[Acoo_solve] Enter MUMPS solver.\n");
  }

  /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
  id.comm_fortran = USE_COMM_WORLD;
  id.par  = 1;   /* host assembles matrix */
  id.sym  = 0;   /* unsymmetric */
  id.job  = -1;  /* initialize */
  zmumps_c(&id);
  if(iproc == 0) {
    printf("[Acoo_solve] MUMPS initialized (job = -1).\n");
    printf("             sym = %d, par = %d, comm_fortran = %d\n",
           id.sym, id.par, id.comm_fortran);
  }
  
  /* Define the problem on the host */
  if (iproc == 0) {
    id.n   = Acoo->nrow;
    id.nnz = (MUMPS_INT8) Acoo->nnz;
    id.irn = Acoo->row;
    id.jcn = Acoo->col;

    /* the row and column indices called by MUMPS are starting from 1 in Fortran */
    for (int k = 0; k < (int)id.nnz; k++) {
      id.irn[k]++;
      id.jcn[k]++;
    }

    id.a   = (mumps_double_complex*) Acoo->val;
    id.rhs = (mumps_double_complex*) emf->be;

    printf("[Acoo_solve] Problem assembled on host (iproc = %d).\n", iproc);
    printf("             Matrix size n   = %d\n", (int)id.n);
    printf("             Number of nnz   = %lld\n", (long long)id.nnz);
    printf("             Pointer A       = %p\n", (void*)id.a);
    printf("             Pointer RHS     = %p\n", (void*)id.rhs);

    int nprint = (id.nnz < 10) ? (int)id.nnz : 10;
    printf("             First %d (irn, jcn, val) entries:\n", nprint);
    for (int k = 0; k < nprint; k++) {
      printf("               k=%d  (%d,%d)  val=(%e,%e)\n",
             k, id.irn[k], id.jcn[k],
             id.a[k].r, id.a[k].i);
    }
  }

#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */

  id.ICNTL(1) = 6;  /* error messages unit (6 通常是 stdout)       */
  id.ICNTL(2) = 6;  /* diagnostic printing unit                     */
  id.ICNTL(3) = 6;  /* global info                                   */
  id.ICNTL(4) = 3;  /* print level: 0–3                */

  if (iproc == 0) {
    printf("[Acoo_solve] ICNTL set to verbose: (1,2,3,4) = (%d,%d,%d,%d)\n",
           id.ICNTL(1), id.ICNTL(2), id.ICNTL(3), id.ICNTL(4));
  }

  /* -------- Phase 1: ANALYSIS -------- */
  if (iproc == 0) printf("[Acoo_solve] Starting ANALYSIS phase (job=1)...\n");
  id.job = 1;
  zmumps_c(&id);

  if (iproc == 0) {
    printf("[Acoo_solve] ANALYSIS finished. INFOG(1)=%d, INFOG(2)=%d\n",
           id.infog[0], id.infog[1]);
  }
  if (id.infog[0] < 0) {
    printf(" (PROC %d) ERROR during ANALYSIS: INFOG(1)=%d INFOG(2)=%d\n",
           iproc, id.infog[0], id.infog[1]);
    error = 1;
  }

  /* -------- Phase 2: FACTORIZATION -------- */
  if (!error) {
    if (iproc == 0) printf("[Acoo_solve] Starting FACTORIZATION phase (job=2)...\n");
    id.job = 2;
    zmumps_c(&id);

    if (iproc == 0) {
      printf("[Acoo_solve] FACTORIZATION finished. INFOG(1)=%d, INFOG(2)=%d\n",
             id.infog[0], id.infog[1]);
    }
    if (id.infog[0] < 0) {
      printf(" (PROC %d) ERROR during FACTORIZATION: INFOG(1)=%d INFOG(2)=%d\n",
             iproc, id.infog[0], id.infog[1]);
      error = 1;
    }
  }

  /* -------- Phase 3: SOLVE -------- */
  if (!error) {
    if (iproc == 0) printf("[Acoo_solve] Starting SOLVE phase (job=3)...\n");
    id.job = 3;
    zmumps_c(&id);

    if (iproc == 0) {
      printf("[Acoo_solve] SOLVE finished. INFOG(1)=%d, INFOG(2)=%d\n",
             id.infog[0], id.infog[1]);

      printf("[Acoo_solve] Some INFOG values after solve:\n");
      for (int i = 0; i < 10; i++) {
        printf("             INFOG(%2d) = %d\n", i+1, id.infog[i]);
      }
      printf("[Acoo_solve] Some RINFOG values after solve:\n");
      for (int i = 0; i < 10; i++) {
        printf("             RINFOG(%2d) = %e\n", i+1, id.rinfog[i]);
      }

      printf("[Acoo_solve] MUMPS solve done, solution stored in RHS.\n");
    }
    if (id.infog[0] < 0) {
      printf(" (PROC %d) ERROR during SOLVE: INFOG(1)=%d INFOG(2)=%d\n",
             iproc, id.infog[0], id.infog[1]);
      error = 1;
    }
  }

  /* Terminate instance. */
  id.job = -2;
  zmumps_c(&id);
  if (iproc == 0) {
    printf("[Acoo_solve] MUMPS terminated (job=-2).\n");
    printf("====================================================\n");
  }
}


