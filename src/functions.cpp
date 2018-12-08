#include <iostream>
#include <stdio.h>
#include <math.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Initialization of arrays for storing Primitive Variables
void InitializeField(double *phi, int row, int col){
for(int i = 0; i<row; ++i){
 for(int j =0; j<col; ++j){
   phi[i*col+j]=0.0;
    }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void updateBoundaryCondition(double* ux, double *uy, double *p,
                             int col,int totCell){
// updateBoudaryCondition for serial mode
    // North Boundary
    for (int i = 0; i<col; i++)
    {
        ux[i]= 1;
        uy[i]= 0;
        p[i]  = p[i+col];
    }
    // South Boundary
    for (int i = totCell-col; i<totCell; i++)
    {
        ux[i]= 0;
        uy[i]= 0;
        p[i]= p[i-col];
    }
    // West Boundary - Left end
    for (int i = 0; i<totCell; i=(i+col))
    {
        ux[i]= 0;
        uy[i]= 0;
        p[i] = p[i+1];
    }
    // East Boundary - Right end
    for (int i = col-1; i<totCell; i=(i+col))
    {
        ux[i]=0;
        uy[i]=0;
        p[i]  = p[i-1];
    }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void refValueUpdate(double* Phi, int row, int col, int refCell){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
   Phi[i*col+j]=Phi[i*col+j]-Phi[refCell];
    }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 void storeOldValue(double *phinew, double *phiOld,int totCell){
   for(int i =0; i<totCell; i++){
     phiOld[i]=phinew[i];
   }
 }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void updateCn(double* Cn,double dt, int col,int row){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
    Cn[i*col+j] = Cn[i*col+j]/dt;
  }
 }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void L2norm(double *Phinew, double *Phiold,double *L2Phi,int  totCell){
*L2Phi = 0;
  for(int i = 0; i<totCell;i++){
    *L2Phi= *L2Phi+pow((Phiold[i]-Phinew[i]),2);
   }
   *L2Phi=sqrt(*L2Phi/totCell);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void normL2(double *Phi1new, double *Phi1old,
	    double *Phi2new, double *Phi2old,
            double *Phi3new, double *Phi3old,
            double *L2Phi,   int  totCell){
    for(int j = 0; j<totCell;j++){
    L2Phi[0]= L2Phi[0]+pow((Phi1old[j]-Phi1new[j]),2);
    L2Phi[1]= L2Phi[1]+pow((Phi2old[j]-Phi2new[j]),2);
    L2Phi[2]= L2Phi[2]+pow((Phi3old[j]-Phi3new[j]),2);
   }
  L2Phi[0]=(L2Phi[0])/totCell;
  L2Phi[1]=(L2Phi[1])/totCell;
  L2Phi[2]=(L2Phi[2])/totCell;
// square root is not performed here.. perform it when you print it
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
