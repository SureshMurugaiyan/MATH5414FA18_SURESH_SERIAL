#include <iostream>
#include <math.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Function Declarations                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void storeOldValue(double *phinew, double *phiOld,int totCell);
void L2norm(double *Phinew, double *Phiold,double *L2Phi,int  totCell);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Main Function-->Poisson Solver for Pressure Finite Volume Solver         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void PoissonPressure(double* Phi, int row, int col,
             double delX,double delY,double* source,
             int totCell){
double lam = 1;
double PhiOld[totCell];
double normPhi = 1.0;
double L2Phi = 1.0;
double L2oPhi= 1.0;
double normTarget  = 1.0e-5;

int itr = 0;
int stop = 0;
while (stop==0){
itr++;
storeOldValue(Phi,PhiOld,totCell);

for(int i=1; i<(row-1); i++){
 for(int j=1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double AP   = (-2*delY/delX)-(2*delX/delY);
   double AS   = (delX/delY);
   double AW   = (delY/delX);
   double AE   = (delY/delX);
   double AN   = (delX/delY);

   double  R     = source[k]- AP*PhiP-AE*PhiE-AW*PhiW-AN*PhiN-AS*PhiS;
   double delPhi = R/AP;
   Phi[k] = Phi[k]+lam*delPhi;
  }
 }
if(itr==1){
  //L2norm(Phi,PhiOld, &L2oPhi,totCell);
 // L2norm(Phi,PhiOld, &L2Phi,totCell);

} else {
 // L2norm(Phi,PhiOld, &L2Phi,totCell);
}

int selectNorm = 1;        // choose 0 for normalized norm or 1 for actual norm
if(selectNorm==0){
normPhi = L2Phi/L2oPhi;    // Normalized norm wrt to initial correction
}else{
normPhi = L2Phi;           // Actual norm without normalizing with initial correction
}

//if(normPhi<(1e-5)){stop=1;}
if(itr>500){stop=1;}
  }

}
