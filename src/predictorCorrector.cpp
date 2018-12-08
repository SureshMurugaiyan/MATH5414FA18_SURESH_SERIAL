#include <iostream>
#include "input.h"
void eulerPredictor(double* ux, double *uy,double delT,
                    double* Cnx, double *Cny,
                    double* Dnx,double* Dny,
                    int col, int row,double vol){
for(int i = 1; i<(row-1); ++i){
    for(int j =1; j<(col-1); ++j){
        ux[i*col+j] = ux[i*col+j]+(delT/vol)*(-Cnx[i*col+j]+(Dnx[i*col+j])/Re);
        uy[i*col+j] = uy[i*col+j]+(delT/vol)*(-Cny[i*col+j]+(Dny[i*col+j])/Re);
    }
}

}


void adamPredictor(double* ux, double *uy,double delT,
                      double* Cnx, double *Cny,
                      double* Dnx,double* Dny,
                      double* CnxOld, double *CnyOld,
                      double* DnxOld,double* DnyOld,
                      int col, int row,double vol){
for(int i = 1; i<(row-1); ++i){
    for(int j =1; j<(col-1); ++j){
        ux[i*col+j] = ux[i*col+j]+(delT/vol)*(
                      1.5*(-Cnx[i*col+j]+(Dnx[i*col+j])/Re)
                     -0.5*(-CnxOld[i*col+j]+(DnxOld[i*col+j])/Re)
                      );
        uy[i*col+j] = uy[i*col+j]+(delT/vol)*(
                      1.5*(-Cny[i*col+j]+(Dny[i*col+j])/Re)
                     -0.5*(-CnyOld[i*col+j]+(DnyOld[i*col+j])/Re)
                      );
    }
}
}


void corrector(double* ux, double* uy,
               double* gradxP,double* gradyP,
               double dt, int col, int row){

for(int i = 1; i<(row-1); ++i){
    for(int j =1; j<(col-1); ++j){
        ux[i*col+j] = ux[i*col+j]-dt*gradxP[i*col+j];
        uy[i*col+j] = uy[i*col+j]-dt*gradyP[i*col+j];
    }
}

}
