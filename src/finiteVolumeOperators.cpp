#include <iostream>
#include <math.h>
using namespace std;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Divergence of a Vector with variable coefficient                         !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void Div(double* Dn, double* Phi, double* U, double* V, int row, int col,double delX,double delY){


for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
      int k    = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ee  = 0.5*(UE*PhiE+UP*PhiP);
   double Ew  = 0.5*(UW*PhiW+UP*PhiP);
   double Fn  = 0.5*(VN*PhiN+VP*PhiP);
   double Fs  = 0.5*(VS*PhiS+VP*PhiP);
   Dn[k]      = delX*(Fn-Fs)+delY*(Ee-Ew);
      }
   }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Divergence of a Vector with No coefficient                               !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void Divergence(double* Dn, double* U, double* V,int row, int col, double delX, double delY){

for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){
    int k = i*col+j;
   double UP = U[k];
   double UE = U[k+1];
   double UW = U[k-1];
   double UN = U[k-col];
   double US = U[k+col];

   double VP = V[k];
   double VE = V[k+1];
   double VW = V[k-1];
   double VN = V[k-col];
   double VS = V[k+col];

   double Ue = 0.5*(UE+UP);
   double Uw = 0.5*(UW+UP);
   double Un = 0.5*(UN+UP);
   double Us = 0.5*(US+UP);

   double Ve = 0.5*(VE+VP);
   double Vw = 0.5*(VW+VP);
   double Vn = 0.5*(VN+VP);
   double Vs = 0.5*(VS+VP);

  Dn[k] = (Ue-Uw)*delY+(Vn-Vs)*delX;
   }
 }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Laplacian of a Scalar                                                    !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void Laplacian(double* Ln, double *Phi, int row, int col, double delX, double delY){
for(int i = 1; i<(row-1); i++){
 for(int j =1; j<(col-1); j++){
   int k = i*col+j;
   double PhiP = Phi[k];
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];

   double Ee  = (PhiE-PhiP)/delX;
   double Ew  = (PhiP-PhiW)/delX;
   double Fn  = (PhiN-PhiP)/delY;
   double Fs  = (PhiP-PhiS)/delY;
   Ln[k]      = delX*(Fn-Fs)+delY*(Ee-Ew);
     }
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Gradient Computation using Finite Volume                                 !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void gradient(double* gradxPhi,double* gradyPhi,
              double* vfacePhi,double* hfacePhi,
			        int cellRow, int cellCol, double delX, double delY,
              int vfaceCol, int hfaceCol){
for(int i = 1; i<(cellRow-1); ++i){
 for(int j =1; j<(cellCol-1); ++j){


   double Phie = vfacePhi[i*vfaceCol+(j+1)];
   double Phiw = vfacePhi[i*vfaceCol+j];
   double Phin = hfacePhi[i*hfaceCol+j];
   double Phis = hfacePhi[(i+1)*hfaceCol+j];

   gradxPhi[i*cellCol+j] = (Phie-Phiw)/delX;
   gradyPhi[i*cellCol+j] = (Phin-Phis)/delY;
    }
  }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Gradient Computation using Finite Volume                                 !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void gradientFV(double* gradxPhi,double* gradyPhi,double* Phi,
                        int row, int col, double delX, double delY){
for(int i = 1; i<(row-1); ++i){
 for(int j =1; j<(col-1); ++j){

   int k = i*col+j;
   double PhiE = Phi[k+1];
   double PhiW = Phi[k-1];
   double PhiN = Phi[k-col];
   double PhiS = Phi[k+col];
   double PhiP = Phi[k];

   double Phie = 0.5*(PhiE + PhiP);
   double Phiw = 0.5*(PhiW + PhiP);
   double Phin = 0.5*(PhiN + PhiP);
   double Phis = 0.5*(PhiS + PhiP);

   gradxPhi[k] = (Phie-Phiw)/delX;
   gradyPhi[k] = (Phin-Phis)/delY;
    }
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Interpolation of Vertex Values                                           !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

void vertexInterpolate(double* vertPhi, double* cellPhi, int vertRow, int vertCol, int cellCol){
// Inner Points
for(int i = 1; i<(vertRow-1); ++i){
 for(int j =1; j<(vertCol-1); ++j){

double Ta = cellPhi[(i-1)*cellCol+(j-1)]; // NW
double Tb = cellPhi[(i-1)*cellCol+j];    // NE
double Tc = cellPhi[(i*cellCol)+(j-1)];      // SW
double Td = cellPhi[i*cellCol+j];      //SE

vertPhi[i*vertCol+j]=0.25*(Ta+Tb+Tc+Td);
   }
 }
//Boundary points - Set it to Boundary Cell points
//North Boundary
    for (int i = 0; i<vertCol; i++){
        vertPhi[i]= cellPhi[i+vertCol];
    }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Interpolation of Vertical FaceCenter Interpolation                       !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void vfaceInterpolate(double* vfacePhi, double* vertPhi,
                      double* cellPhi, int vfaceRow, int vfaceCol,
                      int vertCol, int cellCol){

for(int i = 1; i<(vfaceRow-1); ++i){
 for(int j =1; j<(vfaceCol-1); ++j){

double Ta = cellPhi[i*cellCol+(j-1)]; // West cell
double Tb = cellPhi[i*cellCol+j];     // East cell
double Tc = vertPhi[i*vertCol+j];     // North Vertex
double Td = vertPhi[(i+1)*vertCol+j]; //South Vertex

vfacePhi[i*vfaceCol+j]=0.25*(Ta+Tb+Tc+Td);
  }
}

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
// Interpolation of Horizontal FaceCenter Interpolation                     !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
void hfaceInterpolate(double* hfacePhi, double* vertPhi,
                      double* cellPhi, int hfaceRow, int hfaceCol,
                      int vertCol, int cellCol){
for(int i = 1; i<(hfaceRow-1); ++i){
 for(int j =1; j<(hfaceCol-1); ++j){

double Ta = vertPhi[i*vertCol+j]; // West vertex
double Tb = vertPhi[i*vertCol+(j+1)];     // East vertex
double Tc = cellPhi[(i-1)*cellCol+j];     // North cell
double Td = cellPhi[i*cellCol+j]; //South cell

hfacePhi[i*hfaceCol+j]=0.25*(Ta+Tb+Tc+Td);
 }
}
}


