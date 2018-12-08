/*---------------------*- C++ 2D Incompressible FLow -*-----------------------*
|  Solves the  2D incompressible Fluid Flow in 2D geometry                    |
|  User Input is input.h File                                                 |
|  This is the main file of the solver                                        |
*-----------------------------------------------------------------------------*/
#include <stdlib.h>
#include<stdio.h>
#include "input.h"
#include "inputSerial.h"

/*----------------------------------------------------------------------------*
|                    Function Declarations                                    |
*----------------------------------------------------------------------------*/
void solveSerial();
void mesh2Dsquare(double* XX, double *YY, double *ZZ,
                  int nxcM, int nycM,int nxM, int nyM, int ncM);

/*----------------------------------------------------------------------------*
|                      Main Function                                          |
*----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
double X[nx*ny];
double Y[nx*ny];
double Z[nx*ny];

mesh2Dsquare(X,Y,Z,nxc,nyc,nx,ny,nxc*nyc);

solveSerial();  

return 0;
}
// * * * * * * * * * * END  OF PROGRAM * * * * * * * * * * * * * * * * * * * //
