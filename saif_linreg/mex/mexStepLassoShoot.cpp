#include "mex.h"
#include <cmath>
#include "blas.h"
#include <stdio.h>
#include <string.h>
//#include <matrix.h>

typedef ptrdiff_t intt;


///////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
  double *XX2T, *Xy2, *ibeta, *dp, *dlambda, *dmaxStep;   
  XX2T = (double*)mxGetPr(prhs[0]);
  Xy2 = (double*)mxGetPr(prhs[1]);
  ibeta = (double*)mxGetPr(prhs[2]);
  dp = (double*)mxGetPr(prhs[3]);
  dlambda = (double*)mxGetPr(prhs[4]);
  dmaxStep = (double*)mxGetPr(prhs[5]);
  
 int p = *dp;
 //int dlambda 
 double lambda = *dlambda;
 int maxStep = *dmaxStep;
 double SumZ = 0;
 double *xxpt, *bpt, *xypt; 
 double xxjj;
 
 
 plhs[0]  = mxCreateDoubleMatrix(p, 1, mxREAL);
 double * obeta = (double *)mxGetPr(plhs[0]);

 
 memcpy(obeta, ibeta, p*sizeof(double));

 for(int m = 0; m< maxStep; m++)
 {
    for(int j = 0; j<p; j++)
    {
       SumZ = 0;
       xxjj = XX2T[j + j*p];
       xxpt = &XX2T[0 + j*p];
       bpt = &obeta[0];  
  
       for(int i = 0; i<p; i++) { SumZ += ((*xxpt)*(*bpt)); xxpt++; bpt++;}
       
       SumZ = SumZ - xxjj*obeta[j] - Xy2[j];
       
       if(SumZ > lambda)
        { obeta[j] = (lambda - SumZ)/xxjj; }
       else if( SumZ < -lambda )
        { obeta[j] = (-lambda - SumZ)/xxjj; }
       else 
        { obeta[j] = 0; }
    }
 }
 
// plhs[0] = ibeta;
}






