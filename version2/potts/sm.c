/*==========================================================
 * blocksum.c - example in MATLAB External Interfaces
 *
 *  rowsums a (mnxp)matrix in blocks of m 
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = blocksum( inMatrix, sample space size(m),nuobs n)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"

double rowsum(double *y, double *z,  int p, int m, int obs, int n)
{
    int i,j,k, ind;
    
    
    for(i = 0; i<p;i++){
        ind = i*n + obs - 1;
        z[ind] = 0;
        for(j = 0;j<m;j++){
            k =  i*n*m + (obs -1)*m + j;/*i*m + j + (obs-1)*m*p;*/
            z[ind] = z[ind]+ y[k];
          /* printf("%f\n",z[i]);*/
        }   
    }
 
}



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double blsize;                  /* input scalar */
    double nuobs;
    double *inMatrix;               /* 1xN input matrix */ /*size_t ncols;*/
    
    int nrows,m, ncols,i;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments 
    
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    /* check that number of rows in second input argument is 1 */
  /*  if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }*/
    
   

    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    
    ncols = mxGetN(prhs[0]);
    m = (int)mxGetScalar(prhs[1]);
    nrows = (int)mxGetScalar(prhs[2]);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
    outMatrix = mxGetPr(plhs[0]);

    /* get a pointer to the real data in the output matrix */
     
    int ind;
    ind = nrows + 1;
    /* call the computational routine */
    for( i = 1; i < ind;i++){
     rowsum(inMatrix,outMatrix,(int)ncols,(int) m,(int)i,(int)nrows);
    }
     
}
