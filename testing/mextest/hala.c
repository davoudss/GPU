#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray **plhs,int nrhs, const mxArray **prhs)
{
    #define B_OUT plhs[0]
    #define A_IN prhs[0]
    #define P_IN prhs[1]
            
    double *B, *A, sum = 0;
    int M, N, m, n, p;
    
    if (nrhs < 1 || nrhs >2)
        mexErrMsgTxt("Wrong");
    if (nlhs > 1)
        mexErrMsgTxt("Too many");
    
    M = mxGetM(A_IN);
    N = mxGetN(A_IN);
    A = mxGetPr(A_IN);
    
    if (nrhs == 1)
        p = 1;
    else
        p = mxGetScalar(P_IN);
    
    B_OUT = mxCreateDoubleMatrix(1, 2, mxREAL);
    B = mxGetPr(B_OUT);
        
    for (n=0;n<N;n++)
        for (m=0;m<M;m++)
                sum += A[m + M*n];
    B[0] = sum;
    
    B[1] = pow(p,2);
       
    return;   
}