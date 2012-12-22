#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray **plhs,int nrhs, const mxArray **prhs)
{
    #define B_OUT plhs[0]
    #define A_IN prhs[0]
    #define P_IN prhs[1]
            
    double *B, *A, sumr = 0, sumi = 0, *Ar, *Ai;
    int M, N, m, n, p, flag = 0;
    
    if (nrhs < 1 || nrhs >2)
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Too many input arguments");
    
    M = mxGetM(A_IN);
    N = mxGetN(A_IN);
    A = mxGetPr(A_IN);
    
    if (nrhs == 1)
        p = 1;
    else
        p = mxGetScalar(P_IN);
    
    m = 1; n = 2;
    B_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
    
    //mxSetData(B_OUT, mxMalloc(sizeof(double)*m*n)); //allocating memory
    B = mxGetPr(B_OUT);
    printf("OK");
    if (mxIsComplex(A_IN)){
        Ar = mxGetPr(A_IN);
        Ai = mxGetPi(A_IN);
        flag = 1;
    }
    
    for (n=0 ; n<N ; n++)
        for (m=0 ; m<M ; m++){
            if (flag){
                sumr += Ar[m + M*n];
                sumi += Ai[m + M*n];
            }
            else
                sumr += A[m + M*n];
        }
                
    B[0] = sumr;
    
    B[1] = sumi;
   

    return;   
}
