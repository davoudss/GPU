#include <mex.h>
#include <iostream>
#include <string.h>

void fftshiftC(double* xr, double* xi, int N)
{
  double temp;
  for(int k=0;k<N/2;k++){
      temp = xr[k];
      xr[k] = xr[k+N/2];
      xr[k+N/2] = temp;

      temp = xi[k];
      xi[k] = xi[k+N/2];
      xi[k+N/2] = temp;
  }
}

void DisplayMatrix(char           *Name, 
		   double       *Data_r,
		   double       *Data_i, 
		   int                M, 
		   int                N
		   )
{
  int m, n;
  mexPrintf("%s = \n", Name);
  for(m = 0; m < M; m++, mexPrintf("\n"))
    for(n = 0; n < N; n++)
      mexPrintf("%+8.4f%+8.4fi ", Data_r[m + M*n],Data_i[m+M*n]);
}


// void callFFTSHIFT(double            *data_r, 
// 		  double            *data_i,
// 		  size_t                  K, 
// 		  const mwSize           *N, 
// 		  double            **sol_r, 
// 		  double            **sol_i
// 		  )
// {
//   mxArray *F,*f;
//   mxArray *LHS[1];
//   int N1 = N[0];
//   int N2 = N[1];
  

//   f = mxCreateNumericArray(K,N,mxDOUBLE_CLASS,mxCOMPLEX);
//   memcpy(mxGetPr(f),data_r,sizeof(double)*N1*N2);
//   memcpy(mxGetPi(f),data_i,sizeof(double)*N1*N2);
//   mexCallMATLABWithTrap(1,LHS,1,&f,"fftshift");

//   F = LHS[0];
  
//   //  DisplayMatrix("F",mxGetPr(F),mxGetPi(F),M,N);

//   mxDestroyArray(f);
//   mxDestroyArray(F);


// }


void callFFTN(double          *data, 
	      size_t              K,
	      const mwSize       *N,
	      double        **sol_r,
	      double        **sol_i
	      )
{
  mxArray *F,*f;
  mxArray *LHS[1];
  int N1 = N[0];
  int N2 = N[1];

  f = mxCreateNumericArray(K,N,mxDOUBLE_CLASS,mxCOMPLEX);
  memcpy(mxGetPr(f),data,sizeof(double)*N1*N2);
  mexCallMATLABWithTrap(1,LHS,1,&f,"fftn");

  F = LHS[0];

  *sol_r = (double *)mxGetData(F);
  *sol_i = (double *)mxGetImagData(F);
  
  DisplayMatrix("F",mxGetPr(F),mxGetPi(F),N1,N2);

  mxDestroyArray(f);
  //  mxDestroyArray(F);  
}


void mexFunction(int               nlhs, 
		 mxArray         **plhs,
		 int               nrhs, 
		 const mxArray   **prhs
		 )
{

#define f_IN  prhs[0]
#define f_OUT plhs[0]

  if(nrhs<1 || nrhs>1){
    mexErrMsgTxt("Wrong number of input arguments");
    return;
  }

  double *Fr,*Fi;  
  size_t K = mxGetNumberOfDimensions(f_IN);
  const mwSize *N = mxGetDimensions(f_IN);

  Fr = (double*) calloc(N[0],sizeof(Fr));
  Fi = (double*) calloc(N[1],sizeof(Fi));

  callFFTN(mxGetPr(f_IN),K,N,&Fr,&Fi);
  DisplayMatrix("SOL",Fr,Fi,N[0],N[1]);

  fftshiftC(Fr,Fi,N[1]);
  DisplayMatrix("SOL",Fr,Fi,N[0],N[1]);


  return;
}
