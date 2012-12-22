#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>


float norm (float* x, int N);
double normD (double* x, int N);
double normw(fftw_complex* x,int N);
double normC(double *x_r, double *x_i, int N);
void fftshift(fftw_complex* in, fftw_complex* out,int N);
void fftshiftn(fftw_complex* in, fftw_complex* out,int N1,int N2, int N3,int dim);

/***********************************************************************/
float norm(float* x,int N)
{
  float sum=0.;
  for(int s=0;s<N;s++)
    sum += (float) x[s]*x[s];

  return sqrt(sum);
}

/***********************************************************************/
double normD(double* x,int N)
{
  double sum=0.;
  for(int s=0;s<N;s++)
    sum += (double) x[s]*x[s];

  return sqrt(sum);
}

/***********************************************************************/
double normC(double *x_r, double *x_i, int N)
{
  double sum=0;
  for(int s=0;s<N;s++){
    sum += (double) x_r[s]*x_r[s];
    sum += (double) x_i[s]*x_i[s];
  }
  
  return sqrt(sum);
}

/***********************************************************************/
double normw(fftw_complex* x,int N)
{
  double sum=0;
  for(int s=0;s<N;s++){
    sum += (double) x[s][0]*x[s][0];
    sum += (double) x[s][1]*x[s][1];
  }

  return sqrt(sum);
}

/***********************************************************************/
void fftshift(fftw_complex* x, fftw_complex* y,int N)
{
  fftw_complex *temp;

  int M = (int) ceil((float) N/2.);
  int e = (int) fabs(N-2.*M);
  
  temp    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

  if(e>0){
    for (int k=0;k<=M-2;k++){
      temp[k+M-e][0] = x[k][0];
      temp[k+M-e][1] = x[k][1];
      
      temp[k][0] = x[M+k][0];
      temp[k][1] = x[M+k][1]; 
    }

    temp[N-1][0] = x[M-1][0];
    temp[N-1][1] = x[M-1][1];
  }
  else{
    for (int k=0;k<=M-1;k++){
      temp[k+M-e][0] = x[k][0];
      temp[k+M-e][1] = x[k][1];
      
      temp[k][0] = x[M+k][0];
      temp[k][1] = x[M+k][1];
    }
  }

  for (int k=0;k<N;k++){
    y[k][0]=temp[k][0];
    y[k][1]=temp[k][1];
  }
}



/***********************************************************************/
void fftshiftn(fftw_complex* x, fftw_complex* y,int N,int dim)
{

  int f1=0,f2=0,f3=0;
  int e1=0,e2=0,e3=0;
  int M=N/2;

  if(dim==1){
    f3=1;
    e2=1;
    e3=1;
  }
  else if(dim==2){
    f2=1;
    e2=1;
    e1=1;
  }
  else{
    f1=1;
    e1=1;
    e3=1;
  }
  
  if (x==y){
    double temp;
    
    for (int n1=0;n1<M+e2*M;n1++)
      for (int n2=0;n2<M+e3*M;n2++)
	for (int n3=0;n3<M+e1*M;n3++)
	{
	  int ind1 = (n1+f1*M)*N*N+(n2+f2*M)*N+(n3+f3*M);
	  int ind2 = n1*N*N+n2*N+n3;
	  
	  temp=x[ind1][0];
	  y[ind1][0] = x[ind2][0];
	  y[ind2][0] = temp;
	  
	  temp=x[ind1][1];
	    y[ind1][1] = x[ind2][1];
	    y[ind2][1] = temp;
	}
  }
  else{   
    
    
    for (int n1=0;n1<M+e2*M;n1++)
      for (int n2=0;n2<M+e3*M;n2++)
	for (int n3=0;n3<M+e1*M;n3++)
	{
	  int ind1 = (n1+f1*M)*N*N+(n2+f2*M)*N+(n3+f3*M);
	  int ind2 = n1*N*N+n2*N+n3;
	  
	  y[ind1][0] = x[ind2][0];
	  y[ind2][0] = x[ind1][0];
	  
	  y[ind1][1] = x[ind2][1];
	  y[ind2][1] = x[ind1][1];
	}
  }
  
}

