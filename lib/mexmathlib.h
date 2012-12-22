#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>


float norm (float* x, int N);
double normD (double* x, int N);
double normC(fftw_complex* x,int N);

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
double normC(fftw_complex* x,int N)
{
  double sum=0;
  for(int s=0;s<N;s++){
    sum += (double) x[s][0]*x[s][0];
    sum += (double) x[s][1]*x[s][1];
  }

  return sqrt(sum);
}
