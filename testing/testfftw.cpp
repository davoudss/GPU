#include <stdio.h>
#include <iostream>
//#include <complex>
#include <fftw3.h>

int main(){

  int N=4;
  fftw_complex *in, *out;
  fftw_plan p;


  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  
  in[0][0] = 1;
  in[0][1] = 2;
  in[1][0] = -1;
  in[1][1] = 1;
  in[2][0] = 0;
  in[2][1] = -2;
  in[3][0] = 2;
  in[3][1] = 3;

  for (int k=0;k<4;k++)
    printf("(%f,%f)\n",in[k][0],in[k][1]);

  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);
  for (int k=0;k<4;k++)
    printf("(%f,%f)\n",out[k][0],out[k][1]);


  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return 0;

}
