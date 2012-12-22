#include <iostream>
#include <fftw3.h>
#include "../lib/dmathlib.h"



int main(){
  fftw_complex *x,*y;
  int N=2;
  x   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);
  y   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*N*N);


  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      for (int k=0;k<N;k++){
	x[i*N*N+j*N+k][0] = i*N*N+j*N+k+1;
	x[i*N*N+j*N+k][1] = 0;

	std::cout << x[i*N*N+j*N+k][0] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }

  fftshiftn(x,x,N,1);
  fftshiftn(x,x,N,2);
  fftshiftn(x,x,N,3);

  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      for (int k=0;k<N;k++){
	std::cout << x[i*N*N+j*N+k][0] << " ";
      }
      std::cout << "\n";
    }
    std::cout << "\n\n";
  }


  return 0;
}
  

