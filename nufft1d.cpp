#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "dfftlib.h"
#include "dmathlib.h"
#include "diolib.h"

const double PI = 3.141592653589793;


int main()
{
  
  int M,N,Mr,P,cnt;
  double h,l,Tau,pi,nrm1,sum,nrm2,L;
  double *temp,*xj,*fj,*Ftau,*ftau,*Fexact,*Fapprox,gtau,fjs;

  L = 2.*PI;
  M = 64;
  N = 4;
  Mr = 2*M;
  
  h = (double) 2.*PI/(double)Mr;
  Tau = (double) 12./(double)(M*M);

  P = 12;

  //k = (-M/2:M/2-1)';

  xj     = (double*) malloc(sizeof(double)*N);
  fj     = (double*) malloc(sizeof(double)*N);

  Ftau   = (double*) calloc(2*Mr,sizeof(double));
  ftau   = (double*) calloc(2*Mr,sizeof(double));

  Fexact = (double*) malloc(sizeof(double)*2*M);
  Fapprox= (double*) calloc(2*M,sizeof(double));
  temp= (double*) calloc(2*M,sizeof(double));


  for(int s=0;s<N;s++){
    //    srand(s);
    //    xj[s] = rand() % L;
    //    fj[s] = 10*xj[s];

    xj[s] = (double) (s+1)*PI/N;
    fj[s] = (double) (s+1)*PI/N;
  }

  
  //find the exact sum
  cnt = 0;
  for (int k=-M/2;k<=M/2-1;k++){
    for (int s=0;s<N;s++){
      Fexact[2*cnt]   += fj[s]*cos(-k*xj[s]);
      Fexact[2*cnt+1] += fj[s]*sin(-k*xj[s]);
    }
    cnt++;
  }
 

  //fft ruins the data
  for(int k=0;k<2*M;k++)
    temp[k] = Fexact[k];


  printC(Fexact,0,M);
  nrm1 = norm(Fexact,2*M);
  // printf("\nFexact: %e",nrm);
  // printFile(Fexact,2*M,1);

  //Find approximate sum
  //--------------------------------------------------
  //Gridding on an oversampled mesh
  
  for (int s=0;s<N;s++){
    int m1 = (int) round((xj[s])/h);
    for (int m=m1-P;m<=m1+P;m++){
      gtau = exp(-(2.*PI*m/Mr-xj[s]-2.*PI)*(2.*PI*m/Mr-xj[s]-2.*PI)/(4.*Tau)) 
            +exp(-(2.*PI*m/Mr-xj[s])*(2.*PI*m/Mr-xj[s])/(4.*Tau))
            +exp(-(2.*PI*m/Mr-xj[s]+2.*PI)*(2.*PI*m/Mr-xj[s]+2.*PI)/(4.*Tau));
      if(m<0)
	ftau[2*(m+Mr)] += fj[s]*gtau;
      else if (m>=Mr)
	ftau[2*(m-Mr)] += fj[s]*gtau;
      else
	ftau[2*(m)]    += fj[s]*gtau;
    }
  }


  //Note that ftau filled with garbage values after fft and N is power of two
  double sgn = 1.0;
  fft(ftau,Ftau,Mr,sgn);
  fftshift(Ftau,2*Mr);

  //down sampling points
  cnt = 0;
  for (int s=M/2;s<M+M/2;s++){
    int k = s-M;
    Fapprox[2*cnt]   = sqrt(PI/Tau)*exp(k*k*Tau)*Ftau[2*s]/Mr;
    Fapprox[2*cnt+1] = sqrt(PI/Tau)*exp(k*k*Tau)*Ftau[2*s+1]/Mr;
    cnt++;
    //    printf("%f\n", sqrt(PI/Tau)*exp(k*k*Tau)/Mr);
  }

  printC(Fapprox,0,M);

  nrm2 = norm(Fapprox,2*M);
  printf("\napprox %e\n",nrm2);
  printf("\nexact %e\n",nrm1);
  for(int k=0;k<2*M;k++)
    Fapprox[k] = Fapprox[k]-temp[k];

  // printC(Fapprox,0,2*M);
  nrm2 = norm(Fapprox,2*M);
  printf("\ndiff: %e",nrm2);

  // printf("\n");
  // printFile(Fapprox,2*M,2);

  printf("\n Error: %e\n",nrm2/nrm1);

  free(xj);
  free(fj);
}


