#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>
//#include "dfftlib.h"
#include "dmathlib.h"
#include "diolib.h"
#include "fgglib.h"
#include <fftw3.h>


int main()
{
  
  int M,N,Mr,P,cnt,nt;
  double h,l,Tau,nrm1,sum,nrm2,L,a,b,ktilde,sgn=1.;
  double *xj2,*xj1,*xj,*fj,gtau,fjs;
  fftw_complex *temp,*Ftau,*ftau,*Fexact,*Fapprox;
  double start, end;


  //FFT SETTING
  fftw_plan p;


  printf("Number of threads: ");
  //  scanf("%d",&nt);
  std::cin>>nt;
  omp_set_num_threads(nt);
  nt = omp_get_max_threads();
  printf("\nRunning on %d threads",nt);

  if(omp_get_thread_num()==0)
    start = omp_get_wtime();

  a = 2.;
  b = 7.;
  L = (b-a)/2.;
  M = 16;
  N = 8;
  Mr = 2*M;
  
  h = (double) 2.*pi/(double)Mr;
  Tau = (double) 12.f/(double)(M*M);

  P = 12;

  xj2     = (double*) malloc(sizeof(double)*N);
  xj1     = (double*) malloc(sizeof(double)*N);
  xj      = (double*) malloc(sizeof(double)*N);
  fj      = (double*) malloc(sizeof(double)*N);

  Ftau    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Mr);
  ftau    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Mr);

  Fexact  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
  Fapprox = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);
  temp    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M);


#pragma omp parallel for shared(xj2,xj1,xj,fj)
  for(int s=0;s<N;s++){
    xj2[s] = (double) a+2.*L*cos((double) s+1.)*cos((double) s+1.);
    xj1[s] = (double) xj2[s]-(b+a)/2.;
    xj[s]  = (double) pi*xj1[s]/L;
    fj[s]  = (double) -1.+2.*xj1[s];
  }


  //find the exact sum
  cnt = 0;
#pragma omp parallel for shared(Fexact,fj,xj2)
  for (int k=-M/2;k<=M/2-1;k++){
    ktilde = k*pi/L;
    for (int s=0;s<N;s++){
      Fexact[cnt][0] += fj[s]*cos(-ktilde*xj2[s]);
      Fexact[cnt][1] += fj[s]*sin(-ktilde*xj2[s]);
    }
    cnt++;
  }

  nrm1 = normC(Fexact,M);

  //before fft ruins it
#pragma omp parallel for shared(temp,Fexact)
  for (int k=0;k<M;k++){
    temp[k][0] = Fexact[k][0];
    temp[k][1] = Fexact[k][1];
  }


  //Find approximate sum
  //--------------------------------------------------
  //Gridding on an oversampled mesh
#pragma omp parallel for shared(ftau,fj)
  for (int s=0;s<N;s++){
    int m1 = (int) round((xj[s])/h);
    for (int m=m1-P;m<=m1+P;m++){
      gtau = fggD(xj[s],m,Tau,Mr);
      if(m<0)
	ftau[m+Mr][0] += fj[s]*gtau;
      else if (m>=Mr)
	ftau[m-Mr][0] += fj[s]*gtau;
      else
	ftau[m][0]    += fj[s]*gtau;
    }
  }

  //Note that ftau filled with garbage values after fft and N is power of two
  //  fftD(ftau,Ftau,Mr,sgn);
  //fftshiftD(Ftau,2*Mr);

  p = fftw_plan_dft_1d(Mr,ftau,Ftau,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);
  fftshift(Ftau,Ftau,Mr);

  //down sampling points
  cnt = 0;
#pragma omp parallel for shared(Fapprox)
  for (int s=M/2;s<M+M/2;s++){
    double k = s-M;
    ktilde = k*pi/L;

    double reFtau    =  sqrt(pi/Tau)*exp(k*k*Tau)*Ftau[s][0]/Mr;
    double imFtau    =  sqrt(pi/Tau)*exp(k*k*Tau)*Ftau[s][1]/Mr;
    double reScaling =  cos(ktilde*(b+a)/2.);
    double imScaling = -sin(ktilde*(b+a)/2.);

    Fapprox[cnt][0]  = reFtau*reScaling-imFtau*imScaling;
    Fapprox[cnt][1]  = reFtau*imScaling+imFtau*reScaling;
    cnt++;
  }


  //find the difference
#pragma omp parallel for shared(Fapprox,temp)
  for (int k=0;k<2*M;k++){
    Fapprox[k][0] = Fapprox[k][0] - Fexact[k][0];
    Fapprox[k][1] = Fapprox[k][1] - Fexact[k][1];
  }

  nrm2 = normC(Fapprox,M);
  printf("\n Error: %e\n",nrm2/nrm1);

  if(omp_get_thread_num()==0)
    end = omp_get_wtime();

  printf("Duration: %f\n",end-start);

  free(xj);
  free(fj);
}


