#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "dfftlib.h"
#include "dmathlib.h"
#include "diolib.h"
#include "fgglib.h"


int main()
{
  
  int M,N,Mr,P,cnt,nt;
  double h,l,Tau,nrm1,sum,nrm2,L,a,b,ktilde,sgn=1.;
  double *temp,*xj2,*xj1,*xj,*fj,*Ftau,*ftau,*Fexact,*Fapprox,gtau,fjs;
  double start, end;


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
  M = 16*1024;
  N = 16*128;
  Mr = 2*M;
  
  h = (double) 2.*pi/(double)Mr;
  Tau = (double) 12.f/(double)(M*M);

  P = 12;

  xj2     = (double*) malloc(sizeof(double)*N);
  xj1     = (double*) malloc(sizeof(double)*N);
  xj      = (double*) malloc(sizeof(double)*N);
  fj      = (double*) malloc(sizeof(double)*N);

  Ftau    = (double*) calloc(2*Mr,sizeof(double));
  ftau    = (double*) calloc(2*Mr,sizeof(double));

  Fexact  = (double*) malloc(sizeof(double)*2*M);
  Fapprox = (double*) calloc(2*M,sizeof(double));
  temp  = (double*) malloc(sizeof(double)*2*M);

#pragma omp parallel for shared(xj2,xj1,xj,fj)
  for(int s=0;s<N;s++){
    //    srand(s);
    //    xj[s] = rand() % L;
    //    fj[s] = 10*xj[s];
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
      Fexact[2*cnt]   += fj[s]*cos(-ktilde*xj2[s]);
      Fexact[2*cnt+1] += fj[s]*sin(-ktilde*xj2[s]);
    }
    cnt++;
  }

  //before fft ruins it
#pragma omp parallel for shared(temp,Fexact)
  for (int k=0;k<2*M;k++)
    temp[k] = Fexact[k];

  //  printC(Fexact,0,2*M);
  nrm1 = normD(Fexact,2*M);

   //printf("\nFexact: %e",nrm1);
  // printFile(Fexact,2*M,1);

  //Find approximate sum
  //--------------------------------------------------
  //Gridding on an oversampled mesh
#pragma omp parallel for shared(ftau,fj)
  for (int s=0;s<N;s++){
    int m1 = (int) round((xj[s])/h);
    for (int m=m1-P;m<=m1+P;m++){
      gtau = fggD(xj[s],m,Tau,Mr);
      if(m<0)
	ftau[2*(m+Mr)] += fj[s]*gtau;
      else if (m>=Mr)
	ftau[2*(m-Mr)] += fj[s]*gtau;
      else
	ftau[2*(m)]    += fj[s]*gtau;
    }
  }

  //Note that ftau filled with garbage values after fft and N is power of two
  fftD(ftau,Ftau,Mr,sgn);
  fftshiftD(Ftau,2*Mr);


  //down sampling points
  cnt = 0;
#pragma omp parallel for shared(Fapprox)
  for (int s=M/2;s<M+M/2;s++){
    double k = s-M;
    ktilde = k*pi/L;

    double reFtau    =  sqrt(pi/Tau)*exp(k*k*Tau)*Ftau[2*s]/Mr;
    double imFtau    =  sqrt(pi/Tau)*exp(k*k*Tau)*Ftau[2*s+1]/Mr;
    double reScaling =  cos(ktilde*(b+a)/2.);
    double imScaling = -sin(ktilde*(b+a)/2.);

    Fapprox[2*cnt]   = reFtau*reScaling-imFtau*imScaling;
    Fapprox[2*cnt+1] = reFtau*imScaling+imFtau*reScaling;
    cnt++;
    //    printf("%f\n", sqrt(PI/Tau)*exp(k*k*Tau)/Mr);
  }


  //find the difference
#pragma omp parallel for shared(Fapprox,temp)
  for (int k=0;k<2*M;k++)
    Fapprox[k] = Fapprox[k] - temp[k];

  nrm2 = normD(Fapprox,2*M);
  // printf("\nFapprox: %e",nrm);

  // printf("\n");
  // printFile(Fapprox,2*M,2);

  printf("\n Error: %e\n",nrm2/nrm1);

  if(omp_get_thread_num()==0)
    end = omp_get_wtime();

  printf("Duration: %f\n",end-start);

  free(xj);
  free(fj);
}


