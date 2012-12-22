#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <omp.h>
//#include "dfftlib.h"
#include "lib/dmathlib.h"
#include "lib/diolib.h"
#include "lib/fgglib.h"
#include <fftw3.h>

#define fj(i,j,k)             fj[k+N*(j+N*i)]

#define Fexact_r(i,j,k)       Fexact_r[k+M*(j+M*i)]
#define Fexact_i(i,j,k)       Fexact_i[k+M*(j+M*i)]
#define Fapprox_r(i,j,k)      Fapprox_r[k+M*(j+M*i)]
#define Fapprox_i(i,j,k)      Fapprox_i[k+M*(j+M*i)]


int main()
{
  
  int *u,M,N,Mr,P,cnt,nt,c1,c2,c3;
  double h,Tau,nrm1,nrm2,L,a,b,utilde,vtilde, wtilde;
  double *xj2,*xj1,*xj,*yj2,*yj1,*yj,*zj2,*zj1,*zj,*fj,gtau;
  double *Fapprox_r, *Fapprox_i, *Fexact_r, *Fexact_i;
  double start, end;
  fftw_complex *Ftau,*ftau;
  fftw_plan plan;


  printf("Number of threads: ");
  //  scanf("%d",&nt);
  std::cin>>nt;
  omp_set_num_threads(nt);
  nt = omp_get_max_threads();
  printf("\nRunning on %d threads\n",nt);

  if(omp_get_thread_num()==0)
    start = omp_get_wtime();
 
  a = 2.;
  b = 7.;
  L = (b-a)/2.;
  M = 32;
  N = 4;
  Mr = 2*M;
  
  h = (double) 2.*pi/(double)Mr;
  Tau = (double) 12.f/(double)(M*M);

  P = 12;
  const size_t DATA_SIZE = N*N*N;
  const   size_t SAMPLE_SIZE = M*M*M;
  const  size_t OVERSAMPLE_SIZE = Mr*Mr*Mr;

  xj2     = (double*) malloc(sizeof(double)*N);
  xj1     = (double*) malloc(sizeof(double)*N);
  xj      = (double*) malloc(sizeof(double)*N);
  yj2     = (double*) malloc(sizeof(double)*N);
  yj1     = (double*) malloc(sizeof(double)*N);
  yj      = (double*) malloc(sizeof(double)*N);
  zj2     = (double*) malloc(sizeof(double)*N);
  zj1     = (double*) malloc(sizeof(double)*N);
  zj      = (double*) malloc(sizeof(double)*N);

  fj      = (double*) malloc(sizeof(double)*DATA_SIZE);

  Ftau    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * OVERSAMPLE_SIZE);

  ftau    = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * OVERSAMPLE_SIZE);

  Fexact_r  = (double*) malloc(sizeof(double)*SAMPLE_SIZE);
  Fexact_i  = (double*) malloc(sizeof(double)*SAMPLE_SIZE);

  Fapprox_r = (double*) calloc(SAMPLE_SIZE,sizeof(double));
  Fapprox_i = (double*) calloc(SAMPLE_SIZE,sizeof(double));

  //initialize ftau
  for (int n1=0;n1<Mr;n1++)
    for (int n2=0;n2<Mr;n2++)
      for (int n3=0;n3<Mr;n3++){
	ftau[n1*Mr*Mr+n2*M+n3][0] = 0.;
	ftau[n1*Mr*Mr+n2*M+n3][1] = 0.;
      }


  //Initialize vector of k space values
  u       = (int*)    malloc(sizeof(int)*M);

  for (int cu=0;cu<M;cu++)
    u[cu] = -M/2 + cu;

  //#pragma omp parallel for shared(xj2,xj1,xj,fj)
  for(int n=0;n<N;n++){
	//    srand(s);
	//    xj[s] = rand() % L;
	//    fj[s] = 10*xj[s];
	xj2[n] = (double) a+2.*L*cos((double) n+1.)*cos((double) n+1.);
	xj1[n] = (double) xj2[n]-(b+a)/2.;
	xj[n]  = (double) pi*xj1[n]/L;
	yj2[n] = (double) a+2.*L*cos((double) n+1.)*cos((double) n+1.);
	yj1[n] = (double) yj2[n]-(b+a)/2.;
	yj[n]  = (double) pi*yj1[n]/L;
	zj2[n] = (double) a+2.*L*cos((double) n+1.)*cos((double) n+1.);
	zj1[n] = (double) zj2[n]-(b+a)/2.;
	zj[n]  = (double) pi*zj1[n]/L;
      }


  for(int n1=0;n1<N;n1++){
    for(int n2=0;n2<N;n2++){
      for(int n3=0;n3<N;n3++){
	fj(n1,n2,n3)  = (double) -1.+2.*xj1[n1]+3.*yj1[n2]-2.*zj1[n3];
	printf("%f\t",fj(n1,n2,n3));}
      printf("\n");}
    printf("\n\n");}

  printf("%f\n",normD(xj,N));
  printf("%f\n",normD(fj,N*N*N));

  return 0;


  //find the exact sum
  cnt = 0;
  //#pragma omp parallel for shared(Fexact,fj,xj2)
  for (int cu=0;cu<M;cu++){
    utilde = u[cu]*pi/L;
    for (int cv=0;cv<M;cv++){
      vtilde = u[cv]*pi/L;
      for (int cw=0;cw<M;cw++){
	wtilde = u[cw]*pi/L;
	for (int n1=0;n1<N;n1++){
	  for (int n2=0;n2<N;n2++){
	    for (int n3=0;n3<N;n3++){
	      double val = -utilde*xj2[n1]-vtilde*xj2[n2]-wtilde*xj2[n3];
	      Fexact_r(n1,n2,n3)   += fj(n1,n2,n3)*cos(val);
	      Fexact_i(n1,n2,n3)   += fj(n1,n2,n3)*sin(val);
	    }	   
	  }
	}
      }
    }
  }
      
  //  printC(Fexact,0,2*M);
  nrm1 = normC(Fexact_r,Fexact_i,M*M*M);

   //printf("\nFexact: %e",nrm1);
  // printFile(Fexact,2*M,1);

  //Find approximate sum
  //--------------------------------------------------
  //Gridding on an oversampled mesh
  //#pragma omp parallel for shared(ftau,fj)
  for (int n1=0;n1<N;n1++){
    int cut1 = (int) round(xj[n1]/h);
    for (int n2=0;n2<N;n2++){
      int cut2 = (int) round(yj[n2]/h);
      for (int n3=0;n3<N;n3++){
	int cut3 = (int) round(zj[n3]/h);
	for (int m1=cut1-P;m1<=cut1+P;m1++){
	  for (int m2=cut2-P;m2<=cut2+P;m2++){
	    for (int m3=cut3-P;m3<=cut3+P;m3++){
	      gtau = fgg3(xj[n1],yj[n2],zj[n3],m1,m2,m3,Tau,Mr);
	      if(m1>=0 && m1<Mr)
		c1 = m1;
	      else if (m1>=Mr)
		c1 = m1-Mr;
	      else
		c1 = m1+Mr;
	      
	      if(m2>=0 && m2<Mr)
		c2 = m2;
	      else if (m2>=Mr)
		c2 = m2-Mr;
	      else
		c2 = m2+Mr;
	      
	      if(m3>=0 && m3<Mr)
		c3 = m3;
	      else if (m3>=Mr)
		c3 = m3-Mr;
	      else
		c3 = m3+Mr;
	      
	      ftau[c1*Mr*Mr+c2*Mr+c3][0]   += fj(c1,c2,c3)*gtau;	     
	    }
	  }
	}
      }
    }
  }


  //copy ftau_r & ftau_i (which is not allocated as it is zero) to ftau for fftw

  plan = fftw_plan_dft_3d(Mr,Mr,Mr,ftau,Ftau,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(plan);
  fftshiftn(Ftau,Ftau,Mr,1);
  fftshiftn(Ftau,Ftau,Mr,2);
  fftshiftn(Ftau,Ftau,Mr,3);


  //down sampling points
  cnt = 0;
  //#pragma omp parallel for shared(Fapprox_r,Fapprox_i)
  for (int n1=0;n1<M;n1++){
    utilde = u[n1]*pi/L;
    for (int n2=0;n2<M;n2++){
      vtilde = u[n2]*pi/L;
      for (int n3=0;n3<M;n3++){
	  wtilde = u[n3]*pi/L;
 
	  double mul    =  sqrt(pi/Tau)*sqrt(pi/Tau)*sqrt(pi/Tau)*exp((u[n1]*u[n1]+u[n2]*u[n2]+u[n3]*u[n3])*Tau)/(double) (Mr*Mr*Mr);
	  double reScaling =  cos((utilde+vtilde+wtilde)*(b+a)/2.);
	  double imScaling = -sin((utilde+vtilde+wtilde)*(b+a)/2.);
	  int ind = n1*M*M+n2*M+n3;
	  
	  Fapprox_r(n1,n2,n3)   = (Ftau[ind][0]*reScaling-Ftau[ind][1]*imScaling)*mul;
	  Fapprox_i(n1,n2,n3)   = (Ftau[ind][0]*imScaling+Ftau[ind][1]*reScaling)*mul;
	}
    }
  }


  //find the difference
  //#pragma omp parallel for shared(Fapprox_r,Fapprox_i)
  for (int k=0;k<M*M*M;k++){
    Fapprox_r[k] = Fapprox_r[k] - Fexact_r[k];
    Fapprox_i[k] = Fapprox_i[k] - Fexact_i[k];
  }

  nrm2 = normC(Fapprox_r,Fapprox_i,M*M*M);
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


