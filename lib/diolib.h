#include<fftw3.h>

int printing(fftw_complex* x, int N, int dim);
int printing(float*x, int N, int dim);
int printing(double*x, double *y, int N,int dim);
int printing(double* x, int N, int dim, char ch);
int printing(double*x, int N, int dim);

int printFile(float *x, int N,int n);
int printFileD(double *x, int N,int n);
void DisplayMatrix(char*,double*, double*, int,int);


/***********************************************************************/
int printing(double*x, double *y, int N,int dim)
{
  printf("\n\n");
  if(dim==1)
    for (int s1=0;s1<N;s1++)
      printf("(%.4e,%.4e)\n",x[s1],y[s1]);
  else if(dim==2)
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	printf("(%.4e,%.4e)\n",x[s1*N+s2],y[s1*N+s2]);
      }
      printf("\n");
    }
  else
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	for (int s3=0;s3<N;s3++){
	  printf("(%.4e,%.4e)\n",x[s1*N*N+s2*N+s3],y[s1*N*N+s2*N+s3]);
	}
	printf("\n");
      }
      printf("\n\n");
    }
  
  return 0;
}

/***********************************************************************/
int printing(fftw_complex* x, int N, int dim)
{
  printf("\n\n");
  if(dim==1)
    for (int s1=0;s1<N;s1++)
      printf("(%.4e,%.4e)\n",x[s1][0],x[s1][0]);
  else if(dim==2)
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	printf("(%.4e,%.4e)\n",x[s1*N+s2][0],x[s1*N+s2][1]);
      }
      printf("\n");
    }
  else
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	for (int s3=0;s3<N;s3++){
	  printf("(%.4e,%.4e)\n",x[s1*N*N+s2*N+s3][0],x[s1*N*N+s2*N+s3][1]);
	}
	printf("\n");
      }
      printf("\n\n");
    }
  
  return 0;
}


/***********************************************************************/
int printing(double* x, int N, int dim, char ch)
{
  if(ch=='c' || ch=='C'){

    printf("\n\n");
    if(dim==1)
      for (int s1=0;s1<N;s1++)
	printf("(%.4e,%.4e)\n",x[2*s1],x[2*s1+1]);
    else if(dim==2)
      for (int s1=0;s1<N;s1++){
	for (int s2=0;s2<N;s2++){
	  printf("(%.4e,%.4e)\n",x[2*(s1*N+s2)],x[2*(s1*N+s2)+1]);
	}
	printf("\n");
      }
    else
      for (int s1=0;s1<N;s1++){
	for (int s2=0;s2<N;s2++){
	  for (int s3=0;s3<N;s3++){
	    printf("(%.4e,%.4e)\n",x[2*(s1*N*N+s2*N+s3)],x[2*(s1*N*N+s2*N+s3)+1]);
	  }
	  printf("\n");
	}
	printf("\n\n");
      }
  }
  return 0;
}



/***********************************************************************/
int printing(float *x, int N, int dim)
{
  printf("\n\n");
  if(dim==1)
    for (int s1=0;s1<N;s1++)
      printf("%.4e\n",x[s1]);
  else if(dim==2)
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	printf("%.4e\n",x[s1*N+s2]);
      }
      printf("\n");
    }
  else
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	for (int s3=0;s3<N;s3++){
	  printf("%.4e\n",x[s1*N*N+s2*N+s3]);
	}
	printf("\n");
      }
      printf("\n\n");
    }
  
  return 0;

}

/***********************************************************************/
int printing(double *x, int start, int N, int dim)
//N is the number of elements, start is the starting double number.
{
  printf("\n\n");
  if(dim==1)
    for (int s1=0;s1<N;s1++)
      printf("%.4e\n",x[s1]);
  else if(dim==2)
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	printf("%.4e\n",x[s1*N+s2]);
      }
      printf("\n");
    }
  else
    for (int s1=0;s1<N;s1++){
      for (int s2=0;s2<N;s2++){
	for (int s3=0;s3<N;s3++){
	  printf("%.4e\n",x[s1*N*N+s2*N+s3]);
	}
	printf("\n");
      }
      printf("\n\n");
    }
  
  return 0;


}

/***********************************************************************/
int printFile(float *x, int N,int n)
//N is the number of complex numbers
{
  FILE *fp;
  if(n==1)
    fp = fopen("matrix.m","w");
  else
    fp = fopen("matrix.m","a+");

  fprintf(fp,"F%d= [\n",n);

  for (int s=0;s<2*N;s++)
    fprintf(fp,"%+f%+fi\n",x[2*s],x[2*s+1]);

  fprintf(fp,"];\n");
  if(n==2)
    fprintf(fp,"plot(real(Fexact),imag(Fexact),'.',real(Fapprox),imag(Fapprox),'ro')\n legend('Exact','Appproximation')\n disp(['Relative error:', num2str(abs(norm(Fexact - Fapprox))/norm(Fexact))])");

  fclose(fp);
  return 0;

}

/***********************************************************************/
int printFileD(double *x, int N,int n)
//N is the number of complex numbers
{
  FILE *fp;
  if(n==1)
    fp = fopen("matrix.m","w");
  else
    fp = fopen("matrix.m","a+");

  fprintf(fp,"F%d= [\n",n);

  for (int s=0;s<2*N;s++)
    fprintf(fp,"%+f%+fi\n",x[2*s],x[2*s+1]);

  fprintf(fp,"];\n");
  if(n==2)
    fprintf(fp,"plot(real(Fexact),imag(Fexact),'.',real(Fapprox),imag(Fapprox),'ro')\n legend('Exact','Appproximation')\n disp(['Relative error:', num2str(abs(norm(Fexact - Fapprox))/norm(Fexact))])");

  fclose(fp);
  return 0;

}


void DisplayMatrix(char           *Name, 
		   double       *Data_r,
		   double       *Data_i, 
		   int                M, 
		   int                N
		   )
{
  int m, n;
  printf("%s = \n", Name);
  for(m = 0; m < M; m++, printf("\n"))
    for(n = 0; n < N; n++)
      printf("%+8.4f%+8.4fi ", Data_r[m + M*n],Data_i[m+M*n]);
}
