void ccopy ( int n, float x[], float y[] );
void cfft2 ( int n, float x[], float y[], float w[], float sgn );
void cffti ( int n, float w[] );
void fft (float* x, float* y, int n, float sgn);
void step ( int n, int mj, float a[], float b[], float c[], float d[], float w[], float sgn );
void fftshift(float* x,int N);

void ccopyD ( int n, double x[], double y[] );
void cfft2D ( int n, double x[], double y[], double w[], double sgn );
void cfftiD ( int n, double w[] );
void fftD (double* x, double* y, int n, double sgn);
void stepD ( int n, int mj, double a[], double b[], double c[], double d[], double w[], double sgn );
void fftshiftD(double* x,int N);

/******************************************************************************/

void fftshift(float* x,int N)
{
  float temp;
  for(int k=0;k<N/4;k++){
      temp = x[2*k];
      x[2*k] = x[2*k+N/2];
      x[2*k+N/2] = temp;

      temp = x[2*k+1];
      x[2*k+1] = x[2*k+N/2+1];
      x[2*k+N/2+1] = temp;
  }
}



/******************************************************************************/
void fftshiftD(double* x,int N)
{
  double temp;
  for(int k=0;k<N/4;k++){
      temp = x[2*k];
      x[2*k] = x[2*k+N/2];
      x[2*k+N/2] = temp;

      temp = x[2*k+1];
      x[2*k+1] = x[2*k+N/2+1];
      x[2*k+N/2+1] = temp;
  }
}

/******************************************************************************/

void fft (float* x,float* y, int n, float sgn)

{
  float *w;
  float *z;

  w = ( float * ) malloc (     n * sizeof ( float ) );
  z = ( float * ) malloc ( 2 * n * sizeof ( float ) );  
  
  cffti ( n, w );
  cfft2 ( n, x, y, w, sgn );

  free ( w );
  free ( z );
}


/******************************************************************************/

void fftD (double* x,double* y, int n, double sgn)

{
  double *w;
  double *z;

  w = ( double * ) malloc (     n * sizeof ( double ) );
  z = ( double * ) malloc ( 2 * n * sizeof ( double ) );  
  
  cfftiD ( n, w );
  cfft2D ( n, x, y, w, sgn );

  free ( w );
  free ( z );
}


/******************************************************************************/

void ccopy ( int n, float x[], float y[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    y[i*2+0] = x[i*2+0];
    y[i*2+1] = x[i*2+1];
   }
  return;
}


/******************************************************************************/

void ccopyD ( int n, double x[], double y[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    y[i*2+0] = x[i*2+0];
    y[i*2+1] = x[i*2+1];
   }
  return;
}


/******************************************************************************/
void cfft2 ( int n, float x[], float y[], float w[], float sgn )
{
  int j;
  int m;
  int mj;
  int tgle;

   m = ( int ) ( log ( ( float ) n ) / log ( 1.99 ) );
   mj   = 1;
/*
  Toggling switch for work array.
*/
  tgle = 1;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  if ( n == 2 )
  {
    return;
  }

  for ( j = 0; j < m - 2; j++ )
  {
    mj = mj * 2;
    if ( tgle )
    {
      step ( n, mj, &y[0*2+0], &y[(n/2)*2+0], &x[0*2+0], &x[mj*2+0], w, sgn );
      tgle = 0;
    }
    else
    {
      step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );
      tgle = 1;
    }
  }
/* 
  Last pass through data: move Y to X if needed.
*/
  if ( tgle ) 
  {
    ccopy ( n, y, x );
  }

  mj = n / 2;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  return;
}


/******************************************************************************/
void cfft2D ( int n, double x[], double y[], double w[], double sgn )
{
  int j;
  int m;
  int mj;
  int tgle;

   m = ( int ) ( log ( ( double ) n ) / log ( 1.99 ) );
   mj   = 1;
/*
  Toggling switch for work array.
*/
  tgle = 1;
  stepD ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  if ( n == 2 )
  {
    return;
  }

  for ( j = 0; j < m - 2; j++ )
  {
    mj = mj * 2;
    if ( tgle )
    {
      stepD ( n, mj, &y[0*2+0], &y[(n/2)*2+0], &x[0*2+0], &x[mj*2+0], w, sgn );
      tgle = 0;
    }
    else
    {
      stepD ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );
      tgle = 1;
    }
  }
/* 
  Last pass through data: move Y to X if needed.
*/
  if ( tgle ) 
  {
    ccopyD ( n, y, x );
  }

  mj = n / 2;
  stepD ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  return;
}

/******************************************************************************/
void cffti ( int n, float w[] )
{
  float arg;
  float aw;
  int i;
  int n2;
  const float pi = 3.141592653589793;

  n2 = n / 2;
  aw = 2.0 * pi / ( ( float ) n );

  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( float ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}

/******************************************************************************/
void cfftiD ( int n, double w[] )
{
  double arg;
  double aw;
  int i;
  int n2;
  const double pi = 3.141592653589793;

  n2 = n / 2;
  aw = 2.0 * pi / ( ( double ) n );

  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( double ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}


/******************************************************************************/
void step ( int n, int mj, float a[], float b[], float c[],
  float d[], float w[], float sgn )
{
  float ambr;
  float ambu;
  int j;
  int ja;
  int jb;
  int jc;
  int jd;
  int jw;
  int k;
  int lj;
  int mj2;
  float wjw[2];

  mj2 = 2 * mj;
  lj  = n / mj2;

  //# pragma omp parallel				  \
  //private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw ) \
  //shared ( a, b, c, d, lj, mj, mj2, sgn, w )

  //# pragma omp for

  for ( j = 0; j < lj; j++ )
  {
    jw = j * mj;
    ja  = jw;
    jb  = ja;
    jc  = j * mj2;
    jd  = jc;

    wjw[0] = w[jw*2+0]; 
    wjw[1] = w[jw*2+1];

    if ( sgn > 0.0 ) 
    {
      wjw[1] = - wjw[1];
    }

    for ( k = 0; k < mj; k++ )
    {
      c[(jc+k)*2+0] = a[(ja+k)*2+0] + b[(jb+k)*2+0];
      c[(jc+k)*2+1] = a[(ja+k)*2+1] + b[(jb+k)*2+1];

      ambr = a[(ja+k)*2+0] - b[(jb+k)*2+0];
      ambu = a[(ja+k)*2+1] - b[(jb+k)*2+1];

      d[(jd+k)*2+0] = wjw[0] * ambr - wjw[1] * ambu;
      d[(jd+k)*2+1] = wjw[1] * ambr + wjw[0] * ambu;
    }
  }
  return;
}


/******************************************************************************/
void stepD ( int n, int mj, double a[], double b[], double c[],
  double d[], double w[], double sgn )
{
  double ambr;
  double ambu;
  int j;
  int ja;
  int jb;
  int jc;
  int jd;
  int jw;
  int k;
  int lj;
  int mj2;
  double wjw[2];

  mj2 = 2 * mj;
  lj  = n / mj2;

  //# pragma omp parallel				  \
  //private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw ) \
  //shared ( a, b, c, d, lj, mj, mj2, sgn, w )

  //# pragma omp for

  for ( j = 0; j < lj; j++ )
  {
    jw = j * mj;
    ja  = jw;
    jb  = ja;
    jc  = j * mj2;
    jd  = jc;

    wjw[0] = w[jw*2+0]; 
    wjw[1] = w[jw*2+1];

    if ( sgn > 0.0 ) 
    {
      wjw[1] = - wjw[1];
    }

    for ( k = 0; k < mj; k++ )
    {
      c[(jc+k)*2+0] = a[(ja+k)*2+0] + b[(jb+k)*2+0];
      c[(jc+k)*2+1] = a[(ja+k)*2+1] + b[(jb+k)*2+1];

      ambr = a[(ja+k)*2+0] - b[(jb+k)*2+0];
      ambu = a[(ja+k)*2+1] - b[(jb+k)*2+1];

      d[(jd+k)*2+0] = wjw[0] * ambr - wjw[1] * ambu;
      d[(jd+k)*2+1] = wjw[1] * ambr + wjw[0] * ambu;
    }
  }
  return;
}



