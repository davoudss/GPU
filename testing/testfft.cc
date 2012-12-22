# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <float.h>
# include <omp.h>

int main ( void );
void ccopy ( int n, float x[], float y[] );
void cfft2 ( int n, float x[], float y[], float w[], float sgn );
void cffti ( int n, float w[] );
double cpu_time ( void );
float ggl ( float *ds );
void step ( int n, int mj, float a[], float b[], float c[], float d[], 
  float w[], float sgn );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/* 
  Purpose:

    MAIN is the main program for FFT.

  Discussion:

    The "complex" vector A is actually stored as a float vector B.

    The "complex" vector entry A[I] is stored as:

      B[I*2+0], the real part,
      B[I*2+1], the imaginary part.

  Modified:

    16 July 2008

  Author:

    C original version by Wesley Petersen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.
*/
{
  double ctime;
  double ctime1;
  double ctime2;
  float error;
  int first;
  int flops;
  float fnm1;
  int i;
  int icase;
  int it;
  int ln2;
  double mflops;
  int n;
  int nits = 10000;
  static float seed;
  float sgn;
  float *w;
  float *x;
  float *y;
  float *z;
  float z0;
  float z1;


  seed  = 331.0;
  n = 4;
/*
  LN2 is the log base 2 of N.  Each increase of LN2 doubles N.
*/
  for ( ln2 = 1; ln2 <= 1; ln2++ )
  {
    //  n = 2 * n;
/*
  Allocate storage for the complex arrays W, X, Y, Z.  

  We handle the complex arithmetic,
  and store a complex number as a pair of floats, a complex vector as a doubly
  dimensioned array whose second dimension is 2. 
*/
    w = ( float * ) malloc (     n * sizeof ( float ) );
    x = ( float * ) malloc ( 2 * n * sizeof ( float ) );
    y = ( float * ) malloc ( 2 * n * sizeof ( float ) );
    z = ( float * ) malloc ( 2 * n * sizeof ( float ) );

    first = 1;

    for ( icase = 0; icase < 1; icase++ )
    {
      if ( first )
      {
        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          z0 = ggl ( &seed );
          z1 = ggl ( &seed );
          x[i] = z0;
          z[i] = z0;
          x[i+1] = z1;
          z[i+1] = z1;
        }
      } 
      else
      {
        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          z0 = 0.0;              /* real part of array */
          z1 = 0.0;              /* imaginary part of array */
          x[i] = z0;
          z[i] = z0;           /* copy of initial real data */
          x[i+1] = z1;
          z[i+1] = z1;         /* copy of initial imag. data */
        }
      }
/* 
  Initialize the sine and cosine tables.
*/
      cffti ( n, w );
      printf("x: \n");
      for(int k=0;k<2*n;k++)
	printf("%f\t",x[k]);

      printf("\ny: \n");
      for(int k=0;k<2*n;k++)
	printf("%f\t",y[k]);



/* 
  Transform forward, back 
*/
      if ( first )
      {
        sgn = + 1.0;
        cfft2 ( n, x, y, w, sgn );
	printf("\nx: \n");
	for(int k=0;k<2*n;k++)
	  printf("%f\t",x[k]);

	printf("\ny: \n");
	for(int k=0;k<2*n;k++)
	  printf("%f\t",y[k]);

	printf("\n");
/* 
  Results should be same as the initial data multiplied by N.
*/
        fnm1 = 1.0 / ( float ) n;
        error = 0.0;
        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          error = error 
          + pow ( z[i]   - fnm1 * x[i], 2 )
          + pow ( z[i+1] - fnm1 * x[i+1], 2 );
        }
        error = sqrt ( fnm1 * error );
        first = 0;
      }
      else
      {
        ctime1 = omp_get_wtime ( );
        for ( it = 0; it < nits; it++ )
        {
          sgn = + 1.0;
          cfft2 ( n, x, y, w, sgn );
        }
        ctime2 = omp_get_wtime ( );
        ctime = ctime2 - ctime1;

      }
    }
    if ( ( ln2 % 4 ) == 0 ) 
    {
      nits = nits / 10;
    }
    if ( nits < 1 ) 
    {
      nits = 1;
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  }
  printf ( "\n" );
  printf ( "FFT:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void ccopy ( int n, float x[], float y[] )

/******************************************************************************/
/*
  Purpose:

    CCOPY copies a complex vector.

  Discussion:

    The "complex" vector A[N] is actually stored as a float vector B[2*N].

    The "complex" vector entry A[I] is stored as:

      B[I*2+0], the real part,
      B[I*2+1], the imaginary part.

  Modified:

    26 April 2008

  Author:

    C original version by Wesley Petersen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parameters:

    Input, int N, the length of the vector.

    Input, "complex" X[N], the vector to be copied.

    Output, "complex" Y[N], a copy of X.
*/
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

/******************************************************************************/
/*
  Purpose:

    CFFT2 performs a complex Fast Fourier Transform.

  Modified:

    26 April 2008

  Author:

    C original version by Wesley Petersen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parameters:

    Input, int N, the size of the array to be transformed.

    Input/output, float X[2*N], the data to be transformed.  
    On output, the contents of X have been overwritten by work information.

    Output, float Y[2*N], the forward or backward FFT of X.

    Input, float W[N], a table of sines and cosines.

    Input, float SGN, is +1 for a "forward" FFT and -1 for a "backward" FFT.
*/
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

void cffti ( int n, float w[] )

/******************************************************************************/
/*
  Purpose:

    CFFTI sets up sine and cosine tables needed for the FFT calculation.

  Modified:

    26 April 2008

  Author:

    C original version by Wesley Petersen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parameters:

    Input, int N, the size of the array to be transformed.

    Output, float W[N], a table of sines and cosines.
*/
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
/*******************************************************************************/

double cpu_time ( void )

/*******************************************************************************/
/*
  Purpose:
 
    CPU_TIME reports the total CPU time for a program.

  Modified:

    27 September 2005

  Author:

    John Burkardt

  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

float ggl ( float *seed )

/******************************************************************************/
/* 
  Purpose:

    GGL generates uniformly distributed pseudorandom real numbers in [0,1]. 

  Modified:

    16 July 2008

  Author:

    C original version by Wesley Petersen, M Troyer, I Vattulainen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parameters:

    Input/output, float *SEED, used as a seed for the sequence.

    Output, float GGL, the next pseudorandom value.
*/
{
  double d2 = 0.2147483647e10;
  double t;
  float value;

  t = ( double ) *seed;
  t = fmod ( 16807.0 * t, d2 );
  *seed = ( float ) t;
  value = ( float ) ( ( t - 1.0 ) / ( d2 - 1.0 ) );

  return value;
}
/******************************************************************************/

void step ( int n, int mj, float a[], float b[], float c[],
  float d[], float w[], float sgn )

/******************************************************************************/
/*
  Purpose:

    STEP carries out one step of the workspace version of CFFT2.

  Modified:

    16 July 2008

  Author:

    C original version by Wesley Petersen
    Modifications by John Burkardt

  Reference:

    Wesley Petersen, Peter Arbenz, 
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.
*/
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

# pragma omp parallel \
    private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw ) \
    shared ( a, b, c, d, lj, mj, mj2, sgn, w )

# pragma omp for

  for ( j = 0; j < lj; j++ )
  {
    jw = j * mj;
    ja  = jw;
    jb  = ja;
    jc  = j * mj2;
    jd  = jc;

    wjw[0] = w[jw*2+0]; 
    wjw[1] = w[jw*2+1];

    if ( sgn < 0.0 ) 
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
