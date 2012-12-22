#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.141592653589793;

float fgg(float x,int c, float tau, int Mr);
double fggD(double x,int c, double tau, int Mr);
double fgg3(double x,double y, double z,int c1,int c2, int c3, double tau, int Mr);

float fgg(float x,int c, float tau, int Mr)
{
  float g=0.,l;
  
  g = exp(-(2.f*pi*c/Mr-x-2.f*pi)*(2.f*pi*c/Mr-x-2.f*pi)/(4.f*tau)) 
    +exp(-(2.f*pi*c/Mr-x)*(2.f*pi*c/Mr-x)/(4.f*tau))
    +exp(-(2.f*pi*c/Mr-x+2.f*pi)*(2.f*pi*c/Mr-x+2.f*pi)/(4.f*tau));

  return g;
}



double fggD(double x,int c, double tau, int Mr)
{
  double g=0.,l;
  
  g = exp(-(2.*pi*c/Mr-x-2.*pi)*(2.*pi*c/Mr-x-2.*pi)/(4.*tau)) 
    +exp(-(2.*pi*c/Mr-x)*(2.*pi*c/Mr-x)/(4.*tau))
    +exp(-(2.*pi*c/Mr-x+2.*pi)*(2.*pi*c/Mr-x+2.*pi)/(4.*tau));

  return g;
}


double fgg3(double x, double y, double z, int c1, int c2,int c3, double tau, int Mr)
{
  double g=0.;
  
  for (int l1=-1;l1<2;l1++)
    for (int l2=-1;l2<2;l2++)
      for (int l3=-1;l3<2;l3++)	
        g += exp(-((2.*pi*c1/Mr-x-2.*l1*pi)*(2.*pi*c1/Mr-x-2.*l1*pi)+(2.*pi*c2/Mr-y-2.*l2*pi)*(2.*pi*c2/Mr-y-2.*l2*pi)+(2.*pi*c3/Mr-z-2.*l3*pi)*(2.*pi*c3/Mr-z-2.*l3*pi))/(4.*tau));


  return g;
}

