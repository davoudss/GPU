#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const double pi = 3.141592653589793;

float fgg(float x,int c, float tau, int Mr);
double fggD(double x,int c, double tau, int Mr);

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
  
  g = exp(-(2.f*pi*c/Mr-x-2.f*pi)*(2.f*pi*c/Mr-x-2.f*pi)/(4.f*tau)) 
    +exp(-(2.f*pi*c/Mr-x)*(2.f*pi*c/Mr-x)/(4.f*tau))
    +exp(-(2.f*pi*c/Mr-x+2.f*pi)*(2.f*pi*c/Mr-x+2.f*pi)/(4.f*tau));

  return g;
}

