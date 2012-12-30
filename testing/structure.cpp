#include <iostream>

typedef struct
{
  double *real;
  double *imag;
}COMP;



int main()
{
  COMP z;
  double r[3]={1,2,3};
  double i[3]={5,6,7};

  z.real=r;
  z.imag=i;
  

  std::cout << z.real[1] << z.imag[1] << std::endl;

  return 0;

}
