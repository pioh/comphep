#include<math.h>
double sigmoid(double x)
{
  return 1./(1.+exp(-x));
}
void rnnfun1(double *rin,double *rout)
{
  rout[0] = 1.;
  return;
}

float rnnfun(double *rin)
{
  return 1.;
}
