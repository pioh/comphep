/*
* Copyright (C) 2002-2017, CompHEP Collaboration
* ------------------------------------------------------
  autor Viacheslav Bunichev
*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "service2/include/4_vector.h"
#include "userFactor.h"



static double ym_func(double Mass)
{
  double result, prod12, ym;  

  prod12 = pvect[0]*pvect[4] -pvect[1]*pvect[4+1] -pvect[2]*pvect[4+2] -pvect[3]*pvect[4+3];

  if(prod12==0) prod12=(125.*125.)*0.5;

  ym = 2.*(Mass*Mass)/prod12; 

  return ym;
}

static double real_func(double ym, double keyp)
{
  double result, Fym, as, sqr, logs;  

  as = asin(1./sqrt(fabs(ym)));
  sqr = sqrt(fabs(1.-ym));
  logs = log((1.+sqr)/(1.-sqr));

  
  if(ym >= 1.0) Fym = as*as; 
  else          Fym = -0.25*(logs*logs-9.869587728);  

  if(keyp < 1.5)   result = ym*(1.+(1.-ym)*Fym);
  else             result = -(2. + 3.*ym + 3.*ym*(2.-ym)*Fym);
  
  return result;
}

static double im_func(double ym, double keyp)
{
  double result, Fym, as, sqr, logs;  

  as = asin(1./sqrt(fabs(ym)));
  sqr = sqrt(fabs(1.-ym));
  logs = log((1.+sqr)/(1.-sqr));
  
  if(ym >= 1.0) Fym = 0; 
  else          Fym = 0.5*(logs*3.14159265358);  

  if(keyp < 1.5)  result = ym*(1.-ym)*Fym;
  else            result = -3.*ym*(2.-ym)*Fym;
  
  return result;
}


double userfactorFun(int factorkey)
{  
  double yt=0, loopt=0, Imlt=0, Mtop=173.4, factorFun = 1.;

  if(factorkey==0) return 1.;

  else
   { 
    yt = ym_func(Mtop);
    loopt = real_func(yt,1);
    Imlt = im_func(yt,1);

    factorFun = loopt*loopt + Imlt*Imlt;
     
/*    printf("%e \n",factorFun); */

    return factorFun;
   }

}

