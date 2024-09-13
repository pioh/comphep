/* Author: Slava Bunichev 2012 */

#include "width_func.h"
#include "4_vector.h"
#include "chep_limits.h"
#include "files.h"
#include "unix_utils.h"
#include <math.h>
#include <stdio.h>

#define mod(z) (((z)<0)? -(z):(z)) 


double myfunc1 (double Mass, double MR, double keyp)
{
  double result, prod12, ym, Fym, as, sqr, logs;

  if(MR!=0) ym = 4.*(Mass*Mass)/(MR*MR);
  else
    {
     prod12 = pvect[0]*pvect[4] -pvect[1]*pvect[4+1] -pvect[2]*pvect[4+2] -pvect[3]*pvect[4+3];

     if(prod12==0) prod12=(125.*125.)*0.5;

     ym = 2.*(Mass*Mass)/prod12;
    }
  
  as = asin(1./sqrt(fabs(ym)));
  sqr = sqrt(fabs(1.-ym));
  logs = log((1.+sqr)/(1.-sqr));
  
  if(ym >= 1.0) Fym = as*as; 
  else          Fym = -0.25*(logs*logs-9.869587728);  

  if(keyp < 1.5)   result = ym*(1.+(1.-ym)*Fym);
  else             result = -(2. + 3.*ym + 3.*ym*(2.-ym)*Fym);
  
  return result;
}



double myfunc2 (double Mass)
{
  double result, prod12, ym;  

  prod12 = pvect[0]*pvect[4] -pvect[1]*pvect[4+1] -pvect[2]*pvect[4+2] -pvect[3]*pvect[4+3];

  if(prod12==0) prod12=(125.*125.)*0.5;

  ym = 2.*(Mass*Mass)/prod12; 

  return ym;
}

double myfunc3 (double ym, double keyp)
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

double myfunc4 (double ym, double keyp)
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



static double roundto(double value)
{
 char numstring[20] = {0}; 
 sprintf(numstring,"%.6e",value);
 sscanf(numstring,"%le",&value);
 return value;
}


double width1 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);

/*  inpf=fopen(scat("%swidth1.txt", pathtocomphep),"r"); */

  sprintf(pathstr,"%swidth1.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file width1.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of width1 parameters\n");
  fclose(inpf);
  return 0;             
}

double width2 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);


  sprintf(pathstr,"%swidth2.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file width2.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of width2 parameters\n");
  fclose(inpf);
  return 0;             
}

double width3 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);


  sprintf(pathstr,"%swidth3.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file width3.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of width3 parameters\n");
  fclose(inpf);
  return 0;             
}


double width4 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);


  sprintf(pathstr,"%swidth4.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file width4.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of width4 parameters\n");
  fclose(inpf);
  return 0;             
}


double width5 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);


  sprintf(pathstr,"%swidth5.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file width5.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of width5 parameters\n");
  fclose(inpf);
  return 0;             
}

double koeff1 (double a, double c)
{ 
  FILE * inpf;
  char pathstr[500] = {0}; 
  double x, y, W, delta=1.0e-14;

  a = roundto(a);
  c = roundto(c);


  sprintf(pathstr,"%skoeff1.txt", pathtocomphep);

  inpf=fopen(pathstr,"r");
 
  if( inpf == NULL )
    { fprintf(stderr,"file koeff1.txt not found\n");
      exit(2);
    } 

  while(fscanf(inpf,"%le %le %le ",&x,&y,&W)!=EOF) 
    if( (mod(a-x)<delta) && (mod(c-y)<delta) ) 
      {
       fclose(inpf);
       return W;
      }

  fprintf(stderr,"\nwrong range of koeff1 parameters\n");
  fclose(inpf);
  return 0;             
}

