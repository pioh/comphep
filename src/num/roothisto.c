/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Slava Bunichev
* ------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/files.h"
#include "service2/include/unix_utils.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "chep_crt/include/chep_crt.h"

#include "roothisto.h"

static void
readtext (char *fname)
{
  FILE *txt;

  trim (fname);
  txt = fopen (fname, "r");
  if (txt == NULL)
    {
      warnanykey (10, 10, " File not found ");
      return;
    }
  showtext (28, 10, 70, 15, "", txt);
  fclose (txt);
}


static int
graphcolor (int first)
{
  int color;

  if (first == 1)
    color = 4;
  if (first == 2)
    color = 2;
  if (first == 3)
    color = 3;
  if (first == 4)
    color = 6;
  if (first == 5)
    color = 7;
  if (first == 6)
    color = 9;
  if (first == 7)
    color = 6;
  if (first == 8)
    color = 8;
  if (first == 9)
    color = 3;
  if (first == 10)
    color = 1;
  if (first == 11)
    color = 2;
  if (first == 12)
    color = 4;
  if (first == 13)
    color = 5;
  if (first == 14)
    color = 8;
  if (first == 15)
    color = 9;

  if (first > 15)
    color = 1;

  return (color);
}


static int
graphstyle (int first)
{
  int style;

  if (first == 1)
    style = 1;
  if (first == 2)
    style = 2;
  if (first == 3)
    style = 3;
  if (first == 4)
    style = 4;

  if (first == 5)
    style = 1;
  if (first == 6)
    style = 2;
  if (first == 7)
    style = 3;
  if (first == 8)
    style = 4;

  if (first == 9)
    style = 1;
  if (first == 10)
    style = 2;
  if (first == 11)
    style = 3;
  if (first == 12)
    style = 4;

  if (first == 13)
    style = 1;
  if (first == 14)
    style = 2;
  if (first == 15)
    style = 3;

  if (first > 15)
    style = 1;

  return (style);
}


static void
graph_identity (FILE * outfile, int first)
{
  if (first == 1)
    f_printf (outfile, "//solid, light blue\n");
  if (first == 2)
    f_printf (outfile, "//dashed, light red\n");
  if (first == 3)
    f_printf (outfile, "//doted, dark green\n");
  if (first == 4)
    f_printf (outfile, "//dots and dashes, dark magenta\n");

  if (first == 5)
    f_printf (outfile, "//solid, dark green and blue\n");
  if (first == 6)
    f_printf (outfile, "//dashed, medium blue\n");
  if (first == 7)
    f_printf (outfile, "//doted, light magenta\n");
  if (first == 8)
    f_printf (outfile, "//dots and dashes, medium green\n");

  if (first == 9)
    f_printf (outfile, "//solid, light green\n");
  if (first == 10)
    f_printf (outfile, "//dashed, grey\n");
  if (first == 11)
    f_printf (outfile, "//doted, dark red\n");
  if (first == 12)
    f_printf (outfile, "//dots and dashes, dark blue\n");

  if (first == 13)
    f_printf (outfile, "//solid, dark yellow\n");
  if (first == 14)
    f_printf (outfile, "//dashed, very dark green\n");
  if (first == 15)
    f_printf (outfile, "//doted, very dark blue\n");

  if (first > 15)
    f_printf (outfile, "//solid, black\n");
}

static void
writeroothist (char *menustr, int numstring)
{
  int i, k;
  int combine = 0;
  int dim;
  int old_dim;
  int first = 1;
  int color;
  int style;
  char c;
  char old_xstr[50] = { 0 };
  char upstr[50] = { 0 };
  char x_str[50] = { 0 };
  char y_str[50] = { 0 };
  double old_xMin;
  double old_xMax;
  double xMin, xMax, yMin, yMax;
  double dX;
  double totYMin = 0, totYMax = 0;

  char xaxistitle[128];
  char yaxistitle[128];
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf;

  nextFileName (outname, "combine_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
	  if (infile == NULL)
	    {
	      warnanykey (10, 10, " File not found ");
	      return;
	    }

	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin, &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);

	  if (first != 1)
	    {
	      if (combine != 0)
		{
		  warnanykey (10, 10, " Combined distribution ");
		  fclose (infile);
		  return;
		}
	      if (strncmp (x_str, old_xstr, 2) != 0)
		{
		  warnanykey (10, 10, " Different functions ");
		  fclose (infile);
		  return;
		}

	      if (xMin != old_xMin || xMax != old_xMax)
		{
		  warnanykey (10, 10, " Different X-range ");
		  fclose (infile);
		  return;
		}
	      if (dim != old_dim)
		{
		  warnanykey (10, 10, " Different binning ");
		  fclose (infile);
		  return;
		}
	    }
	  fclose (infile);

	  if (first == 1)
	    {
	      totYMin = yMin;
	      totYMax = yMax;
	    }

	  if (yMin < totYMin)
	    totYMin = yMin;
	  if (yMax > totYMax)
	    totYMax = yMax;

	  strcpy (old_xstr, x_str);

	  for (k = 0; k < 50; k++)
	    {
	      upstr[k] = 0;
	      x_str[k] = 0;
	      y_str[k] = 0;
	    }

	  old_xMin = xMin;
	  old_xMax = xMax;
	  old_dim = dim;
	  first++;
	}
    }

  rf = fopen (outname_C, "w");
  if (NULL == rf)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin,
		  &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);
	  fclose (infile);

	  if (first == 1)
	    {
	      f_printf (rf, "%s\n", upstr);
	      f_printf (rf, "%s\n", x_str);
	      f_printf (rf, "//X from %e to %e N_bins= %d\n", xMin, xMax, dim);
	      f_printf (rf, "%s\n", y_str);
	      f_printf (rf, "//Y from %e to %e\n", yMin, yMax);
	      f_printf (rf, "//combine= 1\n");
	    }
	  else
	    {
	      f_printf (rf, "//GRAPH%d\n", first);
	      f_printf (rf, "%s\n", upstr);
	      f_printf (rf, "%s\n", x_str);
	    }

	  graph_identity (rf, first);

	  for (k = 0; k < 50; k++)
	    {
	      upstr[k] = 0;
	      x_str[k] = 0;
	      y_str[k] = 0;
	    }
	  first++;
	}
    }

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin,
		  &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);

	  if (first == 1)
	    {
	      f_printf (rf, "\n{\n");
	      f_printf (rf, "gROOT->Reset();\n");
	      f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
	      f_printf (rf, "c1->SetFillColor(0);\n");
	      f_printf (rf, "c1->SetBorderMode(0);\n");
	      f_printf (rf, "c1->SetFrameBorderMode(0);\n");
	      f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
	      f_printf (rf, "//c1->SetLogy(1);  //Set Log. Scale;\n\n\n");
	      if (dim != 0) dX = (xMax - xMin) / dim;
	      f_printf (rf, "double Xarray[%d];\n", dim);
	      f_printf (rf, "for(int i=0;i<%d;i++) Xarray[i]=%e+i*%e;\n\n", dim, xMin, dX);
	    }
	  f_printf (rf, "// GRAPH %d\n", first);

	  while ((getc (infile) != 'Y') ||
		 (getc (infile) != 'a') ||
		 (getc (infile) != 'r') ||
		 (getc (infile) != 'r') ||
		 (getc (infile) != 'a') ||
		 (getc (infile) != 'y'))
	    {
	    }
	  while ((c = getc (infile)) != '{')
	    {
	    }

	  f_printf (rf, "double Yarray%d[%d]={", first, dim);
	  while ((c = getc (infile)) != '}')
	      putc (c, rf);
	  fclose (infile);

	  f_printf (rf, "};\n");
	  f_printf (rf, "TGraph *gr%d = new ", first);
	  f_printf (rf, "TGraph(%d,Xarray,Yarray%d);\n", dim, first);
	  if (first == 1)
	    f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
	  f_printf (rf, "gr%d->SetLineWidth(2);\n", first);
	  color = graphcolor (first);
	  f_printf (rf, "gr%d->SetLineColor(%d);\n", first, color);
	  style = graphstyle (first);
	  f_printf (rf, "gr%d->SetLineStyle(%d);\n", first, style);

	  if (first == 1)
	    {
	      sprintf (xaxistitle, "\"%s\"", x_str);
	      sprintf (yaxistitle, "\"#sigma, fb\"");
	      if (!strncmp (x_str, "//Transverse momentum Pt%d", 23))
		{
		  sprintf (xaxistitle, "\"P_{T}, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d P_{T}, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Mass", 6))
		{
		  sprintf (xaxistitle, "\"M, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d M, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Energy E", 10))
		{
		  sprintf (xaxistitle, "\"E, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d E, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Angle", 7))
		{
		  sprintf (xaxistitle, "\"#theta, deg\"");
		  sprintf (yaxistitle, "\"d#sigma/d #theta, fb/%.0f deg\"", dX);
		}

	      if (!strncmp (x_str, "//pseudo-rapidity", 17))
		{
		  sprintf (xaxistitle, "\"#eta\"");
		  sprintf (yaxistitle, "\"d#sigma/d #eta, fb\"");
		}

	      if (!strncmp (x_str, "//Rapidity", 10))
		{
		  sprintf (xaxistitle, "\"#y\"");
		  sprintf (yaxistitle, "\"d#sigma/d #y, fb\"");
		}

	      if (!strncmp (x_str, "//Cosine", 8))
		{
		  sprintf (xaxistitle, "\"cos(#theta)\"");
		  sprintf (yaxistitle, "\"d#sigma/d cos(#theta), fb\"");
		}
	      f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
	      f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);

	      f_printf (rf, "gr1->SetMinimum(%e);\n", totYMin);
	      f_printf (rf, "gr1->SetMaximum(%e);\n", 1.1 * totYMax);
	      f_printf (rf, "gr1->Draw(\"AC\");\n\n");
	    }
	  else
	    {
	      f_printf (rf, "gr%d->Draw(\"C\");\n\n", first);
	    }
	  first++;
	}
    }

  f_printf (rf, "//c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "//c1->Print(\"%s.png\");\n}\n", outname);
  fclose (rf);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



static void
writeroothist_sum (char *menustr, int numstring)
{
  int i,j, k;
  int combine = 0;
  int dim;
  int old_dim;
  int first = 1;
  int color;
  int style;
  char c;
  char old_xstr[50] = { 0 };
  char upstr[50] = { 0 };
  char x_str[50] = { 0 };
  char y_str[50] = { 0 };
  double old_xMin;
  double old_xMax;
  double xMin, xMax, yMin, yMax;
  double dX;
  double totYMin = 0, totYMax = 0;


  char xaxistitle[128];
  char yaxistitle[128];
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf;

  nextFileName (outname, "combine_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')  
	{
	  FILE *infile = fopen (infname, "r");
	  if (infile == NULL)
	    {
	      warnanykey (10, 10, " File not found ");
	      return;
	    }

	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin, &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);

	  if (first != 1)
	    {
	      if (combine != 0)
		{
		  warnanykey (10, 10, " Combined distribution ");
		  fclose (infile);
		  return;
		}
	      if (strncmp (x_str, old_xstr, 2) != 0)
		{
		  warnanykey (10, 10, " Different functions ");
		  fclose (infile);
		  return;
		}

	      if (xMin != old_xMin || xMax != old_xMax)
		{
		  warnanykey (10, 10, " Different X-range ");
		  fclose (infile);
		  return;
		}
	      if (dim != old_dim)
		{
		  warnanykey (10, 10, " Different binning ");
		  fclose (infile);
		  return;
		}
	    }
	  fclose (infile);

	  if (first == 1)
	    {
	      totYMin = yMin;
	      totYMax = yMax;
	    }

	  if (yMin < totYMin)
	    totYMin = yMin;
	  if (yMax > totYMax)
	    totYMax = yMax;

	  strcpy (old_xstr, x_str);

	  for (k = 0; k < 50; k++)
	    {
	      upstr[k] = 0;
	      x_str[k] = 0;
	      y_str[k] = 0;
	    }

	  old_xMin = xMin;
	  old_xMax = xMax;
	  old_dim = dim;
	  first++;
	}
    }

  rf = fopen (outname_C, "w");
  if (NULL == rf)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')  
	{
	  FILE *infile = fopen (infname, "r");
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin,
		  &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);
	  fclose (infile);

	  if (first == 1)
	    {
	      f_printf (rf, "%s\n", upstr);
	      f_printf (rf, "%s\n", x_str);
	      f_printf (rf, "//X from %e to %e N_bins= %d\n", xMin, xMax, dim);
	      f_printf (rf, "%s\n", y_str);
	      f_printf (rf, "//Y from %e to %e\n", yMin, yMax);
	      f_printf (rf, "//combine= 1\n");
	    }
	  else
	    {
	      f_printf (rf, "//GRAPH%d\n", first);
	      f_printf (rf, "%s\n", upstr);
	      f_printf (rf, "%s\n", x_str);
	    }

	  graph_identity (rf, first);

	  for (k = 0; k < 50; k++)
	    {
	      upstr[k] = 0;
	      x_str[k] = 0;
	      y_str[k] = 0;
	    }
	  first++;
	}
    }

  first = 1;
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')  
	{
	  FILE *infile = fopen (infname, "r");
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      upstr[k] = c;
	      k++;
	    }
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      x_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//X from %le to %le N_bins= %d\n", &xMin,
		  &xMax, &dim);
	  k = 0;
	  while ((c = getc (infile)) != '\n')
	    {
	      y_str[k] = c;
	      k++;
	    }
	  fscanf (infile, "//Y from %le to %le", &yMin, &yMax);
	  fscanf (infile, "//combine= %d", &combine);

	  if (first == 1)
	    {
	      f_printf (rf, "\n{\n");
	      f_printf (rf, "gROOT->Reset();\n");
	      f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
	      f_printf (rf, "c1->SetFillColor(0);\n");
	      f_printf (rf, "c1->SetBorderMode(0);\n");
	      f_printf (rf, "c1->SetFrameBorderMode(0);\n");
	      f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
	      f_printf (rf, "//c1->SetLogy(1);  //Set Log. Scale;\n\n\n");
	      if (dim != 0) dX = (xMax - xMin) / dim;
	      f_printf (rf, "double Xarray[%d];\n", dim);
	      f_printf (rf, "for(int i=0;i<%d;i++) Xarray[i]=%e+i*%e;\n\n", dim, xMin, dX);
	    }
	  f_printf (rf, "// GRAPH %d\n", first);

	  while ((getc (infile) != 'Y') ||
		 (getc (infile) != 'a') ||
		 (getc (infile) != 'r') ||
		 (getc (infile) != 'r') ||
		 (getc (infile) != 'a') ||
		 (getc (infile) != 'y'))
	    {
	    }
	  while ((c = getc (infile)) != '{')
	    {
	    }

	  f_printf (rf, "double Yarray%d[%d]={", first, dim);
	  while ((c = getc (infile)) != '}')
	      putc (c, rf);
	  fclose (infile);

	  f_printf (rf, "};\n");


	  first++;
	}
    }


          f_printf (rf, "for(int i=0;i<%d;i++) Yarray1[i] = \n", dim);
         
          for(i=1;i<first;i++) 
            { 
              if (i==1) f_printf (rf, " Yarray%d[i]",i);
              else      f_printf (rf, " + Yarray%d[i]",i);
              if (0 == (i + 1) % 5) f_printf (rf, "\n");
            }         
         f_printf (rf, ";\n\n",i);

         
         f_printf (rf, "// printf(\" Yarray1[%d]={\")\n",dim);
         f_printf (rf, "// for(i=0;i<%d;i++) printf(\" e,\",Yarray1[i]);\n",dim);
         f_printf (rf, "// printf(\" } \") \n\n"); 

	  f_printf (rf, "TGraph *gr1 = new ");
	  f_printf (rf, "TGraph(%d,Xarray,Yarray1);\n", dim);
	  if (first == 1)
	    f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
	  f_printf (rf, "gr1->SetLineWidth(2);\n");
/*	  color = graphcolor (first);   */
	  f_printf (rf, "gr1->SetLineColor(4);\n");
	  style = graphstyle (first);
	  f_printf (rf, "gr1->SetLineStyle(1);\n");  
	    
	      sprintf (xaxistitle, "\"%s\"", x_str);
	      sprintf (yaxistitle, "\"#sigma, fb\"");
	      if (!strncmp (x_str, "//Transverse momentum Pt%d", 23))
		{
		  sprintf (xaxistitle, "\"P_{T}, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d P_{T}, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Mass", 6))
		{
		  sprintf (xaxistitle, "\"M, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d M, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Energy E", 10))
		{
		  sprintf (xaxistitle, "\"E, GeV\"");
		  sprintf (yaxistitle, "\"d#sigma/d E, fb/%.0f GeV\"", dX);
		}

	      if (!strncmp (x_str, "//Angle", 7))
		{
		  sprintf (xaxistitle, "\"#theta, deg\"");
		  sprintf (yaxistitle, "\"d#sigma/d #theta, fb/%.0f deg\"", dX);
		}

	      if (!strncmp (x_str, "//pseudo-rapidity", 17))
		{
		  sprintf (xaxistitle, "\"#eta\"");
		  sprintf (yaxistitle, "\"d#sigma/d #eta, fb\"");
		}

	      if (!strncmp (x_str, "//Rapidity", 10))
		{
		  sprintf (xaxistitle, "\"#y\"");
		  sprintf (yaxistitle, "\"d#sigma/d #y, fb\"");
		}

	      if (!strncmp (x_str, "//Cosine", 8))
		{
		  sprintf (xaxistitle, "\"cos(#theta)\"");
		  sprintf (yaxistitle, "\"d#sigma/d cos(#theta), fb\"");
		}
	      f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
	      f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);

	      f_printf (rf, "//gr1->SetMinimum(%e);\n", totYMin);
	      f_printf (rf, "//gr1->SetMaximum(%e);\n", 1.1 * totYMax);
	      f_printf (rf, "gr1->Draw(\"AC\");\n\n");
	    


  f_printf (rf, "//c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "//c1->Print(\"%s.png\");\n}\n", outname);
  fclose (rf);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


static void
writeroothist3d (char *menustr, int numstring)
{
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  midstr outname_tmp;
  FILE *rf1;
  long POS;
  int flag=0;
  double x,y,W,x_comb,y_comb,W_comb;

  nextFileName (outname, "sum2d_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   fclose (rf1);

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
          rf1 = fopen (outname_C, "r+");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;
              if(flag==0)  fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W);    
              else
                {
                  POS=ftell(rf1);
                  fscanf(rf1,"%le %le %le ",&x_comb,&y_comb,&W_comb); 
                  fseek(rf1,POS,SEEK_SET);
                  W_comb += W;
                  fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W_comb);  
                }
             }
          flag=1;
          
          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}

static void
writeroothist1d (char *menustr, int numstring)
{
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  midstr outname_tmp;
  FILE *rf1;
  long POS;
  int flag=0;
  double x,y,W,x_comb,y_comb,W_comb;

  nextFileName (outname, "sum1d_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   fclose (rf1);

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
          rf1 = fopen (outname_C, "r+");
         
	  while(fscanf(infile,"%le %le ",&x,&W)!=EOF)  
             {
              if(W!=W) W=0;
              if(flag==0)  fprintf(rf1,"% .6E % .6E \n",x,W);    
              else
                {
                  POS=ftell(rf1);
                  fscanf(rf1,"%le %le ",&x_comb,&W_comb); 
                  fseek(rf1,POS,SEEK_SET);
                  W_comb += W;
                  fprintf(rf1,"% .6E % .6E \n",x,W_comb);  
                }
             }
          flag=1;
          
          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


static void
multiplyroothist3d (char *menustr, int numstring)
{
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  midstr outname_tmp;
  FILE *rf1;
  long POS;
  int flag=0;
  double x,y,W,x_comb,y_comb,W_comb;

  nextFileName (outname, "prod2d_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   fclose (rf1);

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
          rf1 = fopen (outname_C, "r+");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;
              if(flag==0)  fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W);    
              else
                {
                  POS=ftell(rf1);
                  fscanf(rf1,"%le %le %le ",&x_comb,&y_comb,&W_comb); 
                  fseek(rf1,POS,SEEK_SET);
                  W_comb *= W;
                  fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W_comb); 
                }
             }
          flag=1;
          
          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


static void
write_sumHistNum (char *menustr, int numstring)
{
  void *pscr0 = NULL;
  int i, k;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1, *infile;
  double x,y,W,W2, koeff=0.0;  



label1:
  if (!correctDouble (45, 14, "koefficient = ", &koeff, 0))
    {
      goto exi;
    }
 
  goto_xy (45, 14);
  clr_eol ();


  nextFileName (outname, "mult_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  infile = fopen (infname, "r");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;

              W2=W+koeff;
               
              fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W2);          
             }

          fclose (rf1);
	}
    }

  fclose (infile);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));


exi:
  put_text (&pscr0);
  goto_xy (1, 23);
  scrcolor (FGmain, BGmain);
  clr_eol ();
}


static void
writeHistNum (char *menustr, int numstring)
{
  void *pscr0 = NULL;
  int i, k;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1, *infile;
  double x,y,W,W2, koeff=1.0;  



label1:
  if (!correctDouble (45, 14, "koefficient = ", &koeff, 0))
    {
      goto exi;
    }
 
  goto_xy (45, 14);
  clr_eol ();


  nextFileName (outname, "mult_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  infile = fopen (infname, "r");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;

              W2=W*koeff;
               
              fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W2);          
             }

          fclose (rf1);
	}
    }

  fclose (infile);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));


exi:
  put_text (&pscr0);
  goto_xy (1, 23);
  scrcolor (FGmain, BGmain);
  clr_eol ();
}


static void
writeRescaleXY (char *menustr, int numstring, int xy)
{
  void *pscr0 = NULL;
  int i, k;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1, *infile;
  double x,y,W,W2, koeff=1.0, delta=0.0;  


if(xy == 1)
{
label1:
  if (!correctDouble (45, 14, "koeff x = ", &koeff, 0))
    {
      goto exi;
    }
  if (!correctDouble (45, 15, "delta x = ", &delta, 0))
    {
      goto_xy (45, 14);
      clr_eol ();
      goto label1;
    }
}


if(xy == 2)
{
label2:
  if (!correctDouble (45, 14, "koeff y = ", &koeff, 0))
    {
      goto exi;
    }
  if (!correctDouble (45, 15, "delta y = ", &delta, 0))
    {
      goto_xy (45, 14);
      clr_eol ();
      goto label2;
    }
}

  goto_xy (45, 14);
  clr_eol ();



  if(xy == 1) nextFileName (outname, "rescalex_", ".txt");
  if(xy == 2) nextFileName (outname, "rescaley_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   

if (xy == 1)
{
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  infile = fopen (infname, "r");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;

              x = x*koeff + delta;
               
              fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W);          
             }

          fclose (rf1);
	}
    }

  fclose (infile);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


if (xy == 2)
{
  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  infile = fopen (infname, "r");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;

              y = y*koeff + delta;
               
              fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W);          
             }

          fclose (rf1);
	}
    }

  fclose (infile);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



exi:
  put_text (&pscr0);
  goto_xy (1, 23);
  scrcolor (FGmain, BGmain);
  clr_eol ();
}


static void
divideroothist3d (char *infname1, char *infname2)
{
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1, *infile1, *infile2;
  double x1,y1,W1,x2,y2,W2;

  nextFileName (outname, "ratio2d_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   fclose (rf1);


    
	  infile1 = fopen (infname1, "r");
          infile2 = fopen (infname2, "r");
          rf1 = fopen (outname_C, "r+");
         
	  while(fscanf(infile1,"%le %le %le ",&x1,&y1,&W1)!=EOF)  
             {
              if(W1!=W1) W1=0;
            
              fscanf(infile2,"%le %le %le",&x2,&y2,&W2);
                                   
              W1 = W1/W2;
              fprintf(rf1,"% .6E % .6E % .6E\n",x1,y1,W1);     
             }
    
          
          fclose (rf1);
	  fclose (infile1);
          fclose (infile2);
	
  

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


static void
writeroothist3dhi2 (char *menustr, int numstring)
{
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  midstr outname_tmp;
  FILE *rf1;
  long POS;
  int flag=0;
  double x,y,W,x_comb,y_comb,W_comb;

  nextFileName (outname, "sumhi2_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   fclose (rf1);

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");
          rf1 = fopen (outname_C, "r+");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;
              if(flag==0)  fprintf(rf1,"%e %e %e \n",x,y,W);    
              else
                {
                  POS=ftell(rf1);
                  fscanf(rf1,"%le %le %le ",&x_comb,&y_comb,&W_comb); 
                  fseek(rf1,POS,SEEK_SET);
                  W_comb += W;
                  fprintf(rf1,"% .6E % .6E % .6E \n",x,y,W_comb);
                }
             }
          flag=1;
          
          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



static void
writeroothi2 (char *menustr, int numstring)
{
  void *pscr0 = NULL;
  int i, k;
  int combine = 0;
  int dim;
  char c;

  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  midstr outname_tmp;
  FILE *rf1, *infile;
  double x,y,W;
  double hi, hi2, sigma_sm=1.0, data=1.0, data_error=1.0;  



label1:
  if (!correctDouble (45, 14, "SM cross section = ", &sigma_sm, 0))
    {
      goto exi;
    }
  if (sigma_sm == 0)
    {
      warnanykey (55, 17, "Range check error");
      goto_xy (55, 15);
      clr_eol ();
      sigma_sm = 1.0;
      goto label1;
    } 
label2:
  if (!correctDouble (45, 15, "Experimental signal strength = ", &data, 0))
    {
      goto_xy (45, 14);
      clr_eol ();
      goto label1;
    }

label3:
  if (!correctDouble (45, 16, "Signal strength error = ", &data_error, 0))
    {
      goto_xy (45, 14);
      clr_eol ();
      goto label1;
    }
  if (data_error == 0)
    {
      warnanykey (55, 17, "Range check error");
      goto_xy (45, 15);
      clr_eol ();
      data_error = 1.0;
      goto label1;
    }  


  goto_xy (45, 14);
  clr_eol ();
  goto_xy (45, 15);
  clr_eol ();
  goto_xy (45, 16);
  clr_eol ();


/*  informline (0, npoints); */



  nextFileName (outname, "hi2_", ".txt");
  strcpy (outname_C, outname);
  strcat (outname_C, ".txt");

  rf1 = fopen (outname_C, "w");
  if (NULL == rf1)
    {
      messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
      return;
    }
   

  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  infile = fopen (infname, "r");
         
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              if(W!=W) W=0;

              hi = (W/sigma_sm - data)/data_error;
              hi2 = hi*hi;
               
              fprintf(rf1,"% .6E % .6E % .6E \n",x,y,hi2);          
             }

          fclose (rf1);
	}
    }

  fclose (infile);
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));


exi:
  put_text (&pscr0);
  goto_xy (1, 23);
  scrcolor (FGmain, BGmain);
  clr_eol ();
}



static void
write_root3dsurface (char *menustr, int numstring)
{
  int i, dim;
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1;
  long POS;
  int count;
  double x,y,W;


  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");

          nextFileName (outname, "surface_", ".C");
          strcpy (outname_C, outname);
          strcat (outname_C, ".C");
          rf1 = fopen (outname_C, "w+");
          if (NULL == rf1)
           {
            messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
            return;
           }
  
           fprintf(rf1,"#include <TMath.h>\n");
           fprintf(rf1,"#include <TGraph2D.h>\n");
           fprintf(rf1,"#include <TRandom.h>\n");
           fprintf(rf1,"#include <TStyle.h>\n");
           fprintf(rf1,"#include <TCanvas.h>\n");
           fprintf(rf1,"#include <TF2.h>\n");
           fprintf(rf1,"#include <TH1.h>\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"void %s()\n",outname);
           fprintf(rf1,"{\n");
           fprintf(rf1,"\n");

           POS=ftell(rf1);           
           fprintf(rf1,"double x[9],y[9],z[9];          \n"); 
            
           fprintf(rf1,"\n");

          count=0; 
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              fprintf(rf1,"x[%d]=% .6E; y[%d]=% .6E; z[%d]=% .6E;\n",count,x,count,y,count,W);    
              count+=1;
             }

           fprintf(rf1,"\n");
           fprintf(rf1,"for(int i=0;i<%d;i++)  z[i]=1000*z[i];\n",count);
           fprintf(rf1,"\n"); 
           fprintf(rf1,"TCanvas *c = new TCanvas(\"c\",\"Graph2D\",600,0,600,600);\n");
           fprintf(rf1,"TGraph2D *dt = new TGraph2D(%d,x,y,z);\n",count);
           fprintf(rf1,"\n");
           fprintf(rf1,"gPad->SetLogz();\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"// gPad->SetTheta(25);\n");
           fprintf(rf1,"// gPad->SetPhi(-110);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->SetTitle(\"\");\n");
           fprintf(rf1,"dt->SetLineWidth(2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitle(\"par1\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitleOffset(2.2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetLabelSize(0.025);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitle(\"par2\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitleOffset(1.9);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetLabelSize(0.025);\n");
           fprintf(rf1,"dt->GetHistogram()->GetZaxis()->SetTitle(\"#sigma, fb\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetZaxis()->SetTitleOffset(1.4);\n");
           fprintf(rf1,"dt->GetHistogram()->GetZaxis()->SetLabelSize(0.03);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->Draw(\"surf1\");\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"c->Print(\"%s.eps\");\n",outname);
           fprintf(rf1,"c->Print(\"%s.png\");\n",outname);
           fprintf(rf1,"}\n");

           fseek(rf1,POS,SEEK_SET);
           fprintf(rf1,"double x[%d],y[%d],z[%d];\n",count,count,count); 

          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}


static void
write_root1dtable (char *menustr, int numstring)
{
  int i, dim;
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1;
  long POS;
  int count;
  double x,y,W;


  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");

          nextFileName (outname, "table1d_", ".C");
          strcpy (outname_C, outname);
          strcat (outname_C, ".C");
          rf1 = fopen (outname_C, "w+");
          if (NULL == rf1)
           {
            messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
            return;
           }
/*  
           fprintf(rf1,"#include <TMath.h>\n");
           fprintf(rf1,"#include <TGraph2D.h>\n");
           fprintf(rf1,"#include <TRandom.h>\n");
           fprintf(rf1,"#include <TStyle.h>\n");
           fprintf(rf1,"#include <TCanvas.h>\n");
           fprintf(rf1,"#include <TF2.h>\n");
           fprintf(rf1,"#include <TH1.h>\n");
           fprintf(rf1,"\n");
*/
           fprintf(rf1,"void %s()\n",outname);
           fprintf(rf1,"{\n");
           fprintf(rf1,"\n");

           POS=ftell(rf1);           
           fprintf(rf1,"double x[9],y[9];          \n"); 
            
           fprintf(rf1,"\n");

          count=0; 
	  while(fscanf(infile,"%le %le ",&x,&W)!=EOF)  
             {
              fprintf(rf1,"x[%d]=% .6E; y[%d]=% .6E;\n",count,x,count,W);    
              count+=1;
             }

           fprintf(rf1,"\n");
           fprintf(rf1,"for(int i=0;i<%d;i++)  y[i]=1000*y[i];\n",count);
           fprintf(rf1,"\n"); 

           fprintf(rf1,"TCanvas *c = new TCanvas(\"c\",\"Graph\",200,100,700,500);\n");
           fprintf(rf1,"TGaxis::SetMaxDigits(3);\n");
           fprintf(rf1, "c->SetFillColor(0);\n");
           fprintf(rf1, "c->SetBorderMode(0);\n");
           fprintf(rf1, "c->SetFrameBorderMode(0);\n");
           fprintf(rf1, "gStyle->SetOptTitle(kFALSE);\n\n");
           fprintf(rf1, "//c->SetLogy(1);  // Set Log. Scale;\n\n\n");
           fprintf(rf1,"TGraph *dt = new TGraph(%d,x,y);\n",count);
           
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->SetTitle(\"\");\n");
           fprintf(rf1,"dt->SetLineWidth(2);\n");
           fprintf(rf1,"dt->SetLineColor(4);\n");
           fprintf(rf1,"dt->SetLineStyle(1);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitle(\"par1\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitleOffset(2.2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetLabelSize(0.025);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitle(\"#sigma, fb\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetLabelSize(0.03);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->Draw(\"AC\");\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"c->Print(\"%s.eps\");\n",outname);
           fprintf(rf1,"c->Print(\"%s.png\");\n",outname);
           fprintf(rf1,"}\n");

           fseek(rf1,POS,SEEK_SET);
           fprintf(rf1,"double x[%d],y[%d];\n",count,count); 

          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



static void
write_constract2d (char *menustr, int numstring)
{
  int i,j, k, dim;
  char infname[STRSIZ],currentinfname[STRSIZ],c;
  char upstr[50] = { 0 };
  char x_str[50] = { 0 };
  char y_str[50] = { 0 };
  midstr outname;
  midstr outname_C;
  FILE *rf1,*infile,*rootfile;
  long POS;
  int count,combine;
  double x,x1,y,W,xMin,xMax,yMin,yMax,delta_x;


  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");

          nextFileName (outname, "histtab_", ".txt");
          strcpy (outname_C, outname);
          strcat (outname_C, ".txt");
          rf1 = fopen (outname_C, "w+");
          if (NULL == rf1)
           {
            messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
            return;
           }

 
          count=1; 
	  while(fscanf(infile,"%le %le ",&x,&W)!=EOF)  
            {

              sprintf(currentinfname,"root_%d.C",count);

    	      FILE *rootfile = fopen (currentinfname, "r");
              if (rootfile == NULL)
               {
                messanykey (10, 12, scat ("Error! I can not open the file\n%s", currentinfname));
                return;
               }
	      k = 0;
	      while ((c = getc (rootfile)) != '\n')
	       {
	         upstr[k] = c;
	         k++;
	       }
	      k = 0;
	      while ((c = getc (rootfile)) != '\n')
	       {
	         x_str[k] = c;
	         k++;
	       }
	      fscanf (rootfile, "//X from %le to %le N_bins= %d\n", &xMin, &xMax, &dim);
	      k = 0;
	      while ((c = getc (rootfile)) != '\n')
	       {
	         y_str[k] = c;
	         k++;
	       }
	      fscanf (rootfile, "//Y from %le to %le", &yMin, &yMax);
	      fscanf (rootfile, "//combine= %d", &combine);

	      while ((getc (rootfile) != 'Y') ||
		     (getc (rootfile) != 'a') ||
		     (getc (rootfile) != 'r') ||
		     (getc (rootfile) != 'r') ||
		     (getc (rootfile) != 'a') ||
		     (getc (rootfile) != 'y'))
	        {
	        }
	      while ((c = getc (rootfile)) != '{')
	        {
	        }
              
              delta_x= (xMax-xMin)/(dim-1);
  
              for(j=0;j<dim;j++) 
                {
                 if(j<dim-1) fscanf (rootfile, "%le, ", &y);
                 else        fscanf (rootfile, "%le", &y);  
                 x1 = xMin + j*delta_x;
                 fprintf(rf1,"% .6E % .6E % .6E \n",x,x1,y);
                }

	      fclose (rootfile);    
              count+=1;
             }

          fclose (rf1);
	  fclose (infile);
	}
    }
  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



static void
write_rootCounturPlot (char *menustr, int numstring)
{
  int i, dim;
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1;
  long POS;
  int count;
  double x,y,W;
  double min_z,mid_z,max_z;


  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");

          nextFileName (outname, "contour_", ".C");
          strcpy (outname_C, outname);
          strcat (outname_C, ".C");
          rf1 = fopen (outname_C, "w");
          if (NULL == rf1)
           {
            messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
            return;
           }
  
           fprintf(rf1,"#include <TMath.h>\n");
           fprintf(rf1,"#include <TGraph2D.h>\n");
           fprintf(rf1,"#include <TRandom.h>\n");
           fprintf(rf1,"#include <TStyle.h>\n");
           fprintf(rf1,"#include <TCanvas.h>\n");
           fprintf(rf1,"#include <TF2.h>\n");
           fprintf(rf1,"#include <TH1.h>\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"void %s()\n",outname);
           fprintf(rf1,"{\n");
           fprintf(rf1,"\n");

           POS=ftell(rf1);           
           fprintf(rf1,"double x[9],y[9],z[9];          \n"); 
              
           fprintf(rf1,"\n");
          
          min_z=0;
          max_z=0; 

          count=0; 
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              fprintf(rf1,"x[%d]=% .6E; y[%d]=% .6E; z[%d]=% .6E;\n",count,x,count,y,count,W);    
              count+=1;
              
              if(count==1) min_z=W;
              if(W > max_z) max_z=W;
              if(W < min_z) min_z=W;
             }

           mid_z=(max_z+min_z)/2.0;

           fprintf(rf1,"\n"); 
           fprintf(rf1,"TCanvas *c = new TCanvas(\"c\",\"Graph2D\",600,0,600,600);\n");
           fprintf(rf1,"TGraph2D *dt = new TGraph2D(%d,x,y,z);\n",count);
           fprintf(rf1,"\n");
           fprintf(rf1,"Double_t contours[5]={%e, %e, %e, %e, %e};\n",min_z, (mid_z+min_z)/2.0, mid_z, (max_z+mid_z)/2.0, max_z);
           fprintf(rf1,"dt->GetHistogram()->SetContour(5,contours);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->SetTitle(\"\");\n");
           fprintf(rf1,"dt->SetLineWidth(2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitle(\"par1\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitle(\"par2\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);\n");
         
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->Draw(\"CONT4\");\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"c->Print(\"%s.eps\");\n",outname);
           fprintf(rf1,"c->Print(\"%s.png\");\n",outname);
           fprintf(rf1,"}\n");

           fseek(rf1,POS,SEEK_SET);
           fprintf(rf1,"double x[%d],y[%d],z[%d];\n",count,count,count); 

          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}




static void
write_rootDeltaHi2 (char *menustr, int numstring)
{
  int i, dim;
  char infname[STRSIZ];
  midstr outname;
  midstr outname_C;
  FILE *rf1;
  long POS;
  int count;
  double x,y,W;
  double min_z,mid_z,max_z;


  for (i = 1; i <= numstring; i++)
    {
      sprintf (infname, "%.14s", menustr + i * 15 - 14);
      trim (infname);
      if (menustr[i * 15] == '+')
	{
	  FILE *infile = fopen (infname, "r");

          nextFileName (outname, "hi2fit_", ".C");
          strcpy (outname_C, outname);
          strcat (outname_C, ".C");
          rf1 = fopen (outname_C, "w");
          if (NULL == rf1)
           {
            messanykey (10, 12, scat ("Error! I can not open the file\n%s", outname_C));
            return;
           }
  
           fprintf(rf1,"#include <TMath.h>\n");
           fprintf(rf1,"#include <TGraph2D.h>\n");
           fprintf(rf1,"#include <TRandom.h>\n");
           fprintf(rf1,"#include <TStyle.h>\n");
           fprintf(rf1,"#include <TCanvas.h>\n");
           fprintf(rf1,"#include <TF2.h>\n");
           fprintf(rf1,"#include <TH1.h>\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"void %s()\n",outname);
           fprintf(rf1,"{\n");
           fprintf(rf1,"\n");

           POS=ftell(rf1);           
           fprintf(rf1,"double x[9],y[9],z[9];          \n"); 
              
           fprintf(rf1,"\n");
          
          min_z=0;
          max_z=0; 

          count=0; 
	  while(fscanf(infile,"%le %le %le ",&x,&y,&W)!=EOF)  
             {
              fprintf(rf1,"x[%d]=% .6E; y[%d]=% .6E; z[%d]=% .6E;\n",count,x,count,y,count,W);    
              count+=1;
              
              if(count==1) min_z=W;
              if(W > max_z) max_z=W;
              if(W < min_z) min_z=W;
             }

           mid_z=(max_z+min_z)/2.0;

           fprintf(rf1,"\n"); 
           fprintf(rf1,"TCanvas *c = new TCanvas(\"c\",\"Graph2D\",600,0,600,600);\n");
           fprintf(rf1,"TGraph2D *dt = new TGraph2D(%d,x,y,z);\n",count);
           fprintf(rf1,"\n");
           fprintf(rf1,"Double_t contours[4]={%e, %e, %e, %e};\n",min_z, min_z+2.1, min_z+4.61, min_z+9.21);
           fprintf(rf1,"dt->GetHistogram()->SetContour(4,contours);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"Int_t colors[4]={4,3,5,0};\n");
           fprintf(rf1,"gStyle->SetPalette(4,colors);\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->SetTitle(\"\");\n");
           fprintf(rf1,"dt->SetLineWidth(2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitle(\"par1\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitle(\"par2\");\n");
           fprintf(rf1,"dt->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);\n");
         
           fprintf(rf1,"\n");
           fprintf(rf1,"dt->Draw(\"CONT4\");\n");
           fprintf(rf1,"\n");
           fprintf(rf1,"c->Print(\"%s.eps\");\n",outname);
           fprintf(rf1,"c->Print(\"%s.png\");\n",outname);
           fprintf(rf1,"}\n");

           fseek(rf1,POS,SEEK_SET);
           fprintf(rf1,"double x[%d],y[%d],z[%d];\n",count,count,count); 

          fclose (rf1);
	  fclose (infile);
	}
    }

  messanykey (10, 12, scat (" You can find results in the file\n%s", outname_C));
}



void
combine_root (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.C");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	menustr[k++] = ' ';
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeroothist (menustr, numstring);
}

void
combine_rootSum (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.C");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
          if(i==15) menustr[k++] = '+';
          else      menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeroothist_sum (menustr, numstring);
}



void
combine_Sum3d (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
          if(i==15) menustr[k++] = '+';
          else      menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeroothist3d (menustr, numstring);
}

void
combine_Sum1d (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
          if(i==15) menustr[k++] = '+';
          else      menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeroothist1d (menustr, numstring);
}



void
multiply_Sum3d (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
          menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  multiplyroothist3d (menustr, numstring);
}


void
multiply_HistNum (void)
{
  int i,j, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeHistNum (menustr, numstring);
}


void
rescale_xy (int xy)
{
  int i,j, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);

  writeRescaleXY (menustr, numstring, xy);
}



void
sum_HistNum (void)
{
  int i,j, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  write_sumHistNum (menustr, numstring);
}



void
divide_Sum3d (void)
{
  int i,j, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char infname1[STRSIZ];
  char infname2[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

 messanykey (15, 15,"Enter dividend");

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';

            sprintf (infname1, "%.14s", f_name);
            trim (infname1);
            printf ("\n infname1 = %s\n",infname1);
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);

  for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 

  messanykey (15, 15,"Enter divider");
 

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';

            sprintf (infname2, "%.14s", f_name);
            trim (infname2);
            printf ("\n infname2 = %s\n",infname2);
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);


  divideroothist3d (infname1, infname2); 
}


void
combine_hi2 (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
          menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);

  writeroothist3dhi2 (menustr, numstring);
}


void
plot_root3d (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  write_root3dsurface (menustr, numstring);
}


void
plot_root1d (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);

  write_root1dtable (menustr, numstring);
}

void
constract_2d (void)
{ 
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);

  write_constract2d (menustr, numstring);
}

void
plot_rootCountur (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  write_rootCounturPlot (menustr, numstring);
}


void
calc_hi2 (void)
{
  int i,j, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
           {
            for(j=0;j<20016;j++) if(menustr[j] == '+') menustr[j] = ' '; 
	    menustr[k * 15] = '+';
           }
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  writeroothi2 (menustr, numstring);
}


void
plot_deltaHi2 (void)
{
  int i, k, numstring = 0;
  int doserror;
  char f_name[STRSIZ];
  char menustr[20016];
  void *pscr2 = NULL;
  searchrec s;

  menustr[0] = 15;

  k = 1;
  strcpy (f_name, "*.txt");
  doserror = find_first (f_name, &s, anyfile);
  while (doserror == 0 && k <= 20000)
    {
      for (i = 1; (i <= strlen (s.name)) && (i <= 15); i++)
	menustr[k++] = s.name[i - 1];
      for (i = strlen (s.name) + 1; i <= 15; i++)
	{ 
         menustr[k++] = ' ';
        } 
      doserror = find_next (&s);
      numstring++;
    }
  menustr[k] = 0;
  if (0 == menustr[1])
    {
      messanykey (10, 15, "directory  is empty");
      return;
    }

  k = 1;
  do
    {
      menu1 (10, 10, "", menustr, "", &pscr2, &k);
      if (k > 0)
	{
	  sprintf (f_name, "%.14s", menustr + k * 15 - 14);
	  readtext (f_name);
	  if (menustr[k * 15] == ' ')
	    menustr[k * 15] = '+';
	  else
	    menustr[k * 15] = ' ';
	}
    }
  while (k != 0);
  write_rootDeltaHi2 (menustr, numstring);
}

