/* 
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Victor Edneral
* ------------------------------------------------------
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"
#include "service2/include/tptcmac.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt.h"
#include "chep_crt/include/crt_util.h"
#include "chep_crt/include/tex_util.h"
#include "num/include/spline.h"

#include "plot.h"

static int X1, Y1, X2, Y2;

static struct
{
  double xmin, xmax, ymin, ymax;
}
grafminmax;
static double xscale, yscale;


static int pictureX = 300;
static int pictureY = 200;
static char letterSize[14] = "normalize";

static int logScale = TRUE;
static int fgcolor = FGmain, bkcolor = BGmain;



static int
nlog10 (double x)
{
  double lg = log10 (x);
  if (lg < 0)
    return lg - 1;
  else
    return lg;
}


static void
axisDesign (double xmin, double xmax, int islog,
	    double *xfirst, double *step, int *nsub)
{
  double dx, dx0, dx100;
  int n, n0, n1;
  char xmintxt[50];

  if (islog)
    {
      n = 1 + log10 (xmax / xmin) / 10;
      *step = pow ((double) 10, (double) n);
      *xfirst = pow ((double) 10, (double) nlog10 (xmin));
      *nsub = 9;
    }
  else
    {
      dx = xmax - xmin;
      n = nlog10 (dx);
      dx0 = pow ((double) 10, (double) n);

      dx100 = 10 * dx / dx0;
      if (dx100 < 15.0)
	{
	  *step = 0.3 * dx0;
	  *nsub = 3;
	}
      else if (dx100 < 24.0)
	{
	  *step = 0.5 * dx0;
	  *nsub = 5;
	}
      else if (dx100 < 30.0)
	{
	  *step = 0.8 * dx0;
	  *nsub = 8;
	}
      else if (dx100 < 45.0)
	{
	  *step = 1 * dx0;
	  *nsub = 10;
	}
      else if (dx100 < 90.0)
	{
	  *step = 2 * dx0;
	  *nsub = 2;
	}
      else
	{
	  *step = 3 * dx0;
	  *nsub = 3;
	}
      if (fabs (xmin) <= (*step) * 10)
	*xfirst = 0;
      else
	{
	  n0 = nlog10 (*step);
	  n1 = nlog10 (fabs (xmin));

	  sprintf (xmintxt, "%.*E", MAX (0, n1 - n0 - 1), xmin);
	  trim (xmintxt);
	  sscanf (xmintxt, "%lf", xfirst);
	}
      while (*xfirst > xmin + (*step) / (*nsub))
	*xfirst -= *step;
      while (*xfirst + (*step) < xmin)
	*xfirst += *step;
    }
}

static long
chpround (double x)
{
  if (x > 0)
    x += 0.5;
  else
    x -= 0.5;
  return (long) x;
}

static void
gminmax (double *f, double *df, int dim)
{
  int i;
  if (dim == 0)
    return;

  grafminmax.ymin = f[0];
  grafminmax.ymax = f[0];
  if (df)
    {
      grafminmax.ymax += df[0];
      grafminmax.ymin -= df[0];
    }

  if (df)
    for (i = 1; i < dim; i++)
      {
	grafminmax.ymin = MIN (grafminmax.ymin, f[i] - df[i]);
	grafminmax.ymax = MAX (grafminmax.ymax, f[i] + df[i]);
      }
  else
    for (i = 1; i < dim; i++)
      {
	grafminmax.ymin = MIN (grafminmax.ymin, f[i]);
	grafminmax.ymax = MAX (grafminmax.ymax, f[i]);
      }
}


static int
scx (double x)
{
  x -= grafminmax.xmin;
  return X1 + chpround (xscale * x);
}


static int
scy (double x)
{
  if (logScale)
    x = log10 (x / grafminmax.ymin);
  else
    x -= grafminmax.ymin;
  return Y1 - chpround (yscale * x);
}


static double
dscx (double x)
{
  return X1 + (xscale * (x - grafminmax.xmin));
}


static double
dscy (double x)
{
  if (logScale)
    {
      if (x > 0)
	x = log10 (x / grafminmax.ymin);
      else
	x = -10000 * log10 (grafminmax.ymax / grafminmax.ymin);
    }
  else
    x -= grafminmax.ymin;
  return Y1 - (yscale * x);
}


static void
doubleToStr (double x, double step, char *s)
{
  int n1, n0;
  char *emark;
  char s1[128], s2[128], sgn[16];
  float f;

  if (fabs (x) < 1.E-2 * step)
    strcpy (s, "0.0");
  else
    {
      n1 = nlog10 (fabs (x));
      n0 = nlog10 (step);
      sprintf (s1, "%.*E", MAX (0, n1 - n0), x);

      emark = strstr (s1, "E");

      emark[0] = 0;
      emark++;
      sgn[0] = 0;
      switch (emark[0])
	{
	case '-':
	  strcpy (sgn, "-");
	  emark++;
	  break;
	case '+':
	  sgn[0] = 0;
	  emark++;
	  break;
	}
      while (emark[0] == '0')
	emark++;
      if (emark[0] == 0)
	strcpy (s2, s1);
      else
	{
	  sscanf (s1, "%f", &f);
	  if (f == 1.0)
	    sprintf (s2, "10^%s%s", sgn, emark);
	  else
	    sprintf (s2, "%s*10^%s%s", s1, sgn, emark);
	}

      sprintf (s1, "%.*f", MAX (0, -n0), x);

      if (strlen (s1) < strlen (s2))
	strcpy (s, s1);
      else
	strcpy (s, s2);
    }
}


static void
gaxes (char *upstr, char *xstr, char *ystr)
{
  double xmax, xmin, ymax, ymin, aa, step;
  int m, naxis, islog, n_sub;

  int th = tg_textheight ("0");
  int tw = tg_textwidth ("0");
  int hash = (th + tw) / 2;
  int texhash = (th * texyscale + tw * texxscale + 1) / 2;
  xmin = grafminmax.xmin;
  xmax = grafminmax.xmax;
  ymin = grafminmax.ymin;
  ymax = grafminmax.ymax;
  if (1.E-4 * (fabs (ymax) + fabs (ymin)) >= fabs (ymax - ymin))
    {
      if (ymin == 0.0)
	{
	  ymin = -0.5;
	  ymax = 0.5;
	}
      else if (ymin < 0.0)
	{
	  ymin *= 2.0;
	  if (ymax < 0.0)
	    ymax *= 0.5;
	  else
	    ymax *= 2.0;
	}
      else
	{
	  ymax *= 2.0;
	  ymin *= 0.5;
	}
    }
  grafminmax.ymin = ymin;
  grafminmax.ymax = ymax;

  if (logScale && (ymin <= 0.0 || log10 (ymax / ymin) < 1))
    logScale = FALSE;

  tg_setlinestyle (SolidLn, NormWidth);

  X1 = 10 * tw;
  X2 = tg_getmaxx () - 2 * tw;
  Y1 = tg_getmaxy () - 5 * th;
  Y2 = 5 * th / 2;

  xscale = (X2 - X1) / (xmax - xmin);
  yscale = (Y1 - Y2) / (logScale ? log10 (ymax / ymin) : ymax - ymin);
  tg_settextjustify (CenterText, TopText);
  tg_outtextxy ((X1 + X2) / 2, CenterText, upstr);

  for (naxis = 0; naxis < 2; naxis++)
    {
      double zmin, zmax, zmin1, zmax1;
      char xy;
      int xend, yend;
      int len;


      if (naxis == 0)
	{
	  xy = 'X';
	  if (texflag)
	    {
	      hash = texhash / texyscale;
	      texhash = -texhash;
	    }
	  islog = FALSE;
	  tg_settextjustify (CenterText, TopText);
	  zmin = xmin;
	  zmax = xmax;
	  xend = X2;
	  yend = Y1;
	}
      else
	{
	  xy = 'Y';
	  if (texflag)
	    {
	      texhash = -texhash;
	      hash = texhash / texxscale;
	    }
	  islog = logScale;
	  tg_settextjustify (RightText, CenterText);
	  zmin = ymin;
	  zmax = ymax;
	  xend = X1;
	  yend = Y2;
	}

      if (texflag)
	f_printf (out_tex,
		  "%% ====================   %c-axis =============\n", xy);

      zmax1 = zmax + fabs (zmax) * 1.E-7;
      zmin1 = zmin - fabs (zmin) * 1.E-7;
      axisDesign (zmin, zmax, islog, &aa, &step, &n_sub);
      if (texflag)
	{
	  double Nd, offset;
	  char axis[10];
	  char d[10];

	  if (islog)
	    {
	      strcpy (axis, "LogAxis");
	      Nd = log10 (zmax / zmin) / log10 (step);
	      offset = 10 * zmin / (step * aa);
	      strcpy (d, "");
	    }
	  else
	    {
	      strcpy (axis, "LinAxis");
	      Nd = (zmax - zmin) / step;
	      offset = -n_sub * (aa - zmin) / step;
	      sprintf (d, "%d,", n_sub);
	      if (fabs (offset) > fabs (offset - n_sub))
		offset -= n_sub;
	    }

	  f_printf (out_tex,
		    "\\%s(%.2f,%.2f)(%.2f,%.2f)(%.3f,%s%d,%.3f,1.5)\n", axis,
		    texX (X1), texY (Y1), texX (xend), texY (yend), Nd, d,
		    (2 * texhash) / 3, offset);
	}
      else
	tg_line (X1, Y1, xend, yend);
      len = 0;
      while (aa <= zmax1)
	{
	  double da;
	  char snum[30];
	  int i;
	  if (aa >= zmin1)
	    {
	      if (naxis == 0)
		{
		  m = scx (aa);
		  if (!texflag)
		    tg_line (m, Y1, m, Y1 + hash);
		}
	      else
		{
		  m = scy (aa);
		  if (!texflag)
		    tg_line (X1 - hash, m, X1, m);
		}
	      if (islog)
		doubleToStr (aa, aa, snum);
	      else
		doubleToStr (aa, step, snum);
	      len = MAX (len, tg_textwidth (snum));
	      if (naxis == 0)
		tg_outtextxy (m, Y1 + hash, snum);
	      else
		tg_outtextxy (X1 - hash, m, snum);
	    }
	  if (islog)
	    da = aa * step - aa;
	  else
	    da = step;
	  da = da / n_sub;
	  for (i = 1; i < n_sub; i++)
	    {
	      aa += da;
	      if (!texflag && aa >= zmin1 && aa <= zmax1)
		{
		  if (naxis == 0)
		    {
		      m = scx (aa);
		      tg_line (m, Y1, m, Y1 + 2 * hash / 3);
		    }
		  else
		    {
		      m = scy (aa);
		      tg_line (X1 - 2 * hash / 3, m, X1, m);
		    }
		}
	    }
	  aa += da;
	}

      if (naxis == 0)
	{
	  tg_settextjustify (RightText, TopText);
	  tg_outtextxy (X2, Y1 + hash + th, xstr);
	}
      else if (texflag)
	{
	  f_printf (out_tex, "\\rText(%.1lf,%.1lf)[tr][l]{%s}\n",
		    texX (X1 - len - 3 * hash / 2), texY (Y2), ystr);

	}
      else
	{
	  tg_settextjustify (LeftText, CenterText);
	  tg_outtextxy (X1, 2 * th, ystr);
	}
    }
  if (texflag)
    f_printf (out_tex, "%% ============== end of axis ============\n");

}


static void
plot_curve (double xMin, double xMax, int dim, double *f)
{
  double x, y, xx, yy;
  int i;
  double step = (xMax - xMin) / (dim - 1);
  double ymax = dscy (grafminmax.ymin);
  double ymin = dscy (grafminmax.ymax);

  {
    for (i = 1; i < dim; i++)
      {
	x = dscx (xMin + (i - 1) * step);
	xx = dscx (xMin + i * step);
	y = dscy (f[i - 1]);
	yy = dscy (f[i]);

	if (yy < y)
	  {
	    double z;
	    z = yy;
	    yy = y;
	    y = z;
	    z = xx;
	    xx = x;
	    x = z;
	  }

	if (yy > ymin && y < ymax)
	  {
	    if (yy > ymax)
	      {
		xx = xx - ((xx - x) * (yy - ymax)) / (yy - y);
		yy = ymax;
	      }
	    if (y < ymin)
	      {
		x = x - ((x - xx) * (y - ymin)) / (y - yy);
		y = ymin;
	      }
	    tg_line ((int) x, (int) y, (int) xx, (int) yy);
	  }
      }
  }
}

static void
plot_spline (double xMin, double xMax, int dim, double *f,
	     double *b, double *c, double *d, int indicator)
{
  int i;
  if (texflag)
    {
      double step = (xMax - xMin) / (dim - 1);
      double ymax = dscy (grafminmax.ymin);
      double ymin = dscy (grafminmax.ymax);

      int npoint = 0;
      int stat = 0;
      double dx, dy;
      f_printf (out_tex, "\n");
      for (i = 0; i < dim; i++)
	{
	  dx = dscx (xMin + i * step);
	  dy = dscy (f[i]);
	  if (dy <= ymax && dy >= ymin)
	    {
	      if (stat == 0)
		{
		  f_printf (out_tex, "\\Curve{");
		  stat = 1;
		  npoint = 1;
		}
	      f_printf (out_tex, "(%.2f,%.2f)", texX (dx), texY (dy));
	      npoint++;
	      if (npoint == 5)
		{
		  f_printf (out_tex, "\n");
		  npoint = 0;
		}
	    }
	  else if (stat == 1)
	    {
	      stat = 0;
	      f_printf (out_tex, "}\n");
	    }
	}
      if (stat == 1)
	f_printf (out_tex, "}\n");
    }
  else
    {
      double step0 = (xMax - xMin) * 2. / (dim + 1);
      double step = 2. / xscale;
      double x1 = xMin + 0.5 * step0;
      double y1 = f[0];
      double x2 = x1 + step;
      double y2;


      if (indicator == 1)
	{
	  progonca (dim - 1, f, b, c, d);
	  step0 = (xMax - xMin) / dim;
	}

      for (; x2 <= (xMax - 0.5 * step0); x1 += step, x2 += step, y1 = y2)
	{
	  y2 =
	    spline_for_graph ((dim - 1) * (x2 - xMin - 0.5 * step0) / (xMax -
								       xMin -
								       step0),
			      f, b, c, d);
	  if (y1 < grafminmax.ymax && y2 < grafminmax.ymax
	      && y1 > grafminmax.ymin && y2 > grafminmax.ymin)
	    tg_line (scx (x1), scy (y1), scx (x2), scy (y2));
	}
    }
}

static void
plot_hist (double xMin, double xMax, int dim, double *f, double *df)
{
  double x, y, yy;
  int i;

  double ymax = dscy (grafminmax.ymin);
  double ymin = dscy (grafminmax.ymax);
  double step = (xMax - xMin) / dim;

  for (i = 0; i < dim; i++)
    {
      y = dscy (f[i]);
      if (y < ymax && y > ymin)
	tg_line ((int) dscx (xMin + i * step), (int) y,
		 (int) dscx (xMin + (i + 1) * step), (int) y);

      y = MIN (dscy (f[i] - df[i]), ymax);
      yy = MAX (dscy (f[i] + df[i]), ymin);
      x = (dscx (xMin + i * step) + dscx (xMin + (i + 1) * step)) / 2;
      if (y > yy)
	tg_line ((int) x, (int) y, (int) x, (int) yy);
    }
}



static void
writetable1 (double xMin, double xMax, int dim, double *f,
	     double *df, char *upstr, char *x_str, char *y_str)
{
  char filename[STRSIZ], buff[STRSIZ];
  FILE *outfile;
  int i;
  double dx;

  nextFileName (filename, "tab_", ".txt");
  strcat (filename, ".txt");

  if (dim > 1)
    {
      outfile = fopen (filename, "w");
      f_printf (outfile, " %s\n", upstr);
      if (df)
	f_printf (outfile, "x-axis: \"%s\"  from %f to %f N_bins= %d\n",
		  x_str, xMin, xMax, dim);
      else
	f_printf (outfile, "x-axis: \"%s\"  from %f to %f N_points= %d\n",
		  x_str, xMin, xMax, dim);
      f_printf (outfile, "%s", y_str);

      dx = (xMax - xMin) / (dim - 1);

      for (i = 0; i < dim; i++)
	{
	  f_printf (outfile, "\n%-12E %-12E", xMin + dx * i, f[i]);
	  if (df)
	    f_printf (outfile, " +/-  %-12E", df[i]);
	}
      f_printf (outfile, "\n");
      fclose (outfile);

      sprintf (buff, " You can find results in the file\n%s", filename);
    }
  else
    strcpy (buff, " Set the number of bins larger then 1\n");

  messanykey (10, 12, buff);
}


static void
writehisto1 (double xMin, double xMax, int dim, double *f,
	     double *df, char *upstr, char *x_str, char *y_str)
{
  char filename[STRSIZ], buff[STRSIZ];
  FILE *outfile;
  int i;
  double dx;

  nextFileName (filename, "tab_", ".txt");
  strcat (filename, ".txt");

  outfile = fopen (filename, "w");
  f_printf (outfile, " %s\n", upstr);
  if (df)
    {
      f_printf (outfile, "x-axis: \"%s\"  from %f to %f N_bins= %d\n", x_str,
		xMin, xMax, dim);
    }
  else
    {
      f_printf (outfile, "x-axis: \"%s\"  from %f to %f N_points= %d\n",
		x_str, xMin, xMax, dim);
    }
  f_printf (outfile, "%s", y_str);

  dx = (xMax - xMin) / dim;

  for (i = 0; i < dim; i++)
    {
      f_printf (outfile, "\n%-12E %-12E", xMin + dx * i + dx / 2., f[i]);
      if (df)
	f_printf (outfile, " +/-  %-12E", df[i]);
    }
  f_printf (outfile, "\n");

  fclose (outfile);
  sprintf (buff, " You can find results in the file\n%s", filename);
  messanykey (10, 12, buff);
}

static void 
set_axis_titles (double dX, char * x_str, char * xaxistitle, char * yaxistitle)
{
  sprintf (xaxistitle, "\"%s\"", x_str);
  sprintf (yaxistitle, "\"#sigma, fb\"");
  if (!strncmp (x_str, "Transverse momentum Pt%d", 21))
    {
      sprintf (xaxistitle, "\"P_{T}, GeV\"");
      sprintf (yaxistitle, "\"d#sigma/d P_{T}, fb/%.0f GeV\"", dX);
    }

  if (!strncmp (x_str, "Mass", 4))
    {
      sprintf (xaxistitle, "\"M, GeV\"");
      sprintf (yaxistitle, "\"d#sigma/d M, fb/%.0f GeV\"", dX);
    }

  if (!strncmp (x_str, "Energy E", 8))
    {
      sprintf (xaxistitle, "\"E, GeV\"");
      sprintf (yaxistitle, "\"d#sigma/d E, fb/%.0f GeV\"", dX);
    }

  if (!strncmp (x_str, "Angle", 5))
    {
      sprintf (xaxistitle, "\"#theta, deg\"");
      sprintf (yaxistitle, "\"d#sigma/d #theta, fb/%.0f deg\"", dX);
    }

  if (!strncmp (x_str, "pseudo-rapidity", 15))
    {
      sprintf (xaxistitle, "\"#eta\"");
      sprintf (yaxistitle, "\"d#sigma/d #eta, fb\"");
    }

  if (!strncmp (x_str, "Rapidity", 8))
    {
      sprintf (xaxistitle, "\"#y\"");
      sprintf (yaxistitle, "\"d#sigma/d #y, fb\"");
    }

  if (!strncmp (x_str, "Cosine", 6))
    {
      sprintf (xaxistitle, "\"cos(#theta)\"");
      sprintf (yaxistitle, "\"d#sigma/d cos(#theta), fb\"");
    }
}


static void
writeroothistoAll (double xMin, double xMax, int dim, double *f,
	        double *df, char *upstr, char *x_str, char *y_str)
{
  int i;
  char xaxistitle[128];
  char yaxistitle[128];
  double yMin = 0.0;
  double yMax = 0.0;
  double dX = (xMax - xMin) / dim;

  midstr outname;
  midstr outname_C;
  FILE *rf;
  nextFileName (outname, "root_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");
  rf = fopen (outname_C, "w");
  if (NULL == rf)
    messanykey (10, 12, scat ("Error! file\n%s can not be open", outname_C));

  if (3 > dim || 200 < dim)
    messanykey (10, 12, scat (" 2 < nPoints < 201, Now Npoints = %i", dim));

  set_axis_titles (dX, x_str, xaxistitle, yaxistitle);
  for (i = 0; i < dim; i++)
    {
      if (f[i] > yMax) yMax = f[i];
      if (f[i] < yMin) yMin = f[i];
    }

  f_printf (rf, "//%s\n//%s\n", upstr, x_str);
  f_printf (rf, "//X from %.2f to %.2f N_bins= %d\n", xMin, xMax, dim);
  f_printf (rf, "//%s\n", y_str);
  f_printf (rf, "//Y from %-12E to %-12E\n", yMin * 1.0e+3, yMax * 1.0e+3);
  f_printf (rf, "//combine= 0\n{\ngROOT->Reset();\n\n");
  f_printf (rf, "\ndouble Yarray[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      f_printf (rf, "%-12E", f[i] * 1.0e+3);
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");

  f_printf (rf, "\ndouble Yerr[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      if (df)
	f_printf (rf, "%-12E", df[i] * 1.0e+3);
      else
	f_printf (rf, "0.0");
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");
  f_printf (rf, "double Xarray[%d];\n", dim);
  f_printf (rf, "double Xerr[%d];\n", dim);
  f_printf (rf, "for(int i=0;i<%d;i++) {\n", dim);
  f_printf (rf, "  Xarray[i]=%e+i*%e;\n  Xerr[i]=%e;\n}\n\n", xMin, dX, dX / 2);
  f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
  f_printf (rf, "c1->SetFillColor(0);\n");
  f_printf (rf, "c1->SetBorderMode(0);\n");
  f_printf (rf, "c1->SetFrameBorderMode(0);\n");
  f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
  f_printf (rf, "//c1->SetLogy(1);  // Set Log. Scale;\n\n\n");
  f_printf (rf, "TGraphErrors *gr1 = new TGraphErrors(%d,Xarray,Yarray,Xerr,Yerr);\n", dim);
  f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
  f_printf (rf, "gr1->SetLineWidth(2);\n");
  f_printf (rf, "gr1->SetLineColor(4);\n");
  f_printf (rf, "gr1->SetLineStyle(1);\n");
  f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
  f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);
  f_printf (rf, "gr1->Draw(\"AC\");\n");
  f_printf (rf, "c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "c1->Print(\"%s.png\");\n", outname);
  f_printf (rf, "}\n");
  fclose (rf);
/*  messanykey (10, 12, scat (" ROOT histogram kept in\n%s", outname_C)); */
}



static void
writeroothisto (double xMin, double xMax, int dim, double *f,
	        double *df, char *upstr, char *x_str, char *y_str)
{
  int i;
  char xaxistitle[128];
  char yaxistitle[128];
  double yMin = 0.0;
  double yMax = 0.0;
  double dX = (xMax - xMin) / dim;

  midstr outname;
  midstr outname_C;
  FILE *rf;
  nextFileName (outname, "root_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");
  rf = fopen (outname_C, "w");
  if (NULL == rf)
    messanykey (10, 12, scat ("Error! file\n%s can not be open", outname_C));

  if (3 > dim || 200 < dim)
    messanykey (10, 12, scat (" 2 < nPoints < 201, Now Npoints = %i", dim));

  set_axis_titles (dX, x_str, xaxistitle, yaxistitle);
  for (i = 0; i < dim; i++)
    {
      if (f[i] > yMax) yMax = f[i];
      if (f[i] < yMin) yMin = f[i];
    }

  f_printf (rf, "//%s\n//%s\n", upstr, x_str);
  f_printf (rf, "//X from %.2f to %.2f N_bins= %d\n", xMin, xMax, dim);
  f_printf (rf, "//%s\n", y_str);
  f_printf (rf, "//Y from %-12E to %-12E\n", yMin * 1.0e+3, yMax * 1.0e+3);
  f_printf (rf, "//combine= 0\n{\ngROOT->Reset();\n\n");
  f_printf (rf, "\ndouble Yarray[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      f_printf (rf, "%-12E", f[i] * 1.0e+3);
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");

  f_printf (rf, "\ndouble Yerr[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      if (df)
	f_printf (rf, "%-12E", df[i] * 1.0e+3);
      else
	f_printf (rf, "0.0");
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");
  f_printf (rf, "double Xarray[%d];\n", dim);
  f_printf (rf, "double Xerr[%d];\n", dim);
  f_printf (rf, "for(int i=0;i<%d;i++) {\n", dim);
  f_printf (rf, "  Xarray[i]=%e+i*%e;\n  Xerr[i]=%e;\n}\n\n", xMin, dX, dX / 2);
  f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
  f_printf (rf, "c1->SetFillColor(0);\n");
  f_printf (rf, "c1->SetBorderMode(0);\n");
  f_printf (rf, "c1->SetFrameBorderMode(0);\n");
  f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
  f_printf (rf, "//c1->SetLogy(1);  // Set Log. Scale;\n\n\n");
  f_printf (rf, "TGraphErrors *gr1 = new TGraphErrors(%d,Xarray,Yarray,Xerr,Yerr);\n", dim);
  f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
  f_printf (rf, "gr1->SetLineWidth(2);\n");
  f_printf (rf, "gr1->SetLineColor(4);\n");
  f_printf (rf, "gr1->SetLineStyle(1);\n");
  f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
  f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);
  f_printf (rf, "gr1->Draw(\"AC\");\n");
  f_printf (rf, "c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "c1->Print(\"%s.png\");\n", outname);
  f_printf (rf, "}\n");
  fclose (rf);
  messanykey (10, 12, scat (" ROOT histogram kept in\n%s", outname_C));
}



static void
writeroottable (double xMin, double xMax, int dim, double *f,
	        double *df, char *upstr, char *x_str, char *y_str)
{
  int i;
  char xaxistitle[128];
  char yaxistitle[128];
  double yMin = 0.0;
  double yMax = 0.0;
  double dX = (xMax - xMin) / (dim - 1);

  midstr outname;
  midstr outname_C;
  FILE *rf;
  nextFileName (outname, "root_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");
  rf = fopen (outname_C, "w");
  if (NULL == rf)
    messanykey (10, 12, scat ("Error! file\n%s can not be open", outname_C));

  if (3 > dim || 200 < dim)
    messanykey (10, 12, scat (" 2 < nPoints < 201, Now Npoints = %i", dim));

  set_axis_titles (dX, x_str, xaxistitle, yaxistitle);
  for (i = 0; i < dim; i++)
    {
      if (f[i] > yMax) yMax = f[i];
      if (f[i] < yMin) yMin = f[i];
    }

  f_printf (rf, "//%s\n//%s\n", upstr, x_str);
  f_printf (rf, "//X from %.2f to %.2f N_bins= %d\n", xMin, xMax, dim);
  f_printf (rf, "//%s\n", y_str);
  f_printf (rf, "//Y from %-12E to %-12E\n", yMin * 1.0e+3, yMax * 1.0e+3);
  f_printf (rf, "//combine= 0\n{\ngROOT->Reset();\n\n");
  f_printf (rf, "\ndouble Yarray[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      f_printf (rf, "%-12E", f[i] * 1.0e+3);
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");

  f_printf (rf, "\ndouble Yerr[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      if (df)
	f_printf (rf, "%-12E", df[i] * 1.0e+3);
      else
	f_printf (rf, "0.0");
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");
  f_printf (rf, "double Xarray[%d];\n", dim);
  f_printf (rf, "double Xerr[%d];\n", dim);
  f_printf (rf, "for(int i=0;i<%d;i++) {\n", dim);
  f_printf (rf, "  Xarray[i]=%e+i*%e;\n  Xerr[i]=0.0;\n}\n\n", xMin, dX);
  f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
  f_printf (rf, "c1->SetFillColor(0);\n");
  f_printf (rf, "c1->SetBorderMode(0);\n");
  f_printf (rf, "c1->SetFrameBorderMode(0);\n");
  f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
  f_printf (rf, "//c1->SetLogy(1);  // Set Log. Scale;\n\n\n");
  f_printf (rf, "TGraphErrors *gr1 = new TGraphErrors(%d,Xarray,Yarray,Xerr,Yerr);\n", dim);
  f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
  f_printf (rf, "gr1->SetLineWidth(2);\n");
  f_printf (rf, "gr1->SetLineColor(4);\n");
  f_printf (rf, "gr1->SetLineStyle(1);\n");
  f_printf (rf, "gr1->SetMarkerColor(2);\n");
  f_printf (rf, "gr1->SetMarkerStyle(2);\n");
  f_printf (rf, "gr1->SetMarkerSize(1.5);\n");
  f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
  f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);
  f_printf (rf, "gr1->Draw(\"AP\");\n");
  f_printf (rf, "c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "c1->Print(\"%s.png\");\n", outname);
  f_printf (rf, "}\n");
  fclose (rf);
  messanykey (10, 12, scat (" ROOT histogram kept in\n%s", outname_C));
}

static void
writeroottable1d (double xMin, double xMax, int dim, double *f,
	        double *df, char *upstr, char *x_str, char *y_str)
{
  int i;
  char xaxistitle[128];
  char yaxistitle[128];
  double yMin = 0.0;
  double yMax = 0.0;
  double dX = (xMax - xMin) / (dim - 1);

  midstr outname;
  midstr outname_C;
  FILE *rf;
  nextFileName (outname, "hist1d_", ".C");
  strcpy (outname_C, outname);
  strcat (outname_C, ".C");
  rf = fopen (outname_C, "w");
  if (NULL == rf)
    messanykey (10, 12, scat ("Error! file\n%s can not be open", outname_C));

  if (3 > dim || 200 < dim)
    messanykey (10, 12, scat (" 2 < nPoints < 201, Now Npoints = %i", dim));

  set_axis_titles (dX, x_str, xaxistitle, yaxistitle);
  for (i = 0; i < dim; i++)
    {
      if (f[i] > yMax) yMax = f[i];
      if (f[i] < yMin) yMin = f[i];
    }

  f_printf (rf, "//%s\n//%s\n", upstr, x_str);
  f_printf (rf, "//X from %.2f to %.2f N_bins= %d\n", xMin, xMax, dim);
  f_printf (rf, "//%s\n", y_str);
  f_printf (rf, "//Y from %-12E to %-12E\n", yMin * 1.0e+3, yMax * 1.0e+3);
  f_printf (rf, "//combine= 0\n{\ngROOT->Reset();\n\n");
  f_printf (rf, "\ndouble Yarray[%d]={\n", dim);
  for (i = 0; i < dim; i++)
    {
      f_printf (rf, "%-12E", f[i] * 1.0e+3);
      if (i < (dim - 1))
	f_printf (rf, ", ");

      if (0 == (i + 1) % 5)
	f_printf (rf, "\n");
    }
  f_printf (rf, "};\n\n");


  f_printf (rf, "double Xarray[%d];\n", dim);

  f_printf (rf, "for(int i=0;i<%d;i++) Xarray[i]=%e+i*%e;\n\n",dim, xMin, dX);
  f_printf (rf, "TCanvas *c1 = new TCanvas(\"c1\",\" \",200,100,700,500);\n");
  f_printf (rf, "c1->SetFillColor(0);\n");
  f_printf (rf, "c1->SetBorderMode(0);\n");
  f_printf (rf, "c1->SetFrameBorderMode(0);\n");
  f_printf (rf, "gStyle->SetOptTitle(kFALSE);\n\n");
  f_printf (rf, "//c1->SetLogy(1);  // Set Log. Scale;\n\n\n");
  f_printf (rf, "TGraph *gr1 = new TGraph(%d,Xarray,Yarray);\n", dim);
  f_printf (rf, "TGaxis::SetMaxDigits(3);\n");
  f_printf (rf, "gr1->SetLineWidth(2);\n");
  f_printf (rf, "gr1->SetLineColor(4);\n");
  f_printf (rf, "gr1->SetLineStyle(1);\n");
  f_printf (rf, "gr1->GetXaxis()->SetTitle(%s);\n", xaxistitle);
  f_printf (rf, "gr1->GetYaxis()->SetTitle(%s);\n", yaxistitle);
  f_printf (rf, "gr1->Draw(\"AC\");\n");
  f_printf (rf, "c1->Print(\"%s.eps\");\n", outname);
  f_printf (rf, "c1->Print(\"%s.png\");\n", outname);
  f_printf (rf, "}\n");
  fclose (rf);
/*  messanykey (10, 12, scat (" ROOT histogram kept in\n%s", outname_C)); */
}


static void
writespline (double xMin, double xMax, int dim, double *f,
	     double *b, double *c, double *d, int SplinePoints, char *upstr,
	     char *x_str, char *y_str)
{
  char filename[STRSIZ], buff[STRSIZ];
  FILE *outfile;
  int i;
  double splinepoint, step;

  step = (2. * dim - 2.) / (SplinePoints - 1.);

  nextFileName (filename, "spline", ".txt");
  strcat (filename, ".txt");
  outfile = fopen (filename, "w");
  f_printf (outfile, " %s\n", upstr);
  f_printf (outfile, "x-axis: \"%s\"  from %f to %f N_points= %d\n", x_str,
	    xMin, xMax, SplinePoints);
  f_printf (outfile, "%s", y_str);


  for (i = 0; i < SplinePoints; i++)
    {
      splinepoint = spline_for_graph (step * i, f, b, c, d);
      f_printf (outfile, "\n%-12E", splinepoint);
    }
  f_printf (outfile, "\n");
  fclose (outfile);

  sprintf (buff, " You can find results in the file\n%s.eps", filename);
  messanykey (10, 12, buff);
}


static void
buildSpline (double xMin, double xMax, int dim, double *f, double *ff,
	     int *spline_on, int *SplinePoints, double *spline_f, double *b,
	     double *c, double *d)
{
  char menustr[STRSIZ];
  int changed = 0;
  int k = 1;
  void *pscr = NULL;
  static double chi2 = 1.;
  int SpPoints;

  SpPoints = (*SplinePoints);

  while (k)
    {
      strcpy (menustr, "\022"
	      " Spline       OFF " " chi^2/N= CHI     " " SplPoints= SP    ");
      improveStr (menustr, "CHI", "%-7.1G", chi2);
      improveStr (menustr, "SP", "%d", SpPoints);

      if (*spline_on)
	improveStr (menustr, "OFF", "ON");

      menu1 (60, 2, "", menustr, NULL, &pscr, &k);
      if (k)
	changed = 1;
      switch (k)
	{
	case 1:
	  *spline_on = !(*spline_on);
	  break;

	case 2:
	  correctDouble (33, 11, "chi^2/N = ", &chi2, 1);
	  if (chi2 > 10)
	    chi2 = 10;
	  if (chi2 < 0.1)
	    chi2 = 0.1;
	  break;

	case 3:
	  correctInt (33, 11, "SplinePoints = ", &SpPoints, 1);
	  if (SpPoints > 10000)
	    SpPoints = 10000;
	  if (SpPoints < 1)
	    SpPoints = 1;
	  break;
	}
    }

  if (spline_on)
    {
      SPLINE (chi2, dim - 1, f, ff, spline_f, b, c, d);
    }

  *SplinePoints = SpPoints;
}


void
plot_histoAll (double xMin, double xMax,
	    int dim, double *f, double *ff, char *upstr, char *xstr,
	    char *ystr)
{
  int k, nCol0, nRow0, key;
  char f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];
  double ymin, ymax;
  void *prtscr;

  double *spline_f = malloc ((2 * dim - 1) * sizeof (double));
  double *b = malloc ((2 * dim - 1) * sizeof (double));
  double *c = malloc ((2 * dim - 1) * sizeof (double));
  double *d = malloc ((2 * dim - 1) * sizeof (double));

  int spline_on = 0;
  int SplinePoints = (2 * dim - 1) * 2;

  get_text (1, 1, maxCol (), maxRow (), &prtscr);
  gminmax (f, ff, dim);

  grafminmax.xmax = xMax;
  grafminmax.xmin = xMin;

  if (1.e-4 * (fabs (grafminmax.xmax) + fabs (grafminmax.xmin))
      >= fabs (grafminmax.xmax - grafminmax.xmin))
    {
      warnanykey (10, 10, "Too short interval in X axis !");
      return;
    }

  ymin = grafminmax.ymin;
  ymax = grafminmax.ymax;

  logScale = (ymin > 0 && grafminmax.ymax / grafminmax.ymin > 10);
  k = 0;

  nCol0 = maxCol ();
  nRow0 = maxRow ();
  clr_scr (fgcolor, bkcolor);

  {
    char *histo_upstr = malloc ((13 + strlen (upstr)) * sizeof (char));
    sprintf (histo_upstr, "Histogram (%s)", upstr);
    gaxes (histo_upstr, xstr, ystr);
    free (histo_upstr);
  }

  if (ff)
    {
      plot_hist (xMin, xMax, dim, f, ff);
      if (spline_on)
	plot_spline (xMin, xMax, 2 * dim - 1, spline_f, b, c, d, 2);
    }
  else if (spline_on)
    plot_spline (xMin, xMax, dim, f, b, c, d, 1);
  else
    plot_curve (xMin, xMax, dim, f);


    writeroothistoAll (xMin, xMax, dim, f, ff, upstr, xstr, ystr);

  free (spline_f);
  free (b);
  free (c);
  free (d);

  clr_scr (FGmain, BGmain);
  put_text (&prtscr);
}




void
plot_histo (double xMin, double xMax,
	    int dim, double *f, double *ff, char *upstr, char *xstr,
	    char *ystr)
{
  int k, nCol0, nRow0, key;
  char f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];
  double ymin, ymax;
  void *prtscr;

  double *spline_f = malloc ((2 * dim - 1) * sizeof (double));
  double *b = malloc ((2 * dim - 1) * sizeof (double));
  double *c = malloc ((2 * dim - 1) * sizeof (double));
  double *d = malloc ((2 * dim - 1) * sizeof (double));

  int spline_on = 0;
  int SplinePoints = (2 * dim - 1) * 2;

  get_text (1, 1, maxCol (), maxRow (), &prtscr);
  gminmax (f, ff, dim);

  grafminmax.xmax = xMax;
  grafminmax.xmin = xMin;

  if (1.e-4 * (fabs (grafminmax.xmax) + fabs (grafminmax.xmin))
      >= fabs (grafminmax.xmax - grafminmax.xmin))
    {
      warnanykey (10, 10, "Too short interval in X axis !");
      return;
    }

  ymin = grafminmax.ymin;
  ymax = grafminmax.ymax;

  logScale = (ymin > 0 && grafminmax.ymax / grafminmax.ymin > 10);
  k = 0;
REDRAW:
  nCol0 = maxCol ();
  nRow0 = maxRow ();
  clr_scr (fgcolor, bkcolor);

  {
    char *histo_upstr = malloc ((13 + strlen (upstr)) * sizeof (char));
    sprintf (histo_upstr, "Histogram (%s)", upstr);
    gaxes (histo_upstr, xstr, ystr);
    free (histo_upstr);
  }

  if (ff)
    {
      plot_hist (xMin, xMax, dim, f, ff);
      if (spline_on)
	plot_spline (xMin, xMax, 2 * dim - 1, spline_f, b, c, d, 2);
    }
  else if (spline_on)
    plot_spline (xMin, xMax, dim, f, b, c, d, 1);
  else
    plot_curve (xMin, xMax, dim, f);

  if (texflag)
    {
      f_printf (out_tex, "\\end{picture}\n");
      texFinish ();
      sprintf (buff, "%s\n%s", "LaTeX output is saved in file", f_name);
      messanykey (35, 15, buff);
      goto contin;
    }
  tg_settextjustify (BottomText, LeftText);
  scrcolor (Red, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  do
    {
      key = inkey ();
    }
  while (key == KB_MOUSE);
  if (nCol0 != maxCol () && nRow0 != maxRow ())
    goto REDRAW;
  scrcolor (bkcolor, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  if (key == KB_ESC)
    goto exi;
contin:
  do
    {
      char sScale[20];
      void *pscr = NULL;

      if (logScale)
	strcpy (sScale, "Log.   ");
      else
	strcpy (sScale, "Lin.   ");

      sprintf (menustr, "%c Y-max = %-9.3G Y-min = %-9.3G Y-scale = %s"
	       " Write table      "
	       " Write ROOT-hist  "
	       " Spline       OFF "
	       " LaTeX            "
	       " Redraw plot      "
	       " EXIT PLOT        ", 18, grafminmax.ymax, grafminmax.ymin,
	       sScale);
      if (spline_on)
	improveStr (menustr, "OFF", "ON");

      menu1 (nCol0 - 20, 2, "", menustr, "n_plot_*", &pscr, &k);

      switch (k)
	{
	case 0:
	  goto REDRAW;
	case 1:
	  correctDouble (33, 11, "Y-max = ", &grafminmax.ymax, TRUE);
	  break;
	case 2:
	  correctDouble (33, 11, "Y-min = ", &grafminmax.ymin, TRUE);
	  break;
	case 3:
	  logScale = !logScale;
	  break;
	case 4:
	  writehisto1 (xMin, xMax, dim, f, ff, upstr, xstr, ystr);
	  if (ff && spline_on)
	    writespline (xMin, xMax, dim, spline_f, b, c, d, SplinePoints,
			 upstr, xstr, ystr);
	  break;
	case 5:
	  writeroothisto (xMin, xMax, dim, f, ff, upstr, xstr, ystr);
	  break;
	case 6:
	  if (ff)
	    buildSpline (xMin, xMax, dim, f, ff, &spline_on, &SplinePoints,
			 spline_f, b, c, d);
	  else
	    spline_on = !spline_on;
	  break;
	case 7:
	case 8:
	  if (grafminmax.ymin >= ymax || grafminmax.ymax <= ymin ||
	      grafminmax.ymin >= grafminmax.ymax)
	    {
	      warnanykey (10, 10, " Wrong Y-range");
	      break;
	    }
	  if (logScale && (grafminmax.ymin <= 0 ||
			   grafminmax.ymax / grafminmax.ymin <= 10))
	    {
	      messanykey (10, 10, " To use the logarithmic scale,\n"
			  " please, set Ymin and Ymax limits\n"
			  " such that Ymax > 10*Ymin > 0");
	      logScale = 0;
	    }
	  if (k == 7)
	    {
	      if (!texmenu (&pictureX, &pictureY, letterSize))
		break;
	      nextFileName (f_name, "plot_", ".tex");
	      strcat (f_name, ".tex");
	      texStart (f_name, upstr, letterSize);
	      texPicture (0, 0, tg_getmaxx (),
			  tg_getmaxy () - tg_textheight ("0"), pictureX,
			  pictureY);
	      f_printf (out_tex, "\\begin{picture}(%d,%d)(0,0)\n", pictureX,
			pictureY);
	    }
	  del_text (&pscr);
	  goto REDRAW;
	case 9:
	  put_text (&pscr);
	  goto exi;
	}
      if (nCol0 != maxCol () && nRow0 != maxRow ())
	goto REDRAW;
    }
  while (TRUE);

exi:
  free (spline_f);
  free (b);
  free (c);
  free (d);

  clr_scr (FGmain, BGmain);
  put_text (&prtscr);
}


void
plot_table (double xMin, double xMax,
	    int dim, double *f, double *ff, char *upstr, char *xstr,
	    char *ystr)
{
  int k, nCol0, nRow0, key;
  char f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];
  double ymin, ymax;
  void *prtscr;

  get_text (1, 1, maxCol (), maxRow (), &prtscr);
  gminmax (f, ff, dim);

  grafminmax.xmax = xMax;
  grafminmax.xmin = xMin;

  if (1.e-4 * (fabs (grafminmax.xmax) + fabs (grafminmax.xmin))
      >= fabs (grafminmax.xmax - grafminmax.xmin))
    {
      warnanykey (10, 10, "Too short interval in X axis !");
      return;
    }

  ymin = grafminmax.ymin;
  ymax = grafminmax.ymax;

  logScale = (ymin > 0 && grafminmax.ymax / grafminmax.ymin > 10);
  k = 0;
REDRAW:
  nCol0 = maxCol ();
  nRow0 = maxRow ();
  clr_scr (fgcolor, bkcolor);

  {
    char *histo_upstr = malloc ((13 + strlen (upstr)) * sizeof (char));
    sprintf (histo_upstr, "Table (%s)", upstr);
    gaxes (histo_upstr, xstr, ystr);
    free (histo_upstr);
  }

  if (ff)
    plot_hist (xMin, xMax, dim, f, ff);
  else
    plot_curve (xMin, xMax, dim, f);

  if (texflag)
    {
      f_printf (out_tex, "\\end{picture}\n");
      texFinish ();
      sprintf (buff, "LaTeX output is saved in file\n%s", f_name);
      messanykey (35, 15, buff);
      goto contin;
    }
  tg_settextjustify (BottomText, LeftText);
  scrcolor (Red, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  do
    key = inkey ();
  while (key == KB_MOUSE);
  if (nCol0 != maxCol () && nRow0 != maxRow ())
    goto REDRAW;
  scrcolor (bkcolor, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  if (key == KB_ESC)
    goto exi;
contin:
  do
    {
      char sScale[20];
      void *pscr = NULL;

      if (logScale)
	strcpy (sScale, "Log.   ");
      else
	strcpy (sScale, "Lin.   ");

      sprintf (menustr, "%c Y-max = %-9.3G Y-min = %-9.3G Y-scale = %s"
	       " Write table      "
	       " Write ROOT table "
	       " LaTeX            "
	       " Redraw plot      "
	       " EXIT PLOT        ",
	       18, grafminmax.ymax, grafminmax.ymin, sScale);

      menu1 (nCol0 - 20, 2, "", menustr, "n_plot_*", &pscr, &k);

      switch (k)
	{
	case 0:
	  goto REDRAW;
	case 1:
	  correctDouble (33, 11, "Y-max = ", &grafminmax.ymax, TRUE);
	  break;
	case 2:
	  correctDouble (33, 11, "Y-min = ", &grafminmax.ymin, TRUE);
	  break;
	case 3:
	  logScale = !logScale;
	  break;
	case 4:
	  writetable1 (xMin, xMax, dim, f, ff, upstr, xstr, ystr);
	  break;
	case 5:
	  writeroottable (xMin, xMax, dim, f, ff, upstr, xstr, ystr);
	  break;
	case 6:
	case 7:
	  if (grafminmax.ymin >= ymax || grafminmax.ymax <= ymin ||
	      grafminmax.ymin >= grafminmax.ymax)
	    {
	      warnanykey (10, 10, " Wrong Y-range");
	      break;
	    }
	  if (logScale && (grafminmax.ymin <= 0 ||
			   grafminmax.ymax / grafminmax.ymin <= 10))
	    {
	      messanykey (10, 10, " To use the logarithmic scale,\n"
			  " please, set Ymin and Ymax limits\n"
			  " such that Ymax > 10*Ymin > 0");
	      logScale = 0;
	    }
	  if (k == 6)
	    {
	      if (!texmenu (&pictureX, &pictureY, letterSize))
		break;
	      nextFileName (f_name, "plot_", ".tex");
	      strcat (f_name, ".tex");
	      texStart (f_name, upstr, letterSize);
	      texPicture (0, 0, tg_getmaxx (),
			  tg_getmaxy () - tg_textheight ("0"), pictureX,
			  pictureY);
	      f_printf (out_tex, "\\begin{picture}(%d,%d)(0,0)\n", pictureX,
			pictureY);
	    }
	  del_text (&pscr);
	  goto REDRAW;
	case 9:
	  put_text (&pscr);
	  goto exi;
	}
      if (nCol0 != maxCol () && nRow0 != maxRow ())
	goto REDRAW;
    }
  while (TRUE);

exi:
  clr_scr (FGmain, BGmain);
  put_text (&prtscr);
}



void
plot_table1d (double xMin, double xMax,
	    int dim, double *f, double *ff, char *upstr, char *xstr,
	    char *ystr)
{
  int k, nCol0, nRow0, key;
  char f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];
  double ymin, ymax;
  void *prtscr;

  get_text (1, 1, maxCol (), maxRow (), &prtscr);
  gminmax (f, ff, dim);

  grafminmax.xmax = xMax;
  grafminmax.xmin = xMin;

  if (1.e-4 * (fabs (grafminmax.xmax) + fabs (grafminmax.xmin))
      >= fabs (grafminmax.xmax - grafminmax.xmin))
    {
      warnanykey (10, 10, "Too short interval in X axis !");
      return;
    }

  ymin = grafminmax.ymin;
  ymax = grafminmax.ymax;

/*
  logScale = (ymin > 0 && grafminmax.ymax / grafminmax.ymin > 10);
  k = 0;
REDRAW:
  nCol0 = maxCol ();
  nRow0 = maxRow ();
  clr_scr (fgcolor, bkcolor);

  {
    char *histo_upstr = malloc ((13 + strlen (upstr)) * sizeof (char));
    sprintf (histo_upstr, "Table (%s)", upstr);
    gaxes (histo_upstr, xstr, ystr);
    free (histo_upstr);
  }

  if (ff)
    plot_hist (xMin, xMax, dim, f, ff);
  else
    plot_curve (xMin, xMax, dim, f);

  if (texflag)
    {
      f_printf (out_tex, "\\end{picture}\n");
      texFinish ();
      sprintf (buff, "LaTeX output is saved in file\n%s", f_name);
      messanykey (35, 15, buff);
      goto contin;
    }
  tg_settextjustify (BottomText, LeftText);
  scrcolor (Red, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  do
    key = inkey ();
  while (key == KB_MOUSE);
  if (nCol0 != maxCol () && nRow0 != maxRow ())
    goto REDRAW;
  scrcolor (bkcolor, bkcolor);
  tg_outtextxy (0, tg_getmaxy (),
		"Press the Esc key to exit plot or other key to get the menu");
  if (key == KB_ESC)
    goto exi;
contin:
*/

	  writeroottable1d (xMin, xMax, dim, f, ff, upstr, xstr, ystr);


exi:
  clr_scr (FGmain, BGmain);
  put_text (&prtscr);
}
