/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <math.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "chep_crt/include/chep_crt.h"
#include "service2/include/unix_utils.h"
#include "service2/include/tptcmac.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"

#include "physics.h"
#include "r_code.h"
#include "draw_ampl.h"
#include "process.h"
#include "process_core.h"

#include "showgrph.h"

static char stc[80];
static long nn;           /* nn - first FilePos of subproces */
static long nm;           /* nm-1 - last  FilePos of subproces */
static int pictureX;
static int pictureY;
static int xn;
static int ynu;

static FILE *diagrFile;
static int diagSize;
static char *mark;
static int squared;

static void 
wrttext (int n, int x, int y)
{
  int ypos = y + ynu - 2;

  tg_settextjustify (LeftText, BottomText);
  switch (mark[n - 1])
    {
    case -1:
      scrcolor (Red, White);
      tg_outtextxy (x + 48, ypos, "DEL");
      break;
    case -2:
      scrcolor (Red, White);
      tg_outtextxy (x + 48, ypos, "Out of memory");
      break;
    case 3:
      scrcolor (Red, White);
      tg_outtextxy (x + 48, ypos, "koeff");
      break;
    case 1:
      scrcolor (LightBlue, White);
      tg_outtextxy (x + 48, ypos, "CALC");
      break;
    case 2:
      scrcolor (LightBlue, White);
      tg_outtextxy (x + 48, ypos, "ZERO ");
    }
  scrcolor (Black, White);
}


static void 
diag_picture (int cur_diag_num, int x, int y)
{
  char buff[sizeof (adiagram)];

  wrttext (cur_diag_num, x, y);
  cur_diag_num += nn - 1;
  fseek (diagrFile, cur_diag_num * diagSize, SEEK_SET);
  fread (buff, diagSize, 1, diagrFile);

  vampl amplDiagram;
  mkverts ((particleNumType *) buff, &amplDiagram);
  drawAmpitudeDiagram (&amplDiagram, x, y);
}

static void 
sq_diag_picture (int cur_diag_num, int x, int y)
{
  char buff[sizeof (csdiagram)];

  wrttext (cur_diag_num, x, y);
  cur_diag_num += nn - 1;
  fseek (diagrFile, cur_diag_num * diagSize, SEEK_SET);
  fread (buff, diagSize, 1, diagrFile);

  amplitudeFrameX (getnout ());
  amplitudeFrameY (getnout ());
  drawSquaredDiagram ((csdiagram *) buff, x, x + xn - 1, y);
}


static void 
texinvoke (char *pname)
{
  int i, dn_;
  char buff[MAX (sizeof (csdiagram), sizeof (adiagram))];
  char dd[30];
  char f_name[STRSIZ];

  if (squared)
    strcpy (dd, "csd_");
  else
    strcpy (dd, "fd_");
  nextFileName (f_name, dd, ".tex");
  strcat (f_name, ".tex");
  texStart (f_name, scat ("diagrams for process %s", pname), "scriptsize");

  texxscale = 0.75;
  texyscale = 0.75;

  setPictureScale (squared, &pictureX, &pictureY);

  texxshift = 0;
  texymax1 = pictureY;

  pictureX *= texxscale;
  pictureY *= texyscale;

  dn_ = 0;
  f_printf (out_tex, "\\hspace{-8mm}\n");
  for (i = nn; i < nm; ++i)
    {
      if (mark[i - nn] != -1)
        {
          fseek (diagrFile, i * diagSize, SEEK_SET);
          fread (buff, diagSize, 1, diagrFile);
          dn_++;
          if (i != nn) f_printf (out_tex, "{} \\qquad\\allowbreak\n");
          f_printf (out_tex, "%%  diagram # %d\n", (i - (int) nn) + 1);
          f_printf (out_tex, "\\begin{picture}(%d,%d)(0,0)\n", pictureX, pictureY);

          if (squared) {
            amplitudeFrameX (getnout ());
            amplitudeFrameY (getnout ());
            drawSquaredDiagram ((csdiagram *) buff, 0, xn - 1, 0);
          } else {
            vampl amplDiagram;
            mkverts ((particleNumType *) buff, &amplDiagram);
            drawAmpitudeDiagram (&amplDiagram, 0, 0);
          }
          f_printf (out_tex, "\\Text(%d,0)[b] {diagr.%d}\n", pictureX / 2, dn_);
          f_printf (out_tex, "\\end{picture} \\ \n");
          f_printf (out_tex, "\\vspace{5mm}\n");
        }
    }

  texFinish ();
  setPictureScale (squared, &pictureX, &pictureY);
  messanykey (20, 14, scat ("LaTeX output is written in file\n %s", f_name));
}


static int 
comList (int n, char key)
{
  int l;

  switch (key)
    {
    case 'L':                   /* LaTex output for all undel. diagrs. */
      texinvoke (stc);
      return 0 /*nothing to redraw */ ;
    case 'D':                   /* Del -   Delete all subprocess  */
      for (l = 0; l < nm - nn; l++)
        {
          if (mark[l] == 0 || mark[l] == -2)
            mark[l] = -1;
        }
      return 2 /* redrawScreen */ ;
    case 'R':                   /*  Ins -  Restore all deleted  */
      for (l = 0; l < nm - nn; l++)
        {
          if (mark[l] == -1)
            mark[l] = 0;
        }
      return 2 /* redrawScreen */ ;
    case 'O':                   /*  Space  Bar   del/ins proceess  */
      if (mark[n - 1] == 0 || mark[n - 1] == -2)
        mark[n - 1] = -1;
      else if (mark[n - 1] == -1)
        mark[n - 1] = 0;
      return 1;
    case 'K':                   /*  Space  Bar   del/ins proceess  */
    /*  printf("\n K - key pressed\n"); */
      if (mark[n - 1] == 0 || mark[n - 1] == -2)
        mark[n - 1] = 3;
      else if (mark[n - 1] == 3)
        mark[n - 1] = 0;
      return 1;
    case 'G':
      {
        FILE *txt;
        makeghostdiagr (nn + n - 1, "view.tmp");
        txt = fopen ("view.tmp", "r");
        showtext (1, 1, 80, 1, "Ghosts", txt);
        fclose (txt);
        unlink ("view.tmp");
        return 0;
      }
    }

  fprintf (stderr, "***** comList: End reach. No int val. ?!\n");
  return 0;
}


int 
showgraphs (char upravl)
{
  int i;
  int ndel, nrest, ncalc, ntot;
  int onlyview = (upravl < 0);
  int delMarkPos;
  char help[80];
  char comLine[80];
  FILE * menuf;
    adiagram varc;
    csdiagram ars;

  if (onlyview)
    upravl = -upravl;

  squared = (upravl > 1);

  if (squared) {
    setPictureScale (squared, &xn, &ynu);
  } else {
    xn = amplitudeFrameX (getnout ());
    ynu = amplitudeFrameY (getnout ());
  }

  if (squared) {
    menuf      = fopen (MENUQ_NAME, "r+b");
    diagrFile  = fopen (DIAGRQ_NAME, "r+b");
    diagSize   = sizeof (csdiagram);
    delMarkPos = (char *) &(ars.status) - (char *) &ars;
  } else {
    menuf      = fopen (MENUP_NAME, "r+b");
    diagrFile  = fopen (DIAGRP_NAME, "r+b");
    diagSize   = sizeof (adiagram);
    delMarkPos = (char *) &(varc.delMark) - (char *) &varc;
  }

  if (1 == upravl) {
    rd_menu (menuf, upravl, nsub, stc, &ndel, &ncalc, &nrest, &nn);
  } else {
    rd_menu (menuf, upravl, nsub, stc, &ndel, &ncalc, &nrest, &nn);
  }
  ntot = ndel + ncalc + nrest;
  nm = nn + ntot;

  if (squared)
    {
      strcpy (help, "s_diag_s");
      strcpy (comLine, "Delete,On/off,Restore,Latex,Ghosts");
    }
  else if (onlyview)
    {
      strcpy (help, "s_diag_v");
      strcpy (comLine, "Latex");
    }
  else
    {
      strcpy (help, "s_diag_e");
      strcpy (comLine, "Delete,On/off,Restore,Latex");
    }

  mark = malloc (ntot);
  for (i = 0; i < ntot; i++)
    mark[i] = 0;

  for (i = nn; i < nm; i++)
    {
      fseek (diagrFile, i * diagSize + delMarkPos, SEEK_SET);
      fread (mark + i - nn, 1, 1, diagrFile);
    }

  if (squared) {
    pictures (ntot, xn, ynu, sq_diag_picture, comLine, comList, help);
  } else {
    pictures (ntot, xn, ynu, diag_picture, comLine, comList, help);
  }

  ndel = 0;
  ncalc = 0;
  for (i = nn; i < nm; i++)
    {
      fseek (diagrFile, i * diagSize + delMarkPos, SEEK_SET);
  /*    if (mark[i - nn] == 3) mark[i - nn] = 0;  */
      fwrite (mark + i - nn, 1, 1, diagrFile);
      if (mark[i - nn] > 0 && mark[i - nn]!= 3)
        ncalc++;
      else if (mark[i - nn] == -1)
        ndel++;
    }

  nrest = ntot - ndel - ncalc;
  if (1 == upravl) {
    wrt_menu (menuf, upravl, nsub, stc, ndel, ncalc, nrest, nn);
  } else {
    wrt_menu (menuf, upravl, nsub, stc, ndel, ncalc, nrest, nn);
  }

  fclose (menuf);
  fclose (diagrFile);
  free (mark);
  return 0;
}
