/* 
* Copyright (C) 2008-2009, CompHEP Collaboration
* Copyright (C) 2008, Slava Bunichev
* ----------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define CASMAXINOUT 20

#include "service2/include/chep_limits.h"
#include "service2/include/4_vector.h"

#include "cascade_cpyth1.h"

static void
new_prod_color_fact (int res_pos, char *color)
{
  int i;
  int c1, c2;
  const char *ptr = color;
  char new_color[64];
  char temp_str[7];

  new_color[0] = 0;
  while (sscanf (ptr, "(%d %d)%n", &c1, &c2, &i) == 2) {
    if (c1 > res_pos + 2) c1 -= 1;
    if (c2 > res_pos + 2) c2 -= 1;

    sprintf (temp_str, "(%d %d)", c1, c2);
    strcat (new_color, temp_str);
    ptr += i;
  }

  strcpy (color, new_color);
}

static void
new_decay_color_fact (int n_prod, char *color)
{
  int i, c1, c2;
  char new_color[64];
  char temp_str[7];
  const char *ptr = color;

  new_color[0] = 0;
  while (sscanf (ptr, "(%d %d)%n", &c1, &c2, &i) == 2) {
    sprintf (temp_str, "(%d %d)", c1 + n_prod, c2 + n_prod);
    strcat (new_color, temp_str);
    ptr += i;
  }

  strcpy (color, new_color);
}

int
cascade_cpyth1 (int regime, const char prd_name[], const char dec_name[], const char out_name[])
{
  int i, j;
  int err;
  int pr_num = 0;
  int stop = 0;
  int n_prod, n_decay;
  int N_events;
  int N_tot_subproc;
  int tot_dec_events, tot_prod_events;
  int res_pos[200];
  int prod_i, deca_i;

  double tmp_dbl;
  double Width, Cross_section;
  double QCD, res_mass;
  double m1[CASMAXINOUT];
  double m2[CASMAXINOUT];
  double decay_mass[200][20];

  char c;
  char part3[5];
  char buff[1024];
  char prod[16][16], decay[16][16];
  char tmp_string[64];
  char CompHEPversion[64];

  FILE * pfile = fopen (prd_name, "r");
  FILE * dfile = fopen (dec_name, "r");
  FILE * ofile  = fopen (out_name, "w+");
  if (NULL == pfile || NULL == dfile || NULL == ofile) {
    return -1;
  }

  fprintf (ofile, "#PEVLIB_v.1.0 =============================================\n");

  err = fscanf (dfile, "#PEVLIB_v.1.0 =============================================\n");
  err = fscanf (dfile, "#CompHEP version %s\n", CompHEPversion);
  err = fscanf (dfile, "#PROCESS %s -> ", part3);
  n_decay = 0;
  do {
    ++n_decay;
    fscanf (dfile, "%s", decay[n_decay]);
  }
  while (getc (dfile) != '\n');
  n_decay--;
  err = fscanf (dfile, "  SQRT(S) %le\n", &tmp_dbl);
  err = fscanf (dfile, "  Rapidity(c.m.s) %le\n", &tmp_dbl);
  err = fscanf (dfile, "  StrFun1:%s\n", tmp_string);
  err = fscanf (dfile, "  StrFun2:%s\n", tmp_string);
  err = fscanf (dfile, "#MASSES  %le", &m2[0]);
  for (i = 0; i < n_decay; ++i) {
    err = fscanf (dfile, " %le", &m2[i + 1]);
  }
  fscanf (dfile, "\n#Cross_section(Width) %le\n", &Width);
  fscanf (dfile, "#Number_of_events      %d\n", &N_events);
  fscanf (dfile,
          " #----------------------------------------------------------\n");
  err = fscanf (dfile, "#Number_of_subprocesses = %d\n", &N_tot_subproc);
  err = fscanf (dfile, "#Total_cross_section_(pb) = %le\n", &tmp_dbl);
  err = fscanf (dfile, "#Events_mixed_and_randomized = %d\n", &tot_dec_events);
  err = fscanf (dfile, "#Nproc ================== Events ==========================\n");
  if (N_events != tot_dec_events) {
    fprintf (stderr, "cascade (warning): strange decay event file, it can be prepared unproperly.\n");
    fprintf (stderr, "                   N_events in subprocesess = %i, N_events in procesess = %i\n", 
    N_events,  tot_dec_events);
  }
  if (1 != N_tot_subproc) {
    fprintf (stderr, "cascade (error): mixing with muli-subprocess decays (N = %i) is still not implemented.\n", N_tot_subproc);
    return -1;
  }

  while ((NULL == fgets (buff, 1024, pfile)) || !stop) {
    stop = sscanf (buff, "#Number_of_subprocesses = %d", &N_tot_subproc);
  }
  if (feof (pfile)) {
    fprintf (stderr, "cascade (error): strange format, can't find #Number_of_subprocesses string. Exit\n");
    return -1;
  }
  if (N_tot_subproc < 1) {
    fprintf (stderr, "cascade (error): strange number of subprocesses. Nproc = %i\n", N_tot_subproc);
    return -1;
  }

  err = -1;
  rewind (pfile);
  err = fscanf (pfile, "#PEVLIB_v.1.0 =============================================\n");
  while (pr_num < 1024) {
    int real_decay = 0;
    long curpos;
    double S, Rapidity;
    char part1[5], part2[5];
    char strfun_string1[64], strfun_string2[64];

    ++pr_num;
    err = fscanf (pfile, "#CompHEP version %s\n", CompHEPversion);
    err = fscanf (pfile, "#PROCESS %s %s -> ", part1, part2);
    n_prod = 0;
    do {
      ++n_prod;
      err = fscanf (pfile, "%s", prod[n_prod]);
    } while (getc (pfile) != '\n');
    n_prod--;
    err = fscanf (pfile, "  SQRT(S) %le\n", &S);
    err = fscanf (pfile, "  Rapidity(c.m.s) %le\n", &Rapidity);
    err = fscanf (pfile, "  StrFun1:%s\n", strfun_string1);
    err = fscanf (pfile, "  StrFun2:%s\n", strfun_string2);
    err = fscanf (pfile, "#MASSES  %le %le", &m1[1], &m1[2]);
    for (i = 0; i < n_prod; ++i) {
      err = fscanf (pfile, " %le", &m1[i + 3]);
    }
    err = fscanf (pfile, "\n#Cross_section(Width) %le\n", &Cross_section);
    err = fscanf (pfile, "#Number_of_events      %d\n", &N_events);
    fscanf (pfile, " #----------------------------------------------------------\n");
    fprintf (ofile, "#CompHEP version %s\n", CompHEPversion);
    fprintf (ofile, "#PROCESS %-4.4s %-4.4s->", part1, part2);
    res_pos[pr_num] = -1;
    for (i = 1; i <= n_prod; ++i) {
      if (!strcmp (prod[i], part3) && !real_decay) {
        real_decay = 1;
        res_pos[pr_num] = i;
        res_mass = m1[i + 2];
      } else {
        fprintf (ofile, " %-4.4s", prod[i]);
      }
    }
    for (i = 0; real_decay && i < n_decay; ++i) {
      fprintf (ofile, " %-4.4s", decay[i + 1]);
    }
    fprintf (ofile, "\n");
    fprintf (ofile, "#Initial_state\n");
    fprintf (ofile, "  SQRT(S) %.6E\n", S);
    fprintf (ofile, "  Rapidity(c.m.s) %.6E\n", Rapidity);
    fprintf (ofile, "  StrFun1:%s\n", strfun_string1);
    fprintf (ofile, "  StrFun2:%s\n", strfun_string2);
    fprintf (ofile, "#MASSES  %.10E %.10E", m1[1], m1[2]);
    for (i = 1; i <= n_prod; ++i) {
      if (i != res_pos[pr_num]) {
        fprintf (ofile, " %.10E", m1[i + 2]);
      } else {
        if (fabs (m2[0] - m1[i + 2]) > 1e-10) {
          fprintf (stderr, "Warining: the particle mass in production (M1 = %e) and in decay (M2 = %e) is different (DM = %e)!\n",
          m2[0], m1[i + 2], fabs (m2[0] - m1[i + 2]));
        }
      }
    }

    for (i = 0; real_decay && i < n_decay; ++i) {
      fprintf (ofile, " %.10E", m2[i + 1]);
      decay_mass[pr_num][i + 1] = m2[i + 1];
    }

    fprintf (ofile, "\n#Cross_section(Width) %.6E\n", Cross_section);
    fprintf (ofile, "#Number_of_events       %d\n", N_events);
    fprintf (ofile, "#----------------------------------------------------------\n");
    curpos = ftell (pfile);
    err = fscanf (pfile, "%s", tmp_string);
    if (!strncmp (tmp_string, "#Number_of_subprocesses", 23)) {
      double tot_cs;
      fscanf (pfile, " = %d\n", &N_tot_subproc);
      fscanf (pfile, "#Total_cross_section_(pb) = %le\n", &tot_cs);
      fscanf (pfile, "#Events_mixed_and_randomized = %d\n", &tot_prod_events);
      fscanf (pfile, "#Nproc ================== Events ==========================\n");

      fprintf (ofile, "#Number_of_subprocesses = %d\n", N_tot_subproc);
      fprintf (ofile, "#Total_cross_section_(pb) = %.6E\n", tot_cs);
      fprintf (ofile, "#Events_mixed_and_randomized = %d\n", tot_prod_events);
      fprintf (ofile, "#Nproc ================== Events ==========================\n");
      break;
    } else {
      fseek (pfile, curpos, SEEK_SET);
    }
  }
  if (pr_num > 1023) {
    fprintf (stderr, "cascade (warning): too many subprocesses in production!\n");
    fprintf (stderr, "                  not all event can be processed...\n");
    fprintf (stderr, "                  change the limit in the cascade code.\n");
  }

  for (prod_i = 0, deca_i = 0; prod_i < tot_prod_events && deca_i < tot_dec_events;) {
    int n_subproc;
    int real_decay = 0;
    double newrf[4];
    double px1[CASMAXINOUT], py1[CASMAXINOUT], pz1[CASMAXINOUT];
    double px2[CASMAXINOUT], py2[CASMAXINOUT], pz2[CASMAXINOUT];
    char color1[64], color2[64];
    ++prod_i;
    fscanf (pfile, " %d     %le  %le  ", &n_subproc, &pz1[0], &pz1[1]);
    if (n_subproc > N_tot_subproc) {
      fprintf (stderr, "cascade (warning): strange subproc number (n = %i, n_tot = %i).\n", n_subproc, N_tot_subproc);
      fprintf (stderr, "                   skip the event...\n");
      continue;
    }
    for (i = 0; i < n_prod; ++i) {
      fscanf (pfile, "%le %le %le ", &px1[i + 2], &py1[i + 2], &pz1[i + 2]);
    }
    fscanf (pfile, "%le", &QCD);
    c = getc (pfile);
    if (getc (pfile) != '\n') {
      c = getc (pfile);
      fscanf (pfile, "%[^\n]", color1);
    } else {
      strcpy (color1, "");
    }
    fprintf (ofile, " %-4d % .10E % .10E ", n_subproc, pz1[0], pz1[1]);

    for (i = 0; i < n_prod; ++i) {
      double QCDdec;
      if (i + 1 == res_pos[n_subproc] && 0 == real_decay) {
        int nnn;
        ++deca_i;
        fscanf (dfile, " %d     ", &nnn);
        for (j = 0; j < n_decay; ++j) {
          fscanf (dfile, "%le %le %le ", &px2[j], &py2[j], &pz2[j ]);
        }
        fscanf (dfile, "%le", &QCDdec);
        c = getc (dfile);
        if (getc (dfile) != '\n') {
          c = getc (dfile);
          fscanf (dfile, "%[^\n]", color2);
        } else {
          strcpy (color2, "");
        }
        real_decay = 1;
        newrf[1] = -px1[i + 2];
        newrf[2] = -py1[i + 2];
        newrf[3] = -pz1[i + 2];
        newrf[0] = sqrt (res_mass * res_mass + newrf[1] * newrf[1] + newrf[2] * newrf[2] + newrf[3] * newrf[3]);
      } else {
        fprintf (ofile, "% .10E % .10E % .10E ", px1[i + 2], py1[i + 2], pz1[i + 2]);
      }
    }

    for (i = 0; real_decay && i < n_decay; ++i) {
      double a[4], b[4];
      a[1] = px2[i];
      a[2] = py2[i];
      a[3] = pz2[i];
      a[0] = sqrt (decay_mass[res_pos[n_subproc]][i + 1] * decay_mass[res_pos[n_subproc]][i + 1] + 
                       a[1] * a[1] + a[2] * a[2] + a[3] * a[3]);
      lorenc (a, newrf, b);
      fprintf (ofile, "% .10E % .10E % .10E ", b[1], b[2], b[3]);
    }

    if (real_decay) {
      new_prod_color_fact (res_pos[n_subproc], color1);
      new_decay_color_fact (n_prod, color2);
      fprintf (ofile, " %.3E   %s%s\n", QCD, color1, color2);
    } else {
      fprintf (ofile, " %.3E   %s\n", QCD, color1);
    }
  }
  fclose (pfile);
  fclose (dfile);
  fclose (ofile);

  return 0;
}
