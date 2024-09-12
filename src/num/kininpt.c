/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Dmitry Kovalenko 
* ------------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "chep_crt/include/crt_util.h"
#include "out_ext.h"

#include "kinaux.h"
#include "kininpt.h"

static void wrt_ (char *mes, char *lvv) {
  int i;

  scrcolor (Blue, BGmain);
  print (mes);
  scrcolor (FGmain, BGmain);
  for (i = 0; lvv[i]; i++) {
    print ("%d", lvv[i]);
  }
}

static int fill_ (char * lvar) {
  int j;
  int nn;
  int ncr;
  midstr rdar;

  rdar[0] = 0;
  if (str_redact (rdar, 1, 10) == KB_ESC) {
    return 1;
  }

  lvar[0] = 0;
  ncr = 0;

  for (j = 0; j < strlen (rdar); ++j) {
    nn = rdar[j] - '0';
    if (nn > nin_ && nn <= nin_ + nout_ && !strchr (lvar, nn)) {
      lvar[ncr] = nn;
      lvar[++ncr] = 0;
    }
  }

  if (0 == ncr) {
    warnanykey (10, 10, "ERROR: 1st cluster is empty!");
    return -1;
  }

  return 0;
}

static int CorrectKinScheme (void) {
  int i, l1, l2, nc, i1;
  vshortstr icmp;
  int i_moth[10];
  int k_moth[10];
  kinmtc_ kinmtc_loc[10];

  for (i = 0; i < 10; ++i) {
    strcpy (kinmtc_loc[i].lvin, kinmtc_1[i].lvin);
    strcpy (kinmtc_loc[i].lvout[0], kinmtc_1[i].lvout[0]);
    strcpy (kinmtc_loc[i].lvout[1], kinmtc_1[i].lvout[1]);
  }
  i = 0;
L1:

/* * Fill the vector of the in-particles for this decay */
  if (0 == i) {
    kinmtc_loc[0].lvin[0] = 1;
    if (2 == nin_) {
      kinmtc_loc[0].lvin[1] = 2;
    }
    kinmtc_loc[0].lvin[nin_] = 0;
  } else {
    strcpy (kinmtc_loc[i].lvin, kinmtc_loc[i_moth[i]].lvout[k_moth[i]]);
  }

L5:
  goto_xy (3, 10 + i);
  print ("%40.40s", "");
  goto_xy (3, 10 + i);
  wrt_ ("in=  ", kinmtc_loc[i].lvin);

  if (i > 0 && kinmtc_loc[i].lvin[2] == 0) {
    kinmtc_loc[i].lvout[0][0] = kinmtc_loc[i].lvin[0];
    kinmtc_loc[i].lvout[1][0] = kinmtc_loc[i].lvin[1];
    kinmtc_loc[i].lvout[0][1] = 0;
    kinmtc_loc[i].lvout[1][1] = 0;
  } else {
    if (i == 0 && nout_ == 2) {
      kinmtc_loc[0].lvout[0][0] = nin_ + 1;
      kinmtc_loc[0].lvout[1][0] = nin_ + 2;
      kinmtc_loc[0].lvout[0][1] = 0;
      kinmtc_loc[0].lvout[1][1] = 0;
    } else {
      int j;
      int stat;
      goto_xy (13, 10 + i);
      scrcolor (Blue, BGmain);
      print ("-> out1= ");
      stat = fill_ (kinmtc_loc[i].lvout[0]);
      switch (stat) {
        case 1:
          return 0;
        case -1:
          if (i > 0) i--;
          goto_xy (1, where_y ());
          clr_eol ();
          goto L1;
        case 0: 
          break;
      }

      nc = 0;
      i1 = 0;
      strcpy (icmp, kinmtc_loc[i].lvin);
      if (i == 0) {
        lvmirr_ (icmp);
      }

      for (j = 0; j < strlen (icmp); j++) {
        if (strchr (kinmtc_loc[i].lvout[0], icmp[j])) {
          nc++;
        } else {
          kinmtc_loc[i].lvout[1][i1++] = icmp[j];
        }
      }

      kinmtc_loc[i].lvout[1][i1] = 0;

      if (nc != strlen (kinmtc_loc[i].lvout[0])) {
        warnanykey (10, 10, "ERROR: particle(s) have to be from out state");
        goto L5;
      }

      if (!strlen (kinmtc_loc[i].lvout[1])) {
        warnanykey (10, 20, "ERROR: 2-nd cluster is empty");
        goto L5;
      }
    }
  }
  goto_xy (13, 10 + i);
  wrt_ ("-> out1= ", kinmtc_loc[i].lvout[0]);
  goto_xy (25, 10 + i);
  wrt_ ("out2= ", kinmtc_loc[i].lvout[1]);

  l1 = strlen (kinmtc_loc[i].lvout[0]);
  l2 = strlen (kinmtc_loc[i].lvout[1]);

  if (l1 > 1) {
    i_moth[i + 1] = i;
    k_moth[i + 1] = 0;
  }

  if (l2 > 1) {
    i_moth[i + l1] = i;
    k_moth[i + l1] = 1;
  }

  if (i < nout_ - 2) {
    ++i;
    goto L1;
  }

  for (i = 0; i < 10; ++i) {
    strcpy (kinmtc_1[i].lvin, kinmtc_loc[i].lvin);
    strcpy (kinmtc_1[i].lvout[0], kinmtc_loc[i].lvout[0]);
    strcpy (kinmtc_1[i].lvout[1], kinmtc_loc[i].lvout[1]);
  }
  return 1;
}


int EnterKinScheme (void) {
  int i;
  int status;

  clrbox (1, 5, 53, 100);
  goto_xy (3, 5);
  scrcolor (Red, BGmain);
  print ("========= Current kinematical scheme =========\n");
  scrcolor (FGmain, BGmain);
  for (i = 0; i < nout_ - 1; ++i) {
    goto_xy (3, 6 + i);
    print ("%40.40s", "");
    goto_xy (3, 6 + i);
    wrt_ ("in= ", kinmtc_1[i].lvin);
    goto_xy (13, 6 + i);
    wrt_ ("-> out1= ", kinmtc_1[i].lvout[0]);
    goto_xy (25, 6 + i);
    wrt_ ("out2= ", kinmtc_1[i].lvout[1]);
  }
  goto_xy (3, 5 + nout_ - 1 + 1);
  scrcolor (Red, BGmain);
  print ("==============================================\n");

  status = CorrectKinScheme ();
  clrbox (1, 5, 60, 20);
  return status;
}

void InitKinScheme (void) {
  int i, j;
  int i2;

  for (i = 0; i < nout_ - 1; ++i) {
    kinmtc_1[i].lvout[0][0] = i + nin_ + 1;
    kinmtc_1[i].lvout[0][1] = 0;

    i2 = nout_ - i - 1;
    for (j = 0; j < i2; ++j) {
      kinmtc_1[i].lvout[1][j] = i + nin_ + j + 2;
    }
    kinmtc_1[i].lvout[1][i2] = 0;
  }

  kinmtc_1[0].lvin[0] = 1;
  if (2 == nin_) {
    kinmtc_1[0].lvin[1] = 2;
  }
  kinmtc_1[0].lvin[nin_] = 0;

  for (i = 1; i < nout_ - 1; ++i) {
    strcpy (kinmtc_1[i].lvin, kinmtc_1[i - 1].lvout[1]);
  }
}


int WriteKinScheme (FILE * nchan) {
  int i;

  fprintf (nchan, "\n");
  for (i = 0; i < nout_ - 1; ++i) {
    int k, c;
    for (k = 0; (c = kinmtc_1[i].lvin[k]); k++) {
      fprintf (nchan, "%d", c);
    }
    fprintf (nchan, " -> ");
    for (k = 0; (c = kinmtc_1[i].lvout[0][k]); ++k) {
      fprintf (nchan, "%d", c);
    }
    fprintf (nchan, " , ");
    for (k = 0; (c = kinmtc_1[i].lvout[1][k]); ++k) {
      fprintf (nchan, "%d", c);
    }
    fprintf (nchan, "\n");
  }
  return 0;
}

int ReadKinScheme (FILE * nchan) {
  int i;
  vshortstr strin;
  vshortstr strout1;
  vshortstr strout2;

  for (i = 0; i < nout_ - 1; ++i) {
    int l, k, c;
    fscanf (nchan, "%s -> %s , %s", strin, strout1, strout2);

    for (k = 0, l = 0; (c = strin[k]); k++)
      if (c != ' ')
        kinmtc_1[i].lvin[l++] = c - '0';
    kinmtc_1[i].lvin[l] = 0;

    for (k = 0, l = 0; (c = strout1[k]); k++)
      if (c != ' ')
        kinmtc_1[i].lvout[0][l++] = c - '0';
    kinmtc_1[i].lvout[0][l] = 0;

    for (k = 0, l = 0; (c = strout2[k]); k++)
      if (c != ' ')
        kinmtc_1[i].lvout[1][l++] = c - '0';
    kinmtc_1[i].lvout[1][l] = 0;
  }

  return 0;
}
