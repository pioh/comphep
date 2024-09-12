/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "service2/include/chep_limits.h"

#include "tag_reader.h"
#include "tag_parser.h"
#include "tag_writer.h"
#include "routines.h"

long comphep_particle_name_base(char * p)
{
  long kf = 0;

  if (!strcmp(p, "Proton") || !strcmp(p, "proton") || !strcmp(p, "p+")) kf = 2212;
  if (!strcmp(p, "AntiProton") || !strcmp(p, "antiProton") || 
      !strcmp(p, "anti-proton") || !strcmp(p, "pbar-")) kf = -2212;

  /* charged leptons */
  if (!strcmp(p, "e1") || !strcmp(p, "e") || !strcmp(p, "electron") || !strcmp(p, "e-")) kf = 11;
  if (!strcmp(p, "E1") || !strcmp(p, "E") || !strcmp(p, "positron") || !strcmp(p, "e+")) kf = -11;
  if (!strcmp(p, "e2") || !strcmp(p, "m")) kf = 13;
  if (!strcmp(p, "E2") || !strcmp(p, "M")) kf = -13;
  if (!strcmp(p, "e3") || !strcmp(p, "l")) kf = 15;
  if (!strcmp(p, "E3") || !strcmp(p, "L")) kf = -15;
  if (!strcmp(p, "e4")) kf = 17;
  if (!strcmp(p, "E4")) kf = -17;

  /* neutrinos */
  if (!strcmp(p, "n1") || !strcmp(p, "ne")) kf = 12;
  if (!strcmp(p, "N1") || !strcmp(p, "Ne")) kf = -12;
  if (!strcmp(p, "n2") || !strcmp(p, "nm")) kf = 14;
  if (!strcmp(p, "N2") || !strcmp(p, "Nm")) kf = -14;
  if (!strcmp(p, "n3") || !strcmp(p, "nl")) kf = 16;
  if (!strcmp(p, "N3") || !strcmp(p, "Nl")) kf = -16;
  if (!strcmp(p, "n4")) kf = 18;
  if (!strcmp(p, "N4")) kf = -18;

  /* hash quarks */
  if (!strcmp(p, "d#")) kf = 1;
  if (!strcmp(p, "D#")) kf = -1;
  if (!strcmp(p, "u#")) kf = 2;
  if (!strcmp(p, "U#")) kf = -2;
  if (!strcmp(p, "q#")) kf = 2;
  if (!strcmp(p, "Q#")) kf = -2;

  /* quarks */
  if (!strcmp(p, "d")) kf = 1;
  if (!strcmp(p, "D")) kf = -1;
  if (!strcmp(p, "u")) kf = 2;
  if (!strcmp(p, "U")) kf = -2;
  if (!strcmp(p, "s")) kf = 3;
  if (!strcmp(p, "S")) kf = -3;
  if (!strcmp(p, "c")) kf = 4;
  if (!strcmp(p, "C")) kf = -4;
  if (!strcmp(p, "b")) kf = 5;
  if (!strcmp(p, "B")) kf = -5;
  if (!strcmp(p, "t")) kf = 6;
  if (!strcmp(p, "T")) kf = -6;
  if (!strcmp(p, "b4")) kf = 7;
  if (!strcmp(p, "B4")) kf = -7;
  if (!strcmp(p, "t4")) kf = 8;
  if (!strcmp(p, "T4")) kf = -8;

  /* gauge bosons */
  if (!strcmp(p, "G") || !strcmp(p, "g"))     kf = 21;
  if (!strcmp(p, "A") || !strcmp(p, "gamma")) kf = 22;
  if (!strcmp(p, "Z") || !strcmp(p, "Z0"))    kf = 23;
  if (!strcmp(p, "W+")) kf = 24;
  if (!strcmp(p, "W-")) kf = -24;

  /* higgses */
  if (!strcmp(p, "H") || !strcmp(p, "h") || !strcmp(p, "h0")) kf = 25;
  if (!strcmp(p, "H0")) kf = 35;
  if (!strcmp(p, "A0") || !strcmp(p, "H3")) kf = 36;
  if (!strcmp(p, "H+")) kf = 37;
  if (!strcmp(p, "H-")) kf = -37;

  /* BSM particles */
  if (!strcmp(p, "Gr")) kf = 39;           /* Graviton */
  if (!strcmp(p, "R0")) kf = 41;
  if (!strcmp(p, "LQ")) kf = 42;
  if (!strcmp(p, "~g")) kf = 1000021;      /* gluino */

  if (!strcmp(p, "~dL")) kf =  1000001;
  if (!strcmp(p, "~DL")) kf = -1000001;
  if (!strcmp(p, "~dR")) kf =  2000001;
  if (!strcmp(p, "~DR")) kf = -2000001;
  if (!strcmp(p, "~uL")) kf =  1000002;
  if (!strcmp(p, "~UL")) kf = -1000002;
  if (!strcmp(p, "~uR")) kf =  2000002;
  if (!strcmp(p, "~UR")) kf = -2000002;
  if (!strcmp(p, "~sL")) kf =  1000003;
  if (!strcmp(p, "~SL")) kf = -1000003;
  if (!strcmp(p, "~sR")) kf =  2000003;
  if (!strcmp(p, "~SR")) kf = -2000003;
  if (!strcmp(p, "~cL")) kf =  1000004;
  if (!strcmp(p, "~CL")) kf = -1000004;
  if (!strcmp(p, "~cR")) kf =  2000004; 
  if (!strcmp(p, "~CR")) kf = -2000004;
  if (!strcmp(p, "~b1")) kf =  1000005;
  if (!strcmp(p, "~B1")) kf = -1000005;
  if (!strcmp(p, "~b2")) kf =  2000005;
  if (!strcmp(p, "~B2")) kf = -2000005;
  if (!strcmp(p, "~t1")) kf =  1000006;
  if (!strcmp(p, "~T1")) kf = -1000006;
  if (!strcmp(p, "~t2")) kf =  2000006;
  if (!strcmp(p, "~T2")) kf = -2000006;
  if (!strcmp(p, "~eL")) kf =  1000011;
  if (!strcmp(p, "~EL")) kf = -1000011;
  if (!strcmp(p, "~eR")) kf =  2000011;
  if (!strcmp(p, "~ER")) kf = -2000011;
  if (!strcmp(p, "~ne")) kf =  1000012;
  if (!strcmp(p, "~Ne")) kf = -1000012;
  if (!strcmp(p, "~mL")) kf =  1000013;
  if (!strcmp(p, "~ML")) kf = -1000013;
  if (!strcmp(p, "~mR")) kf =  2000013;
  if (!strcmp(p, "~MR")) kf = -2000013;
  if (!strcmp(p, "~nm")) kf =  1000014;
  if (!strcmp(p, "~Nm")) kf = -1000014;
  if (!strcmp(p, "~l1")) kf =  1000015;
  if (!strcmp(p, "~L1")) kf = -1000015;
  if (!strcmp(p, "~l2")) kf =  2000015;
  if (!strcmp(p, "~L2")) kf = -2000015;
  if (!strcmp(p, "~nl")) kf =  1000016;
  if (!strcmp(p, "~Nl")) kf = -1000016;

  if (!strcmp(p, "~o1")) kf =  1000022;   // neutralino chi_10
  if (!strcmp(p, "~o2")) kf =  1000023;   // neutralino chi_20
  if (!strcmp(p, "~o3")) kf =  1000025;   // neutralino chi_30
  if (!strcmp(p, "~o4")) kf =  1000035;   // neutralino chi_40
  if (!strcmp(p, "~1+")) kf =  1000024;   // chargino chi_1+
  if (!strcmp(p, "~1-")) kf = -1000024;   // chargino chi_1-
  if (!strcmp(p, "~2+")) kf =  1000037;   // chargino chi_2+
  if (!strcmp(p, "~2-")) kf = -1000037;   // chargino chi_2-

/* other lebel for chargino (from extra-dimention theories) */
  if (!strcmp(p, "KG"))  kf = 1000037;   // chargino chi_2-
  if (!strcmp(p, "~G"))  kf = 1000039;   // gravitino
  if (!strcmp(p, "~sS")) kf = 1000040;   // sgoldstino S
  if (!strcmp(p, "~sP")) kf = 1000041;   // sgoldstino P

return kf;
}

int check_format (FILE * f) {
  int ver = 0;
  midstr buff;

  fgets (buff, 1000, f);
  if ( NULL != strstr(buff, "#CompHEP version") || NULL != strstr(buff, "#PEVLIB_v.1.0")) {
    ver = 1;
  } else {
    if ( NULL != strstr(buff, "##") ) {
      while (NULL != fgets (buff, 1000, f))
        if ( NULL != strstr(buff, "CompHEP") ) {
          ver = 2;
          break;
        }
    }
  }
  rewind (f);
  return ver;
}

int get_pdfsup (char * p) {
  int pdfsup;

  if ( !strcmp (p, "cteq6l1") ) pdfsup = 10042; // LHAPDF
  if ( !strcmp (p, "cteq6l" ) ) pdfsup = 10041; // LHAPDF
  if ( !strcmp (p, "cteq6m" ) ) pdfsup = 10000; // LHAPDF
  if ( !strcmp (p, "cteq6d" ) ) pdfsup = 10200; // ????
  if ( !strcmp (p, "cteq5m1") ) pdfsup = 53;
  if ( !strcmp (p, "cteq5l" ) ) pdfsup = 46;
  if ( !strcmp (p, "cteq4l" ) ) pdfsup = 32;
  if ( !strcmp (p, "cteq4m" ) ) pdfsup = 34;

  return pdfsup;
}
