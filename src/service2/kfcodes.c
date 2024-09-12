/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* Author: Alexander Sherstnev
* ----------------------------------------------
*/
#include <string.h>
#include <stdio.h>

#include "kfcodes.h"

int 
kfbeam (char *bname)
{
  int i;
  if (!strcmp (bname, "Proton") ||
      !strcmp (bname, "proton") ||
      !strcmp (bname, "p+"))
    return 2212;

  if (!strcmp (bname, "AntiProton") ||
      !strcmp (bname, "Antiproton") ||
      !strcmp (bname, "antiProton") ||
      !strcmp (bname, "antiproton") ||
      !strcmp (bname, "anti-proton") ||
      !strcmp (bname, "pbar-"))
    return -2212;

  if (!strcmp (bname, "electron") || !strcmp (bname, "e-"))
    return 11;
  if (!strcmp (bname, "positron") || !strcmp (bname, "e+"))
    return -11;
  if (!strcmp (bname, "gamma"))
    return 22;
  i = kfpart (bname);
  if (i)
    return i;
  fprintf (stderr, "Unknown for CompHEP-PYTHIA Interface beam particle name:%s\n", bname);
  return 0;
}


char * bname (int kf) {
  if (kf ==  2212) return "proton" ;
  if (kf == -2212) return "anti-proton" ;
  if (kf ==    11) return "electron";
  if (kf ==   -11) return "positron";
  if (kf ==    22) return "gamma";

  return "unkown";
}


int 
kfpart (char *pname)
{
  if (!strcmp (pname, "e1"))    return 11;
  if (!strcmp (pname, "E1"))    return -11;
  if (!strcmp (pname, "e"))    return 11;
  if (!strcmp (pname, "E"))    return -11;
  if (!strcmp (pname, "e2"))    return 13;
  if (!strcmp (pname, "E2"))    return -13;
  if (!strcmp (pname, "m"))    return 13;
  if (!strcmp (pname, "M"))    return -13;
  if (!strcmp (pname, "e3"))    return 15;
  if (!strcmp (pname, "E3"))    return -15;
  if (!strcmp (pname, "l"))    return 15;
  if (!strcmp (pname, "L"))    return -15;
  if (!strcmp (pname, "e4"))    return 17;
  if (!strcmp (pname, "E4"))    return -17;
  if (!strcmp (pname, "n1"))    return 12;
  if (!strcmp (pname, "N1"))    return -12;
  if (!strcmp (pname, "ne"))    return 12;
  if (!strcmp (pname, "Ne"))    return -12;
  if (!strcmp (pname, "n2"))    return 14;
  if (!strcmp (pname, "N2"))    return -14;
  if (!strcmp (pname, "nm"))    return 14;
  if (!strcmp (pname, "Nm"))    return -14;
  if (!strcmp (pname, "n3"))    return 16;
  if (!strcmp (pname, "N3"))    return -16;
  if (!strcmp (pname, "nl"))    return 16;
  if (!strcmp (pname, "Nl"))    return -16;
  if (!strcmp (pname, "n4"))    return 18;
  if (!strcmp (pname, "N4"))    return -18;
  if (!strcmp (pname, "d#"))    return 1;
  if (!strcmp (pname, "D#"))    return -1;
  if (!strcmp (pname, "u#"))    return 2;
  if (!strcmp (pname, "U#"))    return -2;
  if (!strcmp (pname, "q#"))    return 2;
  if (!strcmp (pname, "Q#"))    return -2;
  if (!strcmp (pname, "d"))    return 1;
  if (!strcmp (pname, "D"))    return -1;
  if (!strcmp (pname, "u"))    return 2;
  if (!strcmp (pname, "U"))    return -2;
  if (!strcmp (pname, "s"))    return 3;
  if (!strcmp (pname, "S"))    return -3;
  if (!strcmp (pname, "c"))    return 4;
  if (!strcmp (pname, "C"))    return -4;
  if (!strcmp (pname, "b"))    return 5;
  if (!strcmp (pname, "B"))    return -5;
  if (!strcmp (pname, "t"))    return 6;
  if (!strcmp (pname, "T"))    return -6;
  if (!strcmp (pname, "b4"))    return 7;
  if (!strcmp (pname, "B4"))    return -7;
  if (!strcmp (pname, "t4"))    return 8;
  if (!strcmp (pname, "T4"))    return -8;
  if (!strcmp (pname, "G"))    return 21;
  if (!strcmp (pname, "g"))    return 21;
  if (!strcmp (pname, "A"))    return 22;
  if (!strcmp (pname, "gamma"))    return 22;
  if (!strcmp (pname, "Z"))    return 23;
  if (!strcmp (pname, "Z0"))    return 23;
  if (!strcmp (pname, "W+"))    return 24;
  if (!strcmp (pname, "W-"))    return -24;
  if (!strcmp (pname, "H"))    return 25;
  if (!strcmp (pname, "h"))    return 25;
  if (!strcmp (pname, "h0"))    return 25;
  if (!strcmp (pname, "H0"))    return 35;
  if (!strcmp (pname, "A0"))    return 36;
  if (!strcmp (pname, "H3"))    return 36;
  if (!strcmp (pname, "H+"))    return 37;
  if (!strcmp (pname, "H-"))    return -37;
  if (!strcmp (pname, "Gr"))    return 39;                  // Graviton

  if (!strcmp (pname, "R0"))    return 41;
  if (!strcmp (pname, "LQ"))    return 42;

  if (!strcmp (pname, "~g"))    return 1000021;             // gluino

  if (!strcmp (pname, "t3"))    return 32;             // spin-3/2 top-like particle
  if (!strcmp (pname, "T3"))    return -32;             // spin-3/2 top-like anti-particle

  if (!strcmp (pname, "~dL"))    return 1000001;
  if (!strcmp (pname, "~DL"))    return -1000001;
  if (!strcmp (pname, "~dR"))    return 2000001;
  if (!strcmp (pname, "~DR"))    return -2000001;
  if (!strcmp (pname, "~uL"))    return 1000002;
  if (!strcmp (pname, "~UL"))    return -1000002;
  if (!strcmp (pname, "~uR"))    return 2000002;
  if (!strcmp (pname, "~UR"))    return -2000002;
  if (!strcmp (pname, "~sL"))    return 1000003;
  if (!strcmp (pname, "~SL"))    return -1000003;
  if (!strcmp (pname, "~sR"))    return 2000003;
  if (!strcmp (pname, "~SR"))    return -2000003;
  if (!strcmp (pname, "~cL"))    return 1000004;
  if (!strcmp (pname, "~CL"))    return -1000004;
  if (!strcmp (pname, "~cR"))    return 2000004;
  if (!strcmp (pname, "~CR"))    return -2000004;
  if (!strcmp (pname, "~b1"))    return 1000005;
  if (!strcmp (pname, "~B1"))    return -1000005;
  if (!strcmp (pname, "~b2"))    return 2000005;
  if (!strcmp (pname, "~B2"))    return -2000005;
  if (!strcmp (pname, "~t1"))    return 1000006;
  if (!strcmp (pname, "~T1"))    return -1000006;
  if (!strcmp (pname, "~t2"))    return 2000006;
  if (!strcmp (pname, "~T2"))    return -2000006;
  if (!strcmp (pname, "~eL"))    return 1000011;
  if (!strcmp (pname, "~EL"))    return -1000011;
  if (!strcmp (pname, "~eR"))    return 2000011;
  if (!strcmp (pname, "~ER"))    return -2000011;
  if (!strcmp (pname, "~ne"))    return 1000012;
  if (!strcmp (pname, "~Ne"))    return -1000012;
  if (!strcmp (pname, "~mL"))    return 1000013;
  if (!strcmp (pname, "~ML"))    return -1000013;
  if (!strcmp (pname, "~mR"))    return 2000013;
  if (!strcmp (pname, "~MR"))    return -2000013;
  if (!strcmp (pname, "~nm"))    return 1000014;
  if (!strcmp (pname, "~Nm"))    return -1000014;
  if (!strcmp (pname, "~l1"))    return 1000015;
  if (!strcmp (pname, "~L1"))    return -1000015;
  if (!strcmp (pname, "~l2"))    return 2000015;
  if (!strcmp (pname, "~L2"))    return -2000015;
  if (!strcmp (pname, "~nl"))    return 1000016;
  if (!strcmp (pname, "~Nl"))    return -1000016;

  if (!strcmp (pname, "~o1"))    return 1000022;             // neutralino chi_10
  if (!strcmp (pname, "~o2"))    return 1000023;             // neutralino chi_20
  if (!strcmp (pname, "~o3"))    return 1000025;             // neutralino chi_30
  if (!strcmp (pname, "~o4"))    return 1000035;             // neutralino chi_40
  if (!strcmp (pname, "~1+"))    return 1000024;             // chargino chi_1+
  if (!strcmp (pname, "~1-"))    return -1000024;            // chargino chi_1-
  if (!strcmp (pname, "~2+"))    return 1000037;             // chargino chi_2+
  if (!strcmp (pname, "~2-"))    return -1000037;            // chargino chi_2-
  if (!strcmp (pname, "~G"))    return 1000039;             // gravitino
  if (!strcmp (pname, "~sS"))    return 1000040;             // sgoldstino S
  if (!strcmp (pname, "~sP"))    return 1000041;             // sgoldstino P

  return 0;
}


char * 
kfname (int kf)
{
  if (kf ==  11) return "e" ;
  if (kf == -11) return "E" ;
  if (kf ==  12) return "ne";
  if (kf == -12) return "Ne";
  if (kf ==  13) return "m" ;
  if (kf == -13) return "M" ;
  if (kf ==  14) return "nm";
  if (kf == -14) return "Nm";
  if (kf ==  15) return "l" ;
  if (kf == -15) return "L" ;
  if (kf ==  16) return "nl";
  if (kf == -16) return "Nl";
  if (kf ==  17) return "e4";
  if (kf == -17) return "E4";
  if (kf ==  18) return "n4";
  if (kf == -18) return "N4";
  if (kf ==   1) return "d" ;
  if (kf ==  -1) return "D" ;
  if (kf ==   2) return "u" ;
  if (kf ==  -2) return "U" ;
  if (kf ==   3) return "s" ;
  if (kf ==  -3) return "S" ;
  if (kf ==   4) return "c" ;
  if (kf ==  -4) return "C" ;
  if (kf ==   5) return "b" ;
  if (kf ==  -5) return "B" ;
  if (kf ==   6) return "t" ;
  if (kf ==  -6) return "T" ;
  if (kf ==   7) return "b4";
  if (kf ==  -7) return "B4";
  if (kf ==   8) return "t4";
  if (kf ==  -8) return "T4";
  if (kf ==  21) return "G" ;
  if (kf ==  21) return "g" ;
  if (kf ==  22) return "A" ;
  if (kf ==  23) return "Z" ;
  if (kf ==  24) return "W+";
  if (kf == -24) return "W-";
  if (kf ==  25) return "h" ;
  if (kf ==  35) return "H0";
  if (kf ==  36) return "H3";
  if (kf ==  37) return "H+";
  if (kf == -37) return "H-";
  if (kf ==  39) return "Gr";         // Graviton
  if (kf ==  41) return "R0";
  if (kf ==  42) return "LQ";

  if (kf ==  32) return  "t3";
  if (kf == -32) return  "T3";

  if (kf ==  1000021) return  "~g";            // gluino
  if (kf ==  1000001) return  "~dL";
  if (kf == -1000001) return  "~DL";
  if (kf ==  2000001) return  "~dR";
  if (kf == -2000001) return  "~DR";
  if (kf ==  1000002) return  "~uL";
  if (kf == -1000002) return  "~UL";
  if (kf ==  2000002) return  "~uR";
  if (kf == -2000002) return  "~UR";
  if (kf ==  1000003) return  "~sL";
  if (kf == -1000003) return  "~SL";
  if (kf ==  2000003) return  "~sR";
  if (kf == -2000003) return  "~SR";
  if (kf ==  1000004) return  "~cL";
  if (kf == -1000004) return  "~CL";
  if (kf ==  2000004) return  "~cR";
  if (kf == -2000004) return  "~CR";
  if (kf ==  1000005) return  "~b1";
  if (kf == -1000005) return  "~B1";
  if (kf ==  2000005) return  "~b2";
  if (kf == -2000005) return  "~B2";
  if (kf ==  1000006) return  "~t1";
  if (kf == -1000006) return  "~T1";
  if (kf ==  2000006) return  "~t2";
  if (kf == -2000006) return  "~T2";
  if (kf ==  1000011) return  "~eL";
  if (kf == -1000011) return  "~EL";
  if (kf ==  2000011) return  "~eR";
  if (kf == -2000011) return  "~ER";
  if (kf ==  1000012) return  "~ne";
  if (kf == -1000012) return  "~Ne";
  if (kf ==  1000013) return  "~mL";
  if (kf == -1000013) return  "~ML";
  if (kf ==  2000013) return  "~mR";
  if (kf == -2000013) return  "~MR";
  if (kf ==  1000014) return  "~nm";
  if (kf == -1000014) return  "~Nm";
  if (kf ==  1000015) return  "~l1";
  if (kf == -1000015) return  "~L1";
  if (kf ==  2000015) return  "~l2";
  if (kf == -2000015) return  "~L2";
  if (kf ==  1000016) return  "~nl";
  if (kf == -1000016) return  "~Nl";

  if (kf == 1000022) return  "~o1";             // neutralino chi_10
  if (kf == 1000023) return  "~o2";             // neutralino chi_20
  if (kf == 1000025) return  "~o3";             // neutralino chi_30
  if (kf == 1000035) return  "~o4";             // neutralino chi_40

  if (kf ==  1000024) return  "~1+";            // chargino chi_1+
  if (kf == -1000024) return  "~1-";            // chargino chi_1-

  if (kf ==  1000037) return  "~2+";            // chargino chi_2+
  if (kf == -1000037) return  "~2-";            // chargino chi_2-

  if (kf == 1000039) return  "~G";              // gravitino
  if (kf == 1000040) return  "~sS";             // sgoldstino S
  if (kf == 1000041) return  "~sP";             // sgoldstino P

  return "unkown";
}
