/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 2000, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "chep_limits.h"
#include "chep_crt/include/chep_crt.h"

#include "unix_utils.h"
#include "syst.h"
#include "files.h"
#include "lbl.h"

static shortstr chepversion;
static shortstr chepname;

char * version (void) {
return chepversion;
}

void setversion (shortstr ver) {
  strcpy (chepversion, ver);
}

char * getname (void) {
  sprintf(chepname, "CompHEP version %s   ", chepversion);
  return chepname;
}

void 
cheplabel (void)
{
  int key;
  char fname[STRSIZ];
  FILE *f;

  while (1)
    {
      goto_xy (20, 2);
      print ("Skobeltsyn Institute of Nuclear Physics");
      goto_xy (28, 3);
      print ("Moscow State University");
      goto_xy (28, 5);
      print ("%s", getname ());
      goto_xy (17, 6);
      goto_xy (15, 7);
      print ("a package for computation in high energy physics");
      scrcolor (FGmain, BGmain);
      goto_xy (16, 10);
      print ("Copyright (C) 2001-2009, CompHEP Collaboration:");
      goto_xy (15, 11);
      print ("E.Boos, V.Bunichev, M.Dubinin, L.Dudko, V.Edneral");
      goto_xy (15, 12);
      print ("V.Ilyin, A.Kryukov, V.Savrin, A.Sherstnev, A.Semenov");
      goto_xy (25, 14);
      print ("According to the license agreement");
      goto_xy (15, 15);
      print ("publications which result from using the program should");
      goto_xy (15, 16);
      print ("contain a reference to the article describing CompHEP.");
      goto_xy (25, 17);
      print ("Click the field below to get");
      chepbox (19, 18, 39, 20);
      chepbox (41, 18, 61, 20);
      scrcolor (Blue, LightGray);
      goto_xy (21, 19);
      print ("License Agreement");
      goto_xy (46, 19);
      print ("References");
      scrcolor (Black, LightGray);
      goto_xy (1, 22);
      print ("For contacts:");
      goto_xy (1, 23);
      print ("Alexander Kryukov<kryukov@theory.sinp.msu.ru>");
      goto_xy (47, 22);
      print ("Home page:");
      goto_xy (47, 23);
      print ("http://comphep.sinp.msu.ru");

      clearTypeAhead ();
      key = inkey ();
      if (key == KB_MOUSE && mouse_info.row == 19 &&
          mouse_info.col > 19 && mouse_info.col < 61)
        {
          if (mouse_info.col < 39)
            {
              sprintf (fname, "%s%cLicence.txt", pathtocomphep, f_slash);
              f = fopen (fname, "r");
              showtext (5, 1, 75, 1, "", f);
              fclose (f);
            }
          else if (mouse_info.col > 41)
            {
              sprintf (fname, "%s%cReference.txt", pathtocomphep, f_slash);
              f = fopen (fname, "r");
              showtext (5, 1, 75, 1, "", f);
              fclose (f);
            }
          else
            break;
        }
      else
        break;
    }
  clr_scr (FGmain, BGmain);
}
