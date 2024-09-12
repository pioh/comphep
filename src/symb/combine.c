/*
* Copyright (C) 2007-2009, CompHEP Collaboration
* Author: Alexander Sherstnev 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/parser.h"
#include "service2/include/unix_utils.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "polynom/include/polynom.h"
#include "polynom/include/tensor.h"
#include "polynom/include/spinor.h"
#include "chep_crt/include/chep_crt.h"

#include "physics.h"
#include "sos.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "pvars.h"
#include "chess.h"
#include "saveres.h"
#include "pre_read.h"
#include "reader0.h"
#include "reader_s.h"
#include "rfactor.h"
#include "process.h"
#include "process_core.h"

static void
keepdiagram (FILE * fsrc, FILE * target, catrec * cr)
{
  int i;
  int m;
  int delta;
  char * ff;

  fseek (fsrc, cr->factpos, SEEK_SET);
  delta = cr->rnumpos - cr->factpos;
  cr->factpos = ftell(target);
  ff = malloc (delta * sizeof(char));
  fread (ff, delta, 1, fsrc);
  f_write (ff, delta, 1, target);
  free(ff);

  fseek (fsrc, cr->rnumpos, SEEK_SET);
  delta = cr->denompos - cr->rnumpos;
  cr->rnumpos = ftell(target);
  ff = malloc (delta * sizeof(char));
  fread (ff, delta, 1, fsrc);
  f_write (ff, delta, 1, target);
  free(ff);

  fseek (fsrc, cr->denompos, SEEK_SET);
  FREAD1 (denrno, fsrc);
  for (i = 0; i < denrno; i++)
    {
      FREAD1 (denom[i].power, fsrc);
      FREAD1 (denom[i].mass, fsrc);
      FREAD1 (denom[i].width, fsrc);
      m = 0;
      do
        {
          FREAD1 (denom[i].momStr[m], fsrc);
        }
      while (denom[i].momStr[m++]);
    }

  cr->denompos = ftell(target);
  FWRITE1 (denrno, target);
  for (i = 0; i < denrno; i++)
    {
      FWRITE1 (denom[i].power, target);
      FWRITE1 (denom[i].mass, target);
      FWRITE1 (denom[i].width, target);
      m = 0;
      do
        {
          FWRITE1 (denom[i].momStr[m], target);
        }
      while (denom[i].momStr[m++]);
    }
}


int
combine (int nfiletot, FILE * menuq, FILE ** archives, FILE ** diaginfo, FILE ** catalogs)
{
  int ndiagr;
  int nfile;
  int naux, ndel, ncalc, nrest;
  long nrecord;
  csdiagram csd;
  catrec cr;

  midstr txt;

  FILE * archiv = fopen (ARCHIV_NAME, "wb");
  FILE * catalog = fopen (CATALOG_NAME, "wb");

  for (nfile = 0; nfile < nfiletot; ++nfile)
    {
      fseek (catalogs[nfile], 0, SEEK_SET);
    }

  fseek (archiv, 0, SEEK_SET);
  fseek (catalog, 0, SEEK_SET);
  for (nsub = 1; nsub <= subproc_sq; ++nsub)
    {
      rd_menu (menuq, 2, nsub, txt, &ndel, &ncalc, &nrest, &nrecord);
      naux = ndel + ncalc + nrest;
      for (ndiagr = 1; ndiagr <= naux; ++ndiagr)
        {
          for (nfile = 0; nfile < nfiletot; ++nfile)
            {
              fseek (diaginfo[nfile], sizeof (csd) * (nrecord + ndiagr - 1), SEEK_SET);
              FREAD1(csd, diaginfo[nfile]);

              if (2 == csd.status || -1 == csd.status || -2 == csd.status) {
                fseek (diaginfo[0], sizeof (csd) * (nrecord + ndiagr - 1), SEEK_SET);
                FWRITE1 (csd, diaginfo[0]);
                break;
              }

              if (1 == csd.status)
                {
                  int num = 0;
		  do {
		    fseek (catalogs[nfile], sizeof (cr) * num, SEEK_SET);
                    FREAD1(cr, catalogs[nfile]);
		    ++num;
		  } while (cr.ndiagr_ != ndiagr || cr.nsub_ != nsub);

                  if (cr.status != 0)
                    {
                      keepdiagram (archives[nfile], archiv, &cr);

                      cr.status = 1;
                      FWRITE1 (cr, catalog);

                      fseek (diaginfo[0], sizeof (csd) * (nrecord + ndiagr - 1), SEEK_SET);
                      FWRITE1 (csd, diaginfo[0]);

                      break;
                   }
                }
            }
        }
    }
  fclose (catalog);
  fclose (archiv);
  return 0;
}

int
compare_menuq (int nfiletot, int nsub, FILE ** menuqs)
{
  int i;
  int nproc;
  int ndel0, ncalc0, nrest0, recpos0;
  int ndel, ncalc, nrest, recpos;
  int ndiagr0;
  int ndiagr;
  int nfile;
  int ntot = 0;

  midstr txt0, txt;

  for (i = 0; i < nsub; ++i)
    {
      fseek (menuqs[0], (nsub - 1) * (60) + 2, SEEK_SET);
      if (6 != fscanf (menuqs[0], "%4d| %[^|]%*c%d|%d|%d|%d", &nproc, txt0, &ndel0, &ncalc0, &nrest0, &recpos0)) {
        return 0;
      }
      ndiagr0 = ncalc0 + nrest0;
      ntot += ndiagr0;
      for (nfile = 1; nfile < nfiletot; ++nfile)
        {
          fseek (menuqs[nfile], (nsub - 1) * (60) + 2, SEEK_SET);
          if (6 != fscanf (menuqs[nfile], "%4d| %[^|]%*c%d|%d|%d|%d", &nproc, txt, &ndel, &ncalc, &nrest, &recpos)) {
            return 0;
          }
          ndiagr = ncalc + nrest;
          if (strcmp(txt0, txt) || ndel != ndel0 || ndiagr != ndiagr0 || recpos != recpos0) {
            return 0;
          }
        }
    }

  return ndiagr0;
}
