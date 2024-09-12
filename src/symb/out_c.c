/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include <strings.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/unix_utils.h"
#include "service2/include/tptcmac.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"
#include "service2/include/files.h"
#include "service2/include/lbl.h"
#include "chep_crt/include/chep_crt.h"
#include "polynom/include/polynom.h"

#include "physics.h"
#include "procvar.h"
#include "pvars.h"
#include "diaprins.h"
#include "optimise.h"
#include "l_string.h"
#include "reader_c.h"
#include "out_service.h"
#include "saveres.h"
#include "denominators.h"
#include "process.h"
#include "process_core.h"
#include "cweight.h"
#include "out_c.h"

#include "prepdiag.h"

int *diag_hashtopology = NULL;
int *subprocess_type = NULL;

static int myfunc_key = 0; 

#define procinfptr struct procinfrec *
typedef struct procinfrec
{
  procinfptr next;
  unsigned tot;
  unsigned firstdiagpos;
  prtclsarray p_name;
  int p_masspos[MAXINOUT];
}
procinfrec;
#undef procinfptr

typedef struct procinfrec *procinfptr;

static procinfptr inf, inftmp;          /*  information about subProcess  */
static unsigned ndiagrtot, diagrcount;  /*  statictics  */
static int nvars, nfunc;                /*  statictics  */
static marktp heapbeg;                  /*  for RELEASE  */
static int nden_w, nden_0, nsub1;       /* from writesubprocess */
static int cBasisPower;
static int nC, *cChains = NULL;
static long *cCoefN, *cCoefD;
longstr rsltpath = "";


static void
clearstatistic (void)
{
  int i;
  for (i = 17; i < 24; i++)
    {
      goto_xy (1, i);
      clr_eol ();
    }
}

int
get_diagnumber (void)
{
  int ndiag = 0;
  catrec cr;
  FILE *cat;

  cat = fopen (CATALOG_NAME, "r");
  while (FREAD1 (cr, cat) == 1)
    ndiag++;
  fclose (cat);

  return ndiag;
}

static void
init_stat (void)
{
  goto_xy (1, 17);
  scrcolor (Yellow, Blue);
  print (" C Source Codes \n");
  scrcolor (Red, BGmain);
  print (" Process..........\n");
  print (" Total diagrams...\n");
  print (" Processed........\n");
  print (" Current..........\n");
  scrcolor (Yellow, Blue);
  print (" Press Esc to stop    ");
  scrcolor (Black, BGmain);
  goto_xy (20, 18);
  print ("%s", getProcessch ());
  goto_xy (20, 19);
  print ("%4u", ndiagrtot);
  goto_xy (20, 20);
  print ("   0");
  scrcolor (Yellow, BGmain);
  goto_xy (20, 21);
  print ("   1");
  scrcolor (Yellow, BGmain);
}


static void
writestatistic (int ndcur, int ndtot)
{
  scrcolor (Black, BGmain);
  goto_xy (20, 19);
  print ("%4u", ndtot);
  goto_xy (20, 20);
  print ("%2u (%%)", (((ndcur - 1) * 100) / ndtot));
  goto_xy (20, 21);
  print ("%4u", ndcur);
}

/*************  big additions in this function for #-mdl ********/
static int
get_hash_info (unsigned ndiagr)
{
  csdiagram csdiagr;
  int lf, nv, qv, ql, if1, if2;
  int fqf, iqf;
  int xqQ, iqQ;
  char *xqname;
  char qmom;

  fseek (diagrq, ndiagr * sizeof (csdiagr), SEEK_SET);
  FREAD1 (csdiagr, diagrq);
  transfdiagr (&csdiagr, &vcs);
  preperdiagram ();
  xqQ = 0;                      /* not a PDF# model */
  iqQ = 0;                      /* not a qQ#-mdl */
  iqf = 1;                      /* scattering type by default */
  fqf = 0;                      /* no final QCD #-loops */
  for (lf = 1; lf <= nloop; lf++)
    {
      if1 = 0;
      if2 = 0;
      for (nv = 1; nv <= fermloops[lf - 1].len; nv++)
        {
          qv = fermloops[lf - 1].vv[nv - 1];
          ql = fermloops[lf - 1].ll[nv - 1];
          qmom = vcs.vertlist[qv - 1][ql - 1].moment;
          xqname = prtclbase[vcs.vertlist[qv - 1][ql - 1].partcl - 1].name;
          if (strchr (xqname, 35) != 0)
            xqQ = 1;
          if (strcmp (xqname, "q#") == 0 || strcmp (xqname, "Q#") == 0)
            iqQ = 1;
          if (abs (qmom) == 1)
            if1 = 1;
          if (getnin () == 2)
            {
              if (abs (qmom) == 2)
                if2 = 1;
            }
          if (if1 * if2 == 1)
            iqf = 0;
        }
      if ((strchr (xqname, 35) != 0) && if1 == 0 && if2 == 0)
        fqf++;
    }

  if (xqQ)
    {
      if (getnin () == 2)
        writeF ("/* PDF #-mdl: annih/scatt (0/1) type:    %d */\n", iqf);
      writeF ("/* PDF #-mdl: number of final QCD #-loops = %d */\n", fqf);
      if (fqf)
        writeF ("int i;\n");
    }
  subprocess_type[nsub1] = xqQ + iqQ;
  if (!xqQ)
    return 1;
  return 2 * fqf + iqf;
}

static void
writpict (unsigned ndiagr)
{
  csdiagram csdiagr;
  fseek (diagrq, ndiagr * sizeof (csdiagr), SEEK_SET);
  FREAD1 (csdiagr, diagrq);
  transfdiagr (&csdiagr, &vcs);
  writeF ("/*\n");
  DiagramToOutFile (&vcs, 1, ' ');
  writeF ("*/\n");
}

static void
labl (void)
{
  writeF ("/************************************************\n");
  writeF ("*    %s*\n", getname ());
  writeF ("*------------------------------------------------\n");
  writeF ("* Copyright (C) 2001-2009, CompHEP Collaboration*\n");
  writeF ("************************************************/\n");
}

  /* =========== Preliminary  calculations ================ */

static void
calc_nvars_nfunc (void)
{
  int k;

  nvars = 0;
  nfunc = 0;

  for (k = 1; k <= nmodelvar; k++)
    {
      if (vararr[k].used)
        {
          if (modelvars[k].func)
            nfunc++;
          else
            nvars++;
        }
    }
}


static void
prepareprocinform (void)
{
  int i, k;
  int ndel, ncalc, nrest;
  long recpos;
  char txt[STRSIZ];
  csdiagram csd;
  char mass[6];
  int nsubs;
  int numtot = getntot ();

  inf = NULL;
  menuq = fopen (MENUQ_NAME, "rb");
  for (nsubs = 1; nsubs <= subproc_sq; nsubs++)
    {
      inftmp = inf;
      inf = (procinfptr) getmem_ ((unsigned) sizeof (procinfrec));
      inf->next = inftmp;
      rd_menu (menuq, 2, nsubs, txt, &ndel, &ncalc, &nrest, &recpos);
      inf->firstdiagpos = recpos;
      getprtcls (txt, inf->p_name);
      for (i = 0; i < numtot; i++)
        {
          int nn = locateinbase (inf->p_name[i]);
          strcpy (mass, prtclbase[nn - 1].massidnt);
          if (strcmp (mass, "0"))
            for (k = 1; strcmp (modelvars[k].varname, mass); k++);
          else
            k = 0;
          inf->p_masspos[i] = k;
        }
      for (i = numtot; i < MAXINOUT; i++)
        {
          strcpy (inf->p_name[i], "***");
          inf->p_masspos[i] = 0;
        }
      fseek (diagrq, recpos * sizeof (csdiagram), SEEK_SET);
      inf->tot = 0;
      for (i = 1; i <= ndel + ncalc + nrest; i++)
        {
          FREAD1 (csd, diagrq);
          if (1 == csd.status || 0 == csd.status || 4 == csd.status)
            ++(inf->tot);
        }
    }
  nsubs--;
  fclose (menuq);
  revers ((pointer *) & inf);
}


static int 
koeff_va (void)
{
  int l;
 
  for (l = 0; l <= nmodelvar; l++)
    {
      if (strcmp (modelvars[l].varname, "koeff") == 0)
	{
	  return l;          
	}
    }
 
  return 0;
}


void
scanvars (FILE * f, int mode)   /* f - just for the sesssion.dat writing!!! */
{
  int num = 0;
  int k;

  for (k = 1; k <= nmodelvar; k++)
    {
      varlist modl = modelvars + k;
      if (vararr[k].used && !modl->func)
        {
          switch (mode)
            {
            case 2:
              {
                writeF ("\n,%E", modl->varvalue);
                num++;
              }
              break;

            case 3:
              writeF ("\n,\"%s\"", modl->varname);
              break;
            case 4:
              fprintf (f, "%10s = %.15E\n", modl->varname, modl->varvalue);
              break;
            }
        }
    }

  for (k = 1; k <= nmodelvar; k++)
    {
      varlist modl = modelvars + k;
      if (vararr[k].used && modl->func)
        {
          switch (mode)
            {
            case 2:
              writeF ("\n,%E", modl->varvalue);
              num++;
              break;
            case 3:
              writeF ("\n,\"%s\"", modl->varname);
              break;
            }
        }
    }
}


/* =========== Common blocks Emit =========== */
static void
common (char *extern_)
{
  if (strcmp (extern_, "extern"))
    {
      writeF (" double va[%d] ={", nvars + nfunc + 1);
      if (getnin () > 1)
        writeF ("%E", getsqrtS ());
      else
        writeF ("0");
      scanvars (NULL, 2); 
    /*  findkoeffvar (NULL, 3); */
      writeF ("};\n");
    }
  else
    writeF ("%s double va[%d];\n", extern_, nvars + nfunc + 1);
}


static void
writesubroutineinit (void)
{
  int l;
  char *ss;
  char d_type[20]; 

  writeF ("#include\"extern.h\"\n");
  ext_h = fopen (scat ("%s/extern.h", rsltpath), "w");
  strcpy (d_type, "double");
  writeF ("int calcFunc(void)\n{\n");
  writeF ("int err=0;\n");

  for (l = 1; l <= nmodelvar; l++)
    {
      if (vararr[l].used && modelvars[l].func)
        {
          ss = (char *) readExpression (modelvars[l].func, rd_c, act_c, free);
          writeF ("   %s=%s;\n", vararr[l].alias, ss + 3);
          if( (strncmp(ss + 3,"myfunc",6)==0) && myfunc_key==0 ) myfunc_key=l; 
          free (ss);
          writeF ("   if(!finite(%s)) return FUCTION_ERROR;\n",
                  vararr[l].alias);
        }
    }
  writeF ("if(err) return 1; else return 0;\n}\n");
  
/*-------------- for myfunc ------------------------*/
  if(myfunc_key > 0) 
   {
    writeF ("\n");
    writeF ("int calcMomFunc(void)\n{\n");
    writeF ("int err=0;\n");
    
    for (l = myfunc_key; l <= nmodelvar; l++)
     {
      if (vararr[l].used && modelvars[l].func)
        {
          ss = (char *) readExpression (modelvars[l].func, rd_c, act_c, free);
          writeF ("   %s=%s;\n", vararr[l].alias, ss + 3); 
          free (ss);
          writeF ("   if(!finite(%s)) return FUCTION_ERROR;\n",
                  vararr[l].alias);
        }
      }
     writeF ("if(err) return 1; else return 0;\n}\n");
    } 
/*---------------------------------------------------*/
  
  fclose (ext_h);
}


static void
onediagram (FILE * fsrc, deninforec * dendescript, catrec cr)
{
  int i, k, koeffnum;
  int addpr;
  int deg1, deg2, nConst;
  long pos_c;
  marktp bh;
  varptr totnum, totdenum, rnum;
  csdiagram csdiagr;

  mark_ (&bh);
  settmpNameMax (0);
  initfortwriting ('c');
  initdegnames ();

  ++diagrcount;

  outFileOpen (scat ("%s/f%d.c", rsltpath, diagrcount));
  labl ();
  writpict (cr.ndiagr_ + inftmp->firstdiagpos - 1);
  diag_hashtopology[diagrcount] =
    get_hash_info (cr.ndiagr_ + inftmp->firstdiagpos - 1);

  writeF ("#include <math.h>\n");
  writeF ("#include <stdio.h>\n");
  writeF ("extern double *Q0, *Q1, *Q2;\n");
  common ("extern");
  writeF ("#include\"out_ext.h\"\n");
  writeF ("#include\"out_int.h\"\n");
  writeF ("FNN F%d;\n", diagrcount);
  writeF ("double F%d(void)\n{\n", diagrcount);
  writeF ("double TOTNUM,TOTDEN,RNUM,result;\n");
  pos_c = ftell (outFile);
  writeF ("%80s\n", "");

  fseek (fsrc, cr.factpos, SEEK_SET);
  readvardef (fsrc);
  readpolynom (fsrc, &totnum);
  readpolynom (fsrc, &totdenum);
  clearvardef ();

  fseek (fsrc, cr.rnumpos, SEEK_SET);
  readvardef (fsrc);
  readpolynom (fsrc, &rnum);
  clearvardef ();
  if(myfunc_key==0) writeF ("if(calcCoef[%d])\n{\n", cr.nsub_);
  nConst = write_const ();
  if(myfunc_key==0) writeF ("}\n");
  deg1 = cleardegnames ();
  initdegnames ();
  fortwriter ("TOTNUM", totnum);
  fortwriter ("TOTDEN", totdenum);
  fortwriter ("RNUM", rnum);

  if(cr.status == 4) 
   {
    koeffnum = koeff_va();
    writeF ("result=%s*RNUM*(TOTNUM/TOTDEN)",vararr[koeffnum].alias);
   }
  else writeF ("result=RNUM*(TOTNUM/TOTDEN)");

  for (i = 0; i < dendescript->tot_den; ++i)
    {
      int numm = dendescript->denarr[i].width ?
        dendescript->denarr[i].order_num :
        dendescript->denarr[i].order_num + nden_w;
      if (dendescript->denarr[i].power == 1)
        writeF ("*Q1[%d]", numm);
      else
        writeF ("*Q2[%d]", numm);
    }

  for (k = 1; k <= nden_w; k++)
    {
      addpr = TRUE;
      for (i = 0; i < dendescript->tot_den; ++i)
        {
          if (dendescript->denarr[i].width &&
              k == dendescript->denarr[i].order_num)
            addpr = FALSE;
        }
      if (addpr)
        {
          writeF ("*Q0[%d]", k);
        }
    }

  writeF (";\n");
  writeF
    (" if(result>Fmax) Fmax=result; else if(result<-Fmax) Fmax=-result;\n");

  fseek (diagrq, (cr.ndiagr_ + inftmp->firstdiagpos - 1) * sizeof (csdiagr), SEEK_SET);
  FREAD1 (csdiagr, diagrq);

  if (cBasisPower
      && generateColorWeights (&csdiagr, cBasisPower, nC, cChains, cCoefN,
                               cCoefD))
    {
      writeF (" if(color_weights)\n {\n");
      for (k = 0; k < cBasisPower; k++)
        if (cCoefN[k])
          writeF ("  color_weights[%d] += result*(%d)/(%d);\n",
                  k, cCoefN[k], cCoefD[k]);
      writeF (" }\n");
    }

  writeF (" return result;\n");
  writeF ("}\n");
  deg2 = cleardegnames ();
  deg1 = MAX (deg1, deg2);

  if (nConst || deg1 || gettmpNameMax ())
    {
      fseek (outFile, pos_c, SEEK_SET);
      if (nConst)
        writeF ("static double C[%d];", nConst);
      if (deg1)
        writeF ("double S[%d];", deg1);
      if (gettmpNameMax ())
        writeF ("double tmp[%d];", gettmpNameMax ());
    }
  outFileClose ();
  release_ (&bh);
}


static void
writesubprocess (int nsub, FILE * file, int *breaker)
{
  denlist den_, den_tmp;
  int i;
  deninforec dendescript;
  FILE *fd;                     /* file of (deninforec)  */
  char fd_name[STRSIZ];
  marktp mem_start;

  nsub1 = nsub;

  outFileOpen (scat ("%s/d%d.c", rsltpath, nsub));
  labl ();
  writeF ("#include<math.h>\n");
  writeF ("#define real double\n");
  writeF ("#include\"out_int.h\"\n");
  writeF ("#include\"out_ext.h\"\n");

  writeF ("extern real *Q0,*Q1,*Q2;\n");
  common ("extern");

  writeF ("DNN  d_%d;\n", nsub);

  writeF ("int d_%d(real * momenta)\n{", nsub);

  writeF ("int I,err=0;\n");
  writeF ("real s0max=0;\n");

  sprintf (fd_name, "%stmp%cden.inf", pathtouser, f_slash);
  fd = fopen (fd_name, "wb");

  mark_ (&mem_start);
  denominatorStatistic (file, nsub, &nden_w, &nden_0, &den_, fd);
  fclose (fd);

  if (nden_w != 0)
    {
      writeF (" real DMASS[%d],DWIDTH[%d];\n", nden_w + 1, nden_w + 1);
      writeF (" if(Q0!=NULL) free(Q0);\n");
      writeF (" Q0=(real*)malloc(sizeof(real)*%d);\n", nden_w + 1);
    }

  writeF (" for(I=0;I<nin_;I++) s0max+=momenta[4*I];\n");
  writeF ("s0max=computer_eps*s0max*s0max;\n");

  if (nden_w + nden_0 != 0)
    {
      writeF (" if(Q1!=NULL) free(Q1);\n");
      writeF (" if(Q2!=NULL) free(Q2);\n");
      writeF (" Q1=(real*)malloc(sizeof(real)*%d);\n", nden_w + nden_0 + 1);
      writeF (" Q2=(real*)malloc(sizeof(real)*%d);\n", nden_w + nden_0 + 1);
    }

  while (den_ != NULL)
    {
      int m = 0;
      i = den_->order_num;
      if (den_->width)
        {
          writeF ("DMASS[%d]=%s;\n", i, vararr[den_->mass].alias);
          writeF ("DWIDTH[%d]=%s;\n", i, vararr[den_->width].alias);
        }
      else
        i += nden_w;
      writeF ("Q1[%d]=", i);
      if (den_->mass != 0)
        writeF ("%s*%s", vararr[den_->mass].alias, vararr[den_->mass].alias);

      writeF ("- sqrMom(\"");
      while (den_->momStr[m])
        writeF ("\\%o", den_->momStr[m++]);
      writeF ("\",momenta);\n");
      den_tmp = den_;
      den_ = den_->next;
    }

  release_ (&mem_start);
  if (nden_w > 0)
    {
      writeF ("for ( I=1;I<= %d;I++)\n {\n  real q=Q1[I],w;\n", nden_w);
      writeF
        ("  if(gwidth==2) w=(DMASS[I]-Q1[I]/DMASS[I])*DWIDTH[I]; else w=DMASS[I]*DWIDTH[I];\n");
      writeF ("  Q2[I]=1/(q*q+w*w);\n");
      writeF ("  if (gwidth==1) Q0[I]=Q2[I]*q*q;else Q0[I]=1; \n");
      writeF ("  Q1[I]=Q2[I]*q;\n");
      writeF (" }\n");
    }

  if (nden_0 > 0)
    {
      writeF ("for ( I=%d; I<=%d;I++)\n  {\n", 1 + nden_w, nden_w + nden_0);
      writeF
        ("  if((Q1[I]>0? Q1[I]:-Q1[I]) < 10*s0max) err=DENOMINATOR_ERROR;\n");
      writeF ("  if(!Q1[I]) Q1[I]=s0max;\n");
      writeF ("  Q1[I]=1/Q1[I];\n");
      writeF ("  Q2[I]=Q1[I]*Q1[I];\n");
      writeF ("}\n");
    }
  /* ===================== */
  writeF ("return err;\n");
  writeF ("}\n");

  outFileClose ();
  fd = fopen (fd_name, "rb");
  *breaker = FALSE;
  while (FREAD1 (dendescript, fd) == 1)
    {
      catrec cr;
      if (escpressed ())
        {
          *breaker = TRUE;
          break;
        }
      fseek (catalog, dendescript.cr_pos, SEEK_SET);
      FREAD1 (cr, catalog);
      if ( (1 == cr.status) || (4 == cr.status) )
        {
          onediagram (file, &dendescript, cr);
        }
      else
        {
          ++diagrcount;
          diag_hashtopology[diagrcount] = 0;
        }
      writestatistic (diagrcount, ndiagrtot);
    }

  fclose (fd);
  unlink (fd_name);
}


static void
make_pinf (void)
{
  int i;
  int numtot = getntot ();

  writeF ("int pinf_(int nsub,int nprtcl,char* pname, double* pmass)\n{\n");
  writeF ("int n;\n");
  writeF ("char names[%d][%d][7]={\n", subproc_sq, numtot);

  inftmp = inf;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      writeF (" {");
      for (i = 1; i <= numtot; i++)
        {
          writeF ("\"%s\"", inftmp->p_name[i - 1]);
          if (i == numtot)
            if (nsub == subproc_sq)
              writeF ("}\n};\n");
            else
              writeF ("},\n");
          else
            writeF (",");
        }
      inftmp = inftmp->next;
    }
  writeF ("int nvalue[%d][%d]={\n", subproc_sq, numtot);

  inftmp = inf;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      writeF (" {");
      for (i = 1; i <= numtot; i++)
        {
          int k = inftmp->p_masspos[i - 1];
          if (k)
            sscanf (vararr[k].alias, "va[%d]", &k);
          writeF ("%d", k);
          if (i != numtot || nsub != subproc_sq)
            writeF (",");
          if (i == numtot)
            writeF ("},\n");
        }
      if (nsub == subproc_sq)
        writeF ("};\n");
      inftmp = inftmp->next;
    }

  writeF ("if  (nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return 1;\n", subproc_sq, numtot);
  writeF ("if(pname) strcpy(pname,names[nsub-1][nprtcl-1]);\n");
  writeF ("if(pmass)\n{\n");
  writeF ("  n=nvalue[nsub-1][nprtcl-1];\n");
  writeF ("if (n>%d) if (calcFunc()) return FUCTION_ERROR;\n", nvars);
  writeF ("if (n==0) *pmass=0; else *pmass=va[n];\n");
  writeF ("if (*pmass<0) (*pmass)=-(*pmass);\n");
  writeF ("}\n");
  writeF ("return 0;\n");
  writeF ("}\n");
}


static void
make_infbasis (void)
{
  int i, j;
  int pcolor[MAXINOUT];
  int numtot = getntot ();

  writeF ("void cStrings(int nsub,int *nC, int * power, int **  chains)\n{\n");
  writeF ("   switch(nsub)\n   {\n");

  inftmp = inf;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      writeF ("   case %d : ", nsub);
      for (i = 0; i < numtot; i++)
        {
          pcolor[i] = prtclbase[locateinbase (inftmp->p_name[i]) - 1].cdim;
          if (i < getnin ())
            {
              if (pcolor[i] == 3)
                pcolor[i] = -3;
              else if (pcolor[i] == -3)
                pcolor[i] = 3;
            }
        }

      infCbases (numtot, pcolor, &nC, &cBasisPower, &cChains);
      if (cBasisPower)
        {
          writeF ("\n     { static int cc[%d]=\n       {\n", 2 * nC * cBasisPower);
          for (i = 0; i < cBasisPower; i++)
            {
              writeF ("       ");
              for (j = 0; j < nC; j++)
                {
                  writeF (" %d,%d", cChains[2 * i * nC + 2 * j] + 1,
                          cChains[2 * i * nC + 2 * j + 1] + 1);
                  if (i == cBasisPower - 1 && j == nC - 1)
                    writeF ("\n       };");
                  else
                    writeF (",");
                }
              writeF ("\n");
            }

          writeF("       *nC=%d; *power=%d; *chains=cc;\n     }\n     break;\n",
             nC, cBasisPower);
          free (cChains);
          cChains = NULL;
        }
      else
        writeF ("   *nC=0; *power=0; *chains=NULL; break;\n");

      inftmp = inftmp->next;
    }
  writeF ("   default: *nC=0; *power=0; *chains=NULL;\n");
  writeF ("   }\n");

  writeF ("}\n\n");
}


static void
make_fosimple (void)
{
  unsigned i, j, totcount;
  int max_top_class = 0;

  for (i = 1; i <= ndiagrtot; i++)
    {
      if (max_top_class < diag_hashtopology[i])
        max_top_class = diag_hashtopology[i];
    }

  writeF ("\nstatic double smpl(int nsub, double * momenta,int * err)\n{\n");
  writeF (" double q;\n");
  writeF (" double ans, ans0=0.0, ans1=0.0");
  for (i = 2; i <= max_top_class; i++)
    writeF (", ans%i=0", i);
  writeF (";\n");
  inftmp = inf;
  totcount = 0;
  
  if(myfunc_key > 0) writeF ("\n calcMomFunc();\n\n");

  writeF ("switch(nsub)\n{\n");
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      if (inftmp->tot != 0)
        {
          writeF ("  case %d:\n    *err=*err|d_%d(momenta);\n", nsub, nsub);
          writeF ("    sprod_(momenta);");
          for (j = 0; j <= max_top_class; j++)
            {
              writeF ("\n    ans%i=0", j);
              for (i = 1; i <= inftmp->tot; i++)
                {
                  if (j == diag_hashtopology[totcount + i])
                    {
                      writeF ("+F%d()", totcount + i);
                    }
                }
              writeF (";");
            }
          totcount += inftmp->tot;
          if (subprocess_type[nsub])
            {
              int mtc = max_top_class / 2;
              if (mtc)
                writeF ("\n    ans0+=", j);
              for (i = 1; i <= mtc; i++)
                writeF ("+%i*(ans%i", 2 * subprocess_type[nsub], 2 * i);
              for (i = 1; i <= mtc; i++)
                writeF (")");
              mtc = (max_top_class - 1) / 2;
              if (mtc)
                writeF (";\n    ans1+=", j);
              for (i = 1; i <= mtc; i++)
                writeF ("+%i*(ans%i", 2 * subprocess_type[nsub], 2 * i + 1);
              for (i = 1; i <= mtc; i++)
                writeF (")");
              writeF (";\n");
            }
          writeF ("\n    break;\n");
        }
      inftmp = inftmp->next;
    }
  writeF ("}\n");
  if (getnin () == 2)
    {
      writeF ("  if (!strfun_calc) {\n");
      writeF ("    strfun_calc = 1;\n");
      writeF ("    q = qcd_Scale_chep();\n");
      writeF ("    if (ans0) xstr0 = strfun_(0, xbjo[0], xbjo[1], q);\n");
      writeF ("    if (ans1) xstr1 = strfun_(1, xbjo[0], xbjo[1], q);\n");
      writeF ("  }\n");
    }
  else
    {
      writeF ("  xstr0 = 1.0;\n");
      writeF ("  xstr1 = 1.0;\n");
    }
  writeF ("  ans = ans0 + ans1;\n");
  writeF ("  if(!(*err) && 10000 * Fmax * computer_eps > (ans > 0 ? ans : -ans)) * err = 1;\n");
  writeF ("  return ans0 * xstr0 + ans1 * xstr1;\n}\n");
}


static void
make_asgn (void)
{
  writeF ("int asgn_ (int numvar, double newval)\n{\n");
  if (nvars > 0)
    {
      writeF ("  if (numvar < 0|| numvar>%d) return 1;\n", nvars);
      writeF ("    va[numvar]=newval;\n");
    }
  writeF ("   return 0;\n");
  writeF ("}\n\n");
}


static void
make_vinf (void)
{
  writeF ("int vinf_ (int numvar, char * name, double * val)\n{\n");
  writeF ("char names[%d][10]={\n\"Sqrt(S)\"", nvars + nfunc + 1);
  scanvars (NULL, 3);
  writeF ("};\n");
  writeF ("   if (numvar < 0 || numvar > %d) return 1;\n", nvars + nfunc);
  writeF ("   if (name) strcpy (name, names[numvar]);\n");
  writeF ("   if (val) *val=va[numvar];\n");
  writeF ("   return 0;\n}\n\n");
}


static void
zeroHeep (void)
{
  goto_xy (1, 1);
  print ("Heep is empty!!!");
  inkey ();
  exit (0);
}

void
get_first_subprocname (midstr subprocname)
{
  int i;
  int thenin = 1;
  int numin = getnin ();
  int numtot = getntot ();
  vshortstr tmp;
  FILE * archiv = fopen (ARCHIV_NAME, "rb");

  outputLanguage = 'c';
  catalog = fopen (CATALOG_NAME, "rb");
  diagrq = fopen (DIAGRQ_NAME, "rb");
  fseek (catalog, 0, SEEK_END);
  ndiagrtot = ftell (catalog) / sizeof (catrec);

  memerror = zeroHeep;
  mark_ (&heapbeg);
  initvararray (archiv, 0, outputLanguage);

  firstVar = nmodelvar;
  if (!strcmp (modelvars[firstVar].varname, strongconst))
    {
      firstVar--;
    }

  prepareprocinform ();
  calc_nvars_nfunc ();
  inftmp = inf;
  nsub = 1;
  if (1 == numin)
    {
      sprintf (subprocname, "%s -> ", inftmp->p_name[0]);
    }
  else
    {
      sprintf (subprocname, "%s,%s -> ", inftmp->p_name[0], inftmp->p_name[1]);
      thenin = 2;
    }
  for (i = thenin; i < numtot - 1; ++i)
    {
      sprintf (tmp, "%s,", inftmp->p_name[i]);
      strcat (subprocname, tmp);
    }
  sprintf (tmp, "%s", inftmp->p_name[i]);
  strcat (subprocname, tmp);

  fclose (catalog);
  fclose (diagrq);
  fclose (archiv);
}

int
c_prog (void)
{
  int breaker;
  int i;
  int numtot = getntot ();
  FILE * archiv = fopen (ARCHIV_NAME, "rb");

  outputLanguage = 'c';
  catalog = fopen (CATALOG_NAME, "rb");
  diagrq = fopen (DIAGRQ_NAME, "rb");
  fseek (catalog, 0, SEEK_END);
  ndiagrtot = ftell (catalog) / sizeof (catrec);

  memerror = zeroHeep;
  mark_ (&heapbeg);
  initvararray (archiv, 0, outputLanguage);

/* Initialisation part */
  firstVar = nmodelvar;
  if (!strcmp (modelvars[firstVar].varname, strongconst))
    firstVar--;
  prepareprocinform ();
  calc_nvars_nfunc ();
  sprintf (rsltpath, "%sresults", pathtouser);

/* service.c */
  outFileOpen (scat ("%s/service.c", rsltpath));
  labl ();
  writeF ("#include<math.h>\n");
  writeF ("#include\"out_int.h\"\n");
  writeF ("#include\"out_ext.h\"\n");
  writeF ("int gwidth = 0;\n");
  common ("");
  writeF ("const int nin_ = %d;\n\n", getnin ());
  writeF ("const int nout_ = %d;\n\n", getnout ());
  writeF ("const int nprc_ = %d;\n\n", subproc_sq);
  make_pinf ();
  writeF ("const int nvar_ = %d;\n\n", nvars);
  writeF ("const int nfunc_ = %d;\n\n", nfunc);
  make_vinf ();
  make_asgn ();
  make_infbasis ();
  writesubroutineinit ();
  outFileClose ();

/* f*.c files */
  diag_hashtopology = malloc ((ndiagrtot + 1) * sizeof (int));
  subprocess_type = malloc ((subproc_sq + 1) * sizeof (int));
  diagrcount = 0;
  inftmp = inf;
  init_stat ();
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      int colors[MAXINOUT];

      if (inftmp->tot != 0)     /*  this subprocess in Archive  */
        {
          for (i = 0; i < numtot; i++)
            {
              colors[i] = prtclbase[locateinbase (inftmp->p_name[i]) - 1].cdim;
            }
          for (i = 0; i < getnin (); ++i)
            if (colors[i] == 3)
              colors[i] = -3;
            else if (colors[i] == -3)
              colors[i] = 3;
          infCbases (numtot, colors, &nC, &cBasisPower, &cChains);
          if (cBasisPower)
            {
              cCoefN = malloc (cBasisPower * sizeof (long));
              cCoefD = malloc (cBasisPower * sizeof (long));
            }
          writesubprocess (nsub, archiv, &breaker);
          if (breaker)
            goto exi;
          if (cBasisPower)
            {
              if (cChains)
                {
                  free (cChains);
                  cChains = NULL;
                }
              free (cCoefN);
              free (cCoefD);
            }
        }
      inftmp = inftmp->next;
    }

/* sqme.c */
  outFileOpen (scat ("%s/sqme.c", rsltpath));
  labl ();
  writeF ("#include <math.h>\n");
  writeF ("#include <stdio.h>\n");
  common ("extern");
  writeF ("#include\"out_int.h\"\n");
  writeF ("#include\"out_ext.h\"\n");
  if (getnin () == 2)
    {
      writeF ("#include\"../src/num/include/alphas2.h\"\n");    /* #-mdl */
      writeF ("#include\"../src/num/include/alphas_menu.h\"\n");/* #-mdl */
      writeF ("#include\"../src/num/include/strfun.h\"\n");     /* #-mdl */
      writeF ("#include\"../src/num/include/q_kin.h\"\n");      /* #-mdl */
    }
  writeF ("char  processch[] = \"%s\";\n", getProcessch ());
  writeF ("double DP[%d];\n", ((MAXINOUT * (MAXINOUT - 1) / 2)));
  writeF ("static double xstr0,xstr1;\n");
  writeF ("static int strfun_calc;\n");
  writeF ("extern DNN ");
  for (i = 1; i < subproc_sq; i++)
    writeF ("d_%d,", i);
  writeF ("d_%d;\n", subproc_sq);
  fseek (catalog, 0, SEEK_END);
  ndiagrtot = ftell (catalog) / sizeof (catrec);
  writeF ("extern FNN ");
  for (i = 1; i < ndiagrtot; i++)
    writeF ("F%d,", i);
  writeF ("F%d;\n", ndiagrtot);
  writeF ("static void sprod_(double*);\n");
  make_fosimple ();
  writeF ("#include\"sqme0.c\"\n");
  outFileClose ();

exi:
  free (diag_hashtopology);
  clearstatistic ();
  fclose (catalog);
  fclose (archiv);
  fclose (diagrq);
  release_ (&heapbeg);
  return !breaker;
}


/*--------------------From FORM ---------------------------------------------------*/
/*--------------------From FORM ---------------------------------------------------*/
/*--------------------From FORM ---------------------------------------------------*/
/*--------------------From FORM ---------------------------------------------------*/
/*--------------------From FORM ---------------------------------------------------*/
/*--------------------From FORM ---------------------------------------------------*/
static void
calc_nvars_nfunc_F (void)
{
  int k;
  nvars = 0;
  nfunc = 0;
  for (k = 1; k <= nmodelvar; k++)
    {
      {
        if (modelvars[k].func)
          nfunc++;
        else
          nvars++;
      }
    }
}



static void
scanvars_F (int mode)
{
  int num = 0;
  int k;

  for (k = 0; k < nmodelvar; ++k)
    {
      varlist modl = modelvars + k + 1;
      if (!modl->func)
        {
          switch (mode)
            {
            case 2:
              {
                writeF ("\n,%E", modl->varvalue);
                num++;
              }
              break;

            case 3:
              writeF ("\n,\"%s\"", modl->varname);
              break;
            }
        }
    }

  for (k = 0; k < nmodelvar; ++k)
    {
      varlist modl = modelvars + k + 1;
      if (modl->func)
        {
          switch (mode)
            {
            case 2:
              {
                writeF ("\n,%E", modl->varvalue);
                num++;
              }
              break;

            case 3:
              writeF ("\n,\"%s\"", modl->varname);
              break;
            }
        }
    }

}


static void
common_F (char *extern_)
{
  if (strcmp (extern_, "extern"))
    {
      writeF (" real va[%d] ={", nvars + nfunc + 1);
      if (getnin () > 1)
        writeF ("%E", getsqrtS ());
      else
        writeF ("0");
      scanvars_F (2);
      writeF ("};\n");

    }
  else
    writeF ("%s real va[%d];\n", extern_, nvars + nfunc + 1);

}


void *
rd_cF (char *s)
{
  char *p;
  int l;
  p = m_alloc (40);
  if ('0' <= s[0] && s[0] <= '9')
    sprintf (p, "MN|(real)%s", s);
  else
    {
      for (l = 1; l <= nmodelvar; l++)
        {
          if (!strcmp (s, modelvars[l].varname))
            {
              sprintf (p, "MR|%s", modelvars[l].varname);
              return (void *) p;
            }
        }
    }
  return (void *) p;
}


static void
writesubroutineinit_F (void)
{
  int l;
  char *ss;
  char d_type[20];

  writeF ("#include\"extern.h\"\n");
  ext_h = fopen (scat ("%s/extern.h", rsltpath), "w");

  strcpy (d_type, "real");
  writeF ("int calcFunc(void)\n{\n");
  writeF ("int err=0;\n");

  for (l = 1; l <= nmodelvar; l++)
    {
      if (modelvars[l].func)
        {
          ss =
            (char *) readExpression (modelvars[l].func, rd_cF, act_c, free);
          writeF ("   %s=%s;\n", modelvars[l].varname, ss + 3);
          free (ss);
          writeF ("   if(!finite(%s)) return FUCTION_ERROR;\n",
                  modelvars[l].varname);
        }
    }

  writeF ("if(err) return 1; else return 0;\n}\n");

  fclose (ext_h);
}


static void
make_pinf_F (void)
{
  int i;
  int ntot = getntot ();

  writeF ("int pinf_(int nsub,int nprtcl,char* pname, double* pmass)\n{\n");
  writeF ("int n;\n");
  writeF (" char names[%d][%d][7] =\n{", subproc_sq, ntot);

  inftmp = inf;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      for (i = 1; i <= ntot; i++)
        {
          writeF ("\"%s\"", inftmp->p_name[i - 1]);
          if (i == ntot)
            if (nsub == subproc_sq)
              writeF ("\n};\n");
            else
              writeF ("\n,");
          else
            writeF (",");
        }
      inftmp = inftmp->next;
    }
  writeF ("int nvalue[%d][%d]={\n", subproc_sq, ntot);

  inftmp = inf;
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      for (i = 1; i <= ntot; i++)
        {
          int k = inftmp->p_masspos[i - 1];

          /*if(k) sscanf(vararr[k].alias,"va[%d]",&k); */
          if (k > nvars)
            k++;
          writeF ("%d", k);
          if (i != ntot || nsub != subproc_sq)
            writeF (",");
          if (i == ntot)
            writeF ("\n");
        }
      if (nsub == subproc_sq)
        writeF ("};\n");
      inftmp = inftmp->next;
    }

  writeF ("if  (nsub<0 ||nsub>%d||nprtcl<0||nprtcl>%d) return 1;\n",
          subproc_sq, ntot);
  writeF ("if(pname) strcpy(pname,names[nsub-1][nprtcl-1]);\n");

  writeF ("if(pmass)\n{\n");
  writeF ("  n=nvalue[nsub-1][nprtcl-1];\n");

  writeF ("if (n>%d) if (calcFunc()) return FUCTION_ERROR;\n", nvars);
  writeF ("if (n==0) *pmass=0; else *pmass=va[n];\n");
  writeF ("if (*pmass<0) (*pmass)=-(*pmass);\n");
  writeF ("}\n");
  writeF ("}\n");
}



static void
make_vinf_F ()
{
  writeF ("int vinf_(int numvar,char * name, double * val)\n{\n");

  writeF ("char names[%d][10]={\"Sqrt(S)\"", nvars + nfunc + 1);
  scanvars_F (3);
  writeF ("};\n");
  writeF ("   if(numvar<0||numvar>%d  ) return 1;\n", nvars + nfunc);
  writeF ("   if(name) strcpy(name,names[numvar]);\n");
  writeF ("   if(val) *val=va[numvar];\n");

  writeF ("   return 0;\n}\n\n");
}


static void
make_fosimple_F (int *info_for_sqme)
{
  unsigned i, totcount;

  writeF ("static real smpl(int nsub, double * momenta,int * err)\n{\n");
  writeF (" real ans=0;\n");

  inftmp = inf;
  totcount = 0;

  writeF ("switch(nsub)\n{\n");
  for (nsub = 1; nsub <= subproc_sq; nsub++)
    {
      if (info_for_sqme[nsub] != 0)
        {
          writeF ("case %d: sprod_(momenta); ans= ", nsub);
          for (i = 1; i <= info_for_sqme[nsub]; i++)
            writeF ("+F%d()", totcount + i);
          totcount += info_for_sqme[nsub];
          writeF (";\n     break;\n");
        }
      /*inftmp = inftmp->next; */
    }
  writeF ("}\nreturn ans;\n}\n");
}


int
service_sqme (int NumberDiagrams, int *info_for_sqme)
{
  int breaker;
  int i;

  outputLanguage = 'c';
  diagrq = fopen (DIAGRQ_NAME, "rb");
  memerror = zeroHeep;
  mark_ (&heapbeg);
  sprintf (rsltpath, "%sresults", pathtouser);

  firstVar = nmodelvar;
  if (!strcmp (modelvars[firstVar].varname, strongconst))
    firstVar--;
  prepareprocinform ();
  calc_nvars_nfunc_F ();


  outFileOpen (scat ("%s/service.c", rsltpath));
  labl ();
  writeF ("#include<math.h>\n");
  writeF ("#define real double\n");
  writeF ("#include\"out_int.h\"\n");
  writeF ("#include\"out_ext.h\"\n");
  writeF ("#include\"var_def.h\"\n");
  writeF ("int gwidth=0;\n");
  writeF ("int rwidth=0;\n");
  common_F ("");
  writeF ("const int nin_ = %d;\n\n", getnin ());
  writeF ("const int nout_ = %d;\n\n", getnout ());
  writeF ("const int nprc_ = %d;\n\n", subproc_sq);
  make_pinf_F ();
  writeF ("const int nvar_ = %d;\n\n", nvars);
  writeF ("const int nfunc_ = %d;\n\n", nfunc);
  make_vinf_F ();
  make_asgn ();
  make_infbasis ();
  writesubroutineinit_F ();
  outFileClose ();


  outFileOpen (scat ("%s/sqme.c", rsltpath));
  labl ();
  writeF ("#include<math.h>\n");
  writeF ("#define real double\n");
  common_F ("extern");
  writeF ("#include\"out_int.h\"\n");
  writeF ("#include\"out_ext.h\"\n");
  writeF ("char  processch[] = \"%s\";\n", getProcessch ());
  writeF ("real DP[%d];\n", ((MAXINOUT * (MAXINOUT - 1) / 2)));

  ndiagrtot = NumberDiagrams;
  writeF ("extern FNN ");
  for (i = 1; i < ndiagrtot; i++)
    writeF ("F%d,", i);
  writeF ("F%d;\n", ndiagrtot);
  writeF ("static void sprod_(real*);\n");
  make_fosimple_F (info_for_sqme);

  writeF ("#include\"sqme0.c\"\n");
  outFileClose ();

  diagrcount = 0;
  inftmp = inf;
  init_stat ();

  clearstatistic ();
  fclose (diagrq);
  release_ (&heapbeg);

  return !breaker;
}
