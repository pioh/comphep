/*
 * Copyright (C) 2001-2009, CompHEP Collaboration
 * Copyright (C) 1997, Alexander Pukhov
 * --------------------------------------------------
 */
#include "chep_crt/include/chep_crt.h"
#include "out_ext.h"
#include "plot/include/plot.h"
#include "service2/include/chep_limits.h"
#include "service2/include/files.h"
#include "service2/include/syst.h"
#include "service2/include/tptcmac.h"
#include "service2/include/unix_utils.h"

#include "alphas_menu.h"
#include "core_data.h"
#include "cs_22.h"
#include "cut.h"
#include "err_code.h"
#include "evnt_menu.h"
#include "histogram.h"
#include "kininpt.h"
#include "mc_menu.h"
#include "num_serv.h"
#include "param.h"
#include "q_kin.h"
#include "regul.h"
#include "runVegas.h"
#include "rw_sess.h"
#include "strfun.h"
#include "strfun_par.h"
#include "subproc.h"
#include "vegas.h"

static void
write3dhisto(int nprm1, int nprm2, double xMin1, double xMax1, double xMin2, double xMax2, int dim1, int dim2, double* f2d)
{
    int i, j, totaldim;

    midstr outname;
    midstr outname_C;
    FILE* rf;
    double stepprm1, stepprm2, prmval1, prmval2;

    nextFileName(outname, "hist2d_", ".txt");
    strcpy(outname_C, outname);
    strcat(outname_C, ".txt");
    rf = fopen(outname_C, "w");
    if (NULL == rf)
        messanykey(10, 12, scat("Error! file\n%s can not be open", outname_C));

    totaldim = dim1 * dim2;

    stepprm1 = (xMax1 - xMin1) / (dim1 - 1);
    stepprm2 = (xMax2 - xMin2) / (dim2 - 1);

    prmval1 = xMin1;
    prmval2 = xMin2;

    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim2; j++) {
            f_printf(rf, "% .6E % .6E % .6E\n", prmval1, prmval2, f2d[i * dim2 + j]);
            prmval2 += stepprm2;
        }
        prmval2 = xMin2;
        prmval1 += stepprm1;
    }

    fclose(rf);
}

static void
write1dhisto(int nprm, double xMin, double xMax, int dim, double* f)
{
    int i, j;

    midstr outname;
    midstr outname_C;
    FILE* rf;
    double stepprm, prmval;

    nextFileName(outname, "hist1d_", ".txt");
    strcpy(outname_C, outname);
    strcat(outname_C, ".txt");
    rf = fopen(outname_C, "w");
    if (NULL == rf)
        messanykey(10, 12, scat("Error! file\n%s can not be open", outname_C));

    stepprm = (xMax - xMin) / (dim - 1);

    prmval = xMin;

    for (i = 0; i < dim; i++) {
        f_printf(rf, "% .6E % .6E\n", prmval, f[i]);
        prmval += stepprm;
    }

    fclose(rf);
}

void paramtable2(void)
{
    FILE* tablefile;
    int subprocnum = 1;
    void* pscr0 = NULL;
    char prmname[10];
    int nprm = 0, nprm1, nprm2;
    int parameternum = 1;
    int dummyflag1 = 0, dummyflag2 = 0;
    double minprm, maxprm, minprm1, maxprm1, minprm2, maxprm2;
    double memprm, memprm1, memprm2, stepprm, stepprm1, stepprm2;
    int npoints, npoints1, npoints2;
    unsigned count1 = 1, count2 = 1;

label5:

    if (parameternum == 1)
        messanykey(15, 15, "Enter first parameter");
    if (parameternum == 2)
        messanykey(15, 15, "Enter second parameter");

    selectParam(&nprm, 54, 11, nin_ == 2, 1, 0, "Choose parameter", NULL);

    if (nprm > nvar_) {
        memprm = 0;
        minprm = 0;
        maxprm = 0;
    } else {
        if (nprm < 0)
            return;
        vinf_(nprm, prmname, &memprm);
        minprm = memprm;
        maxprm = memprm;
    }

label1:

    if (!correctDouble(55, 14, "min= ", &minprm, 0)) {
        nprm = 0;
        goto exi;
    }

label2:
    if (!correctDouble(55, 15, "max= ", &maxprm, 0)) {
        goto_xy(55, 14);
        clr_eol();
        goto label1;
    }
    if (maxprm < minprm) {
        warnanykey(55, 17, "Range check error");
        goto_xy(55, 15);
        clr_eol();
        goto label2;
    }

label4:
    npoints = 100;
    if (correctInt(55, 16, "Number of points= ", &npoints, 0)) {
        /*      if (npoints < 3) */
        if (npoints < 1) {
            warnanykey(55, 17, "Too few points!");
            goto label4;
        }
        if (npoints > 100) {
            warnanykey(55, 17, "Too many points!");
            goto label4;
        }
    } else {
        goto_xy(55, 15);
        clr_eol();
        goto label2;
    }

    goto_xy(55, 14);
    clr_eol();
    goto_xy(55, 15);
    clr_eol();
    goto_xy(55, 16);
    clr_eol();

    informline(0, npoints);

    err_code = 0;

    if (parameternum == 1) {
        npoints1 = npoints;
        minprm1 = minprm;
        maxprm1 = maxprm;
        nprm1 = nprm;

        if (nprm1 > nvar_)
            dummyflag1 = 1;
        else
            dummyflag1 = 0;

        parameternum = 2;
        goto label5;
    }

    if (parameternum == 2) {
        npoints2 = npoints;
        minprm2 = minprm;
        maxprm2 = maxprm;
        nprm2 = nprm;

        if (nprm2 > nvar_)
            dummyflag2 = 1;
        else
            dummyflag2 = 0;
    }

    tablefile = fopen("tablecalc.dat", "w");

    fprintf(tablefile, "subprocnum= %3d\n", subprocnum);
    fprintf(tablefile, "count1= %3d\n", count1);
    fprintf(tablefile, "count2= %3d\n", count2);
    fprintf(tablefile, "dummyflag1= %d\n", dummyflag1);
    fprintf(tablefile, "param1num= %d\n", nprm1);
    fprintf(tablefile, "minprm1= %e\n", minprm1);
    fprintf(tablefile, "maxprm1= %e\n", maxprm1);
    fprintf(tablefile, "npoints1= %d\n", npoints1);

    fprintf(tablefile, "dummyflag2= %d\n", dummyflag2);
    fprintf(tablefile, "param2num= %d\n", nprm2);
    fprintf(tablefile, "minprm2= %e\n", minprm2);
    fprintf(tablefile, "maxprm2= %e\n", maxprm2);
    fprintf(tablefile, "npoints2= %d\n", npoints2);

    fclose(tablefile);

    tablefile = fopen("table_batch", "w");
    fprintf(tablefile, "#!/bin/sh\n");
    fprintf(tablefile, "./n_comphep -blind \"]]]]]]]]]}]]]]]]]]}{{9} &\"\n");
    fclose(tablefile);
    system("chmod 755 table_batch");

exi:
    put_text(&pscr0);
    goto_xy(1, 23);
    scrcolor(FGmain, BGmain);
    clr_eol();
}

/* void paramdependence (r_func ff, char *procname, char *resultname) */
void paramdependence(double (*ff)(void), char* procname, char* resultname)
{
    FILE* tablefile;
    long POS1, POS2, POS3;
    int subprocnum;
    void* pscr0 = NULL;
    char prmname[10];
    int nprm = 0, nprm1, nprm2;
    int parameternum = 1;
    int dummyflag1 = 0, dummyflag2 = 0;
    double minprm, maxprm, minprm1, maxprm1, minprm2, maxprm2;
    int npoints, npoints1, npoints2;
    double memprm, memprm1, memprm2, stepprm, stepprm1, stepprm2;
    double memprm1tot, memprm2tot;
    unsigned count, count1, count2;
    double f[200], f2d[10000];

    double prmval, prmval1, prmval2;

    tablefile = fopen("tablecalc.dat", "r+");
    if (tablefile == NULL) {
        fprintf(stderr, "\n file tablecalc.dat not found\n");
        goto exi;
    }

    POS3 = ftell(tablefile);
    fscanf(tablefile, "subprocnum= %d\n", &subprocnum);
    POS1 = ftell(tablefile);
    fscanf(tablefile, "count1= %d\n", &count1);
    POS2 = ftell(tablefile);
    fscanf(tablefile, "count2= %d\n", &count2);

    fscanf(tablefile, "dummyflag1= %d\n", &dummyflag1);
    fscanf(tablefile, "param1num= %d\n", &nprm1);
    fscanf(tablefile, "minprm1= %le\n", &minprm1);
    fscanf(tablefile, "maxprm1= %le\n", &maxprm1);
    fscanf(tablefile, "npoints1= %d\n", &npoints1);

    fscanf(tablefile, "dummyflag2= %d\n", &dummyflag2);
    fscanf(tablefile, "param2num= %d\n", &nprm2);
    fscanf(tablefile, "minprm2= %le\n", &minprm2);
    fscanf(tablefile, "maxprm2= %le\n", &maxprm2);
    fscanf(tablefile, "npoints2= %d\n", &npoints2);

    stepprm1 = (maxprm1 - minprm1) / (npoints1 - 1);
    prmval1 = minprm1;
    memprm1 = minprm1;
    vinf_(nprm1, prmname, &memprm1tot);

    stepprm2 = (maxprm2 - minprm2) / (npoints2 - 1);
    prmval2 = minprm2;
    memprm2 = minprm2;
    vinf_(nprm2, prmname, &memprm2tot);

    /*------------------------------------------------
    label5:

      if(parameternum==1) messanykey (15, 15,"Enter first parameter");
      if(parameternum==2) messanykey (15, 15,"Enter second parameter");

      selectParam (&nprm, 54, 11, nin_ == 2, 1, 0, "Choose parameter", NULL);

      if(nprm > nvar_)
       {
         memprm = 0;
         minprm = memprm;
         maxprm = minprm;
       }
      else
       { if (nprm < 0)
          return;
         vinf_ (nprm, prmname, &memprm);
         minprm = memprm;
         maxprm = minprm;
       }

    label1:

      if (!correctDouble (55, 14, "min= ", &minprm, 0))
        {
          nprm = 0;
          goto exi;
        }

    label2:
      if (!correctDouble (55, 15, "max= ", &maxprm, 0))
        {
          goto_xy (55, 14);
          clr_eol ();
          goto label1;
        }
      if (maxprm < minprm)
        {
          warnanykey (55, 17, "Range check error");
          goto_xy (55, 15);
          clr_eol ();
          goto label2;
        }

    label4:npoints = 100;
      if (correctInt (55, 16, "Number of points= ", &npoints, 0))
        {
          if (npoints < 3)
            {
              warnanykey (55, 17, "Too few points!");
              goto label4;
            }
          if (npoints > 100)
            {
              warnanykey (55, 17, "Too many points!");
              goto label4;
            }
        }
      else
        {
          goto_xy (55, 15);
          clr_eol ();
          goto label2;
        }

      goto_xy (55, 14);
      clr_eol ();
      goto_xy (55, 15);
      clr_eol ();
      goto_xy (55, 16);
      clr_eol ();


      informline (0, npoints);

      stepprm = (maxprm - minprm) / (npoints - 1);
      prmval = minprm;
      err_code = 0;

      if(parameternum==1)
       { npoints1=npoints;
         stepprm1=stepprm;
         minprm1=minprm;
         maxprm1=maxprm;
         nprm1=nprm;
         prmval1=prmval;
         memprm1=prmval1;

         if(nprm1 > nvar_) dummyflag1=1;
         else dummyflag1=0;

         parameternum=2;
         goto label5;
        }

      if(parameternum==2)
       { npoints2=npoints;
         stepprm2=stepprm;
         minprm2=minprm;
         maxprm2=maxprm;
         nprm2=nprm;
         prmval2=prmval;
         memprm2=prmval2;

         if(nprm2 > nvar_) dummyflag2=1;
         else dummyflag2=0;
        }


    ------------------------------------------------*/

    for (subprocnum = 1; subprocnum <= nprc_; subprocnum++) {
        fseek(tablefile, POS3, SEEK_SET);
        fprintf(tablefile, "subprocnum= %3d", subprocnum);

        proces_1.nsub = subprocnum;
        ComposeSubprocessString();
        if (nin_ == 2) {
            vshortstr name[2];
            pinf_(subprocnum, 1, name[0], NULL);
            pinf_(subprocnum, 2, name[1], NULL);
            initStrFun(name[0], name[1]);
        }
        /* fprintf(stderr,"\nprocnum= %d\n",subprocnum); */

        for (count1 = 1; count1 <= npoints1; count1++) {

            asgn_(nprm1, prmval1);

            /*   fprintf(stderr,"\n par1= %e\n",prmval1); */
            fseek(tablefile, POS1, SEEK_SET);
            fprintf(tablefile, "count1= %3d", count1);

            if (calcFunc())
                err_code = 5;
            else {
                for (count2 = 1; count2 <= npoints2; count2++) {
                    asgn_(nprm2, prmval2);

                    /*   fprintf(stderr," par2= %e ",prmval2); */
                    fseek(tablefile, POS2, SEEK_SET);
                    fprintf(tablefile, "count2= %3d", count2);

                    if (calcFunc())
                        err_code = 5;
                    else {
                        if (count1 > 1 && dummyflag1 == 1)
                            f2d[(count1 - 1) * npoints2 + (count2 - 1)] = f2d[(count1 - 2) * npoints2 + (count2 - 1)];
                        else {
                            if (count2 > 1 && dummyflag2 == 1)
                                f2d[(count1 - 1) * npoints2 + (count2 - 1)] = f2d[(count1 - 1) * npoints2 + (count2 - 2)];
                            else {
                                fillRegArray();
                                f2d[(count1 - 1) * npoints2 + (count2 - 1)] = (*ff)();
                            }
                        }

                        /*             showHistAll(); */
                    }
                    if (err_code) {
                        errormessage();
                        break;
                    } else {
                        informline(count2, npoints2);
                        prmval2 += stepprm2;
                    }
                }
                /*  fprintf(stderr,"\n"); */
                asgn_(nprm2, memprm2);
                prmval2 = memprm2;
            }

            if (err_code) {
                errormessage();
                break;
            } else {
                informline(count1, npoints1);
                prmval1 += stepprm1;
            }
        }
        asgn_(nprm1, memprm1);
        prmval1 = memprm1;

        if (err_code == 0) {
            scrcolor(FGmain, BGmain);
            goto_xy(1, 23);
            clr_eol();

            /*  fprintf(stderr,"\n"); */

            write3dhisto(nprm1, nprm2, minprm1, maxprm1, minprm2, maxprm2, npoints1, npoints2, f2d);
        }
    }

    asgn_(nprm1, memprm1tot);
    asgn_(nprm2, memprm2tot);

    calcFunc();

    /*  fclose(tablefile);*/

    /*
      clearHists ();
      init_vegas_integral ();
      ClearEventMax ();
      ClearVegasGrid ();
      set_nsession (get_nsession () + 1);
      write_session ();
    */

exi:
    put_text(&pscr0);
    goto_xy(1, 23);
    scrcolor(FGmain, BGmain);
    clr_eol();
}

void paramdependence1(double (*ff)(void), char* procname, char* resultname)
{
    int subprocnum;
    void* pscr0 = NULL;
    char prmname[10];
    int nprm = 0, nprm1, nprm2;
    int parameternum = 1;
    double minprm, maxprm, minprm1, maxprm1, minprm2, maxprm2;
    int npoints, npoints1, npoints2;
    double memprm, memprm1, memprm2, stepprm, stepprm1, stepprm2;
    unsigned count, count1, count2;
    double f[200], f2d[10000];

    double prmval, prmval1, prmval2;

label5:

    if (parameternum == 1)
        messanykey(15, 15, "Enter parameter");

    selectParam(&nprm, 54, 11, nin_ == 2, 1, 0, "Choose parameter", NULL);
    if (nprm < 0)
        return;
    vinf_(nprm, prmname, &memprm);
    minprm = memprm;
    maxprm = minprm;

label1:

    if (!correctDouble(55, 14, "min= ", &minprm, 0)) {
        nprm = 0;
        goto exi;
    }

label2:
    if (!correctDouble(55, 15, "max= ", &maxprm, 0)) {
        goto_xy(55, 14);
        clr_eol();
        goto label1;
    }
    if (maxprm <= minprm) {
        warnanykey(55, 17, "Range check error");
        goto_xy(55, 15);
        clr_eol();
        goto label2;
    }

label4:
    npoints = 200;
    if (correctInt(55, 16, "Number of points= ", &npoints, 0)) {
        if (npoints < 3) {
            warnanykey(55, 17, "Too few points!");
            goto label4;
        }
        if (npoints > 200) {
            warnanykey(55, 17, "Too many points!");
            goto label4;
        }
    } else {
        goto_xy(55, 15);
        clr_eol();
        goto label2;
    }

    goto_xy(55, 14);
    clr_eol();
    goto_xy(55, 15);
    clr_eol();
    goto_xy(55, 16);
    clr_eol();

    informline(0, npoints);

    stepprm = (maxprm - minprm) / (npoints - 1);
    prmval = minprm;
    err_code = 0;

    printf("\n npoints= %d", npoints);
    printf("\n memprm= %e", memprm);

    for (subprocnum = 1; subprocnum <= nprc_; subprocnum++) {
        proces_1.nsub = subprocnum;
        ComposeSubprocessString();
        if (nin_ == 2) {
            vshortstr name[2];
            pinf_(subprocnum, 1, name[0], NULL);
            pinf_(subprocnum, 2, name[1], NULL);
            initStrFun(name[0], name[1]);
        }
        printf("\n procnum= %d\n", subprocnum);

        for (count1 = 1; count1 <= npoints; count1++) {
            clearHists();
            correctHistList(0);

            asgn_(nprm, prmval);
            printf("\n par= %e\n", prmval);

            if (calcFunc())
                err_code = 5;
            else {
                fillRegArray();
                f[count1 - 1] = (*ff)();
            }

            if (err_code) {
                errormessage();
                break;
            } else {
                informline(count1, npoints);
                prmval += stepprm;
            }

            showHistAll();
        }
        asgn_(nprm, memprm);
        prmval = memprm;

        if (err_code == 0) {
            scrcolor(FGmain, BGmain);
            goto_xy(1, 23);
            clr_eol();

            /*   plot_table1d (minprm, maxprm, npoints, f, NULL, procname, prmname, resultname); */
            write1dhisto(nprm, minprm, maxprm, npoints, f);
        }
    }

    calcFunc();

exi:
    put_text(&pscr0);
    goto_xy(1, 23);
    scrcolor(FGmain, BGmain);
    clr_eol();
}
