/*
 * Copyright (C) 2001-2009, CompHEP Collaboration
 * Copyright (C) 1997, Alexander Pukhov
 * ------------------------------------------------------
 */
#include "chep_crt/include/crt_util.h"
#include "chep_crt/include/edittab.h"
#include "chep_crt/include/file_scr.h"
#include "out_ext.h"
#include "plot/include/plot.h"
#include "service2/include/chep_limits.h"
#include "service2/include/read_func.h"
#include "service2/include/syst.h"
#include "service2/include/tptcmac.h"

#include "histogram.h"
#include "phys_val.h"
#include "rd_num.h"
#include "roothisto.h"
#include "subproc.h"

table histTab = { "*** Table ***", "Distributions",
    "  Parameter  |> Min bound <|> Max bound <|> Rest Frame <|", NULL };

typedef struct histRec {
    struct histRec* next;
    midstr mother;
    long nPoints;
    char key;
    char plist[10];
    char restfr[10];
    double hMin;
    double hMax;
    double f[300];
    double ff[300];
} histRec;

static histRec* histPtr = NULL;

void clearHists(void)
{
    histRec* hptr = histPtr;
    for (; hptr != NULL; hptr = hptr->next) {
        int j;
        hptr->nPoints = 0;
        for (j = 0; j < 300; ++j) {
            hptr->f[j] = 0;
            hptr->ff[j] = 0;
        }
    }
}

void fillHists(double w)
{
    histRec* hists = histPtr;
    int i;
    double z;
    while (hists) {
        z = calcPhysVal(hists->key, hists->plist, hists->restfr);

        if (z > hists->hMin && z < hists->hMax) {
            i = 300 * (z - hists->hMin) / (hists->hMax - hists->hMin);
            hists->f[i] += w;
            hists->ff[i] += w * w;
        }
        hists->nPoints++;
        hists = hists->next;
    }
}

int correctHistList(int mode)
{
    int i, n, j;
    int lineNum = 0;
    int num_tab_list = 0;
    int* histlistentry_used;
    char keyChar;
    double rMin;
    double rMax;
    char fieldName[50];

    midstr histStr;
    midstr minStr;
    midstr maxStr;
    midstr restStr;

    linelist ln = histTab.strings;
    histRec* hptr = histPtr;

    while (ln != NULL) {
        ++num_tab_list;
        ln = ln->next;
    }
    while (hptr != NULL) {
        ++num_tab_list;
        hptr = hptr->next;
    }
    histlistentry_used = malloc(num_tab_list * sizeof(int));
    for (i = 0; i < num_tab_list; ++i)
        histlistentry_used[i] = 0;

    ln = histTab.strings;
    while (ln != NULL) {
        int add_new_histo = 1;
        char lv[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        char restnumbers[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        if (4 != sscanf(ln->line, "%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]%*c", histStr, minStr, maxStr, restStr))
            goto errorExit;
        lineNum++;

        trim(minStr);
        trim(histStr);
        trim(maxStr);
        trim(restStr);

        strcpy(fieldName, "Wrong field 'Min. bound'");
        if (!minStr[0] || calcExpression(minStr, rd_num, &rMin))
            goto errorExit;

        strcpy(fieldName, "Wrong field 'Max bound'");
        if (!maxStr[0] || calcExpression(maxStr, rd_num, &rMax))
            goto errorExit;

        strcpy(fieldName, "Wrong field 'Parameter'");
        if (!checkPhysVal(histStr, &keyChar, lv))
            goto errorExit;

        if (rMax < rMin) {
            double tmp = rMax;
            rMax = rMin;
            rMin = tmp;
        }

        strcpy(fieldName, "Wrong field 'Rest frame'");
        if (isalnum(restStr[0])) {
            i = 0;
            while (restStr[i] && restStr[i] != ' ') {
                n = restStr[i] - '0';
                if (n <= 0 || n > nin_ + nout_)
                    goto errorExit;
                for (j = 0; j < i; ++j)
                    if (restnumbers[j] == n)
                        goto errorExit;
                restnumbers[i] = n;
                i++;
            }

            while (i < 10) {
                if (restStr[i] && (restStr[i] != ' '))
                    goto errorExit;
                i++;
            }
        }

        for (i = 0, hptr = histPtr; hptr != NULL; ++i, hptr = hptr->next) {
            if (!histlistentry_used[i])
                if (!strcmp(hptr->mother, ln->line)) {
                    histlistentry_used[i] = 1;
                    add_new_histo = 0;
                }
        }
        if (add_new_histo) {
            hptr = malloc(sizeof(histRec));
            hptr->next = histPtr;
            strcpy(hptr->mother, ln->line);
            strcpy(hptr->restfr, restnumbers);
            hptr->key = keyChar;
            strcpy(hptr->plist, lv);
            hptr->hMin = rMin;
            hptr->hMax = rMax;
            for (j = 0; j < 300; ++j) {
                hptr->f[j] = 0;
                hptr->ff[j] = 0;
            }
            hptr->nPoints = 0;
            histPtr = hptr;
            histlistentry_used[i] = 1;
        }

        ln = ln->next;
    }

    hptr = histPtr;
    histPtr = NULL;
    for (i = 0; hptr != NULL; ++i) {
        histRec* hptr_ = hptr;
        hptr = hptr->next;
        if (!histlistentry_used[i])
            free(hptr_);
        else {
            hptr_->next = histPtr;
            histPtr = hptr_;
        }
    }

    return 0;

errorExit:
    sprintf(errorText, " Error in line %d .\nWrong value %s.", lineNum,
        fieldName);
    warnanykey(2, 10, errorText);
    if (mode)
        return 1;
    return 0;
}

static void
writeDistributions(FILE* iprt)
{
    histRec* hists = histPtr;
    int i;
    linelist rec;

    for (rec = histTab.strings; rec; rec = rec->next) {
        for (hists = histPtr; hists; hists = hists->next)
            if (!strcmp(hists->mother, rec->line)) {
                fprintf(iprt, " nPoints=%ld key=%c", hists->nPoints, hists->key);
                if ('U' != hists->key)
                    for (i = 0; i < 10 && hists->plist[i]; i++)
                        fprintf(iprt, "%d", hists->plist[i]);
                else
                    fprintf(iprt, "%s", hists->plist);
                fprintf(iprt, " hMin=%-12E hMax=%-12E", hists->hMin,
                    hists->hMax);
                fprintf(iprt, "\n  f: ");
                for (i = 0; i < 300; i++)
                    fprintf(iprt, " %-12E", hists->f[i]);
                fprintf(iprt, "\n ff: ");
                for (i = 0; i < 300; i++)
                    fprintf(iprt, " %-12E", hists->ff[i]);
            }
        fprintf(iprt, "\n");
    }
}

static void
readDistributions(FILE* iprt)
{
    linelist rec;
    for (rec = histTab.strings; rec; rec = rec->next) {
        int i;
        long tmp_nPoints;
        char tmp_key;
        double tmp_hMin;
        double tmp_hMax;
        char tmp_plist[10];
        double tmp_f[300];
        double tmp_ff[300];

        if (2 != fscanf(iprt, " nPoints=%ld key=%c", &tmp_nPoints, &tmp_key))
            return;
        for (i = 0; i < 10; i++) {
            fscanf(iprt, "%c", tmp_plist + i);
            if (tmp_plist[i] == ' ') {
                tmp_plist[i] = 0;
                break;
            } else if ('U' != tmp_key)
                tmp_plist[i] -= '0';
        }
        if (2 != fscanf(iprt, "hMin=%lf hMax=%lf\n", &tmp_hMin, &tmp_hMax))
            return;
        fscanf(iprt, "  f: ");
        for (i = 0; i < 300; i++)
            fscanf(iprt, " %lf", tmp_f + i);
        fscanf(iprt, " ff: ");
        for (i = 0; i < 300; i++)
            fscanf(iprt, " %lf", tmp_ff + i);
        {
            histRec* hist = malloc(sizeof(histRec));
            hist->next = histPtr;
            histPtr = hist;
            strcpy(hist->mother, rec->line);
            hist->nPoints = tmp_nPoints;
            hist->key = tmp_key;
            hist->hMin = tmp_hMin;
            hist->hMax = tmp_hMax;
            for (i = 0; i < 10; ++i)
                hist->plist[i] = tmp_plist[i];
            for (i = 0; i < 300; ++i) {
                hist->f[i] = tmp_f[i];
                hist->ff[i] = tmp_ff[i];
            }
        }
    }
}

int WriteHistograms(FILE* nchan)
{
    fprintf(nchan, "\n");
    writetable0(&histTab, nchan);
    writeDistributions(nchan);
    return 0;
}

int ReadHistograms(FILE* nchan)
{
    fscanf(nchan, "\n");
    readtable0(&histTab, nchan);
    readDistributions(nchan);
    correctHistList(0);
    return 0;
}

static int
nBinMenu(void)
{

    static int kBinMenu = 3;
    char strmen[] = "\015"
                    " 300         "
                    " 150         "
                    " 100         "
                    "  75         "
                    "  60         "
                    "  50         "
                    "  30         "
                    "  25         "
                    "  20         "
                    "  15         "
                    "  12         "
                    "  10         "
                    "  6          "
                    "  5          "
                    "  4          "
                    "  3          "
                    "  2          ";

    void* pscr = NULL;

    int n;
    if (!kBinMenu)
        kBinMenu = 3;
    for (;;) {
        menu1(54, 14, "number of bins", strmen, "", &pscr, &kBinMenu);
        if (kBinMenu) {
            sscanf(strmen + 1 + strmen[0] * (kBinMenu - 1), "%d", &n);
            put_text(&pscr);
            return n;
        }
        return 0;
    }
}

void showHistAll(void)
{
    char histStr[STRSIZ], minStr[STRSIZ], maxStr[STRSIZ], restStr[STRSIZ];
    ;
    linelist ln = histTab.strings;
    char* menutxt;
    void* pscr = NULL;
    int mode = 0;

    int npos = 0;
    int ii;
    int width, width2;
    char restnumbers[10] = { ' ' };

    while (ln) {
        npos++;
        ln = ln->next;
    }
    if (!npos)
        return;

    sscanf(histTab.format, "%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]%*c", histStr,
        minStr, maxStr, restStr);

    width = strlen(histStr);
    width2 = strlen(restStr);

    menutxt = malloc(2 + (width + width2 + 3) * npos);
    menutxt[0] = width + width2 + 3;
    menutxt[1] = 0;

    ln = histTab.strings;
    while (ln) {
        sscanf(ln->line, "%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]%*c", histStr, minStr,
            maxStr, restStr);
        strcat(menutxt, histStr);
        if ((restStr[0] == '|') || (restStr[0] == ' ') || (restStr[0] == '\n'))
            strcat(menutxt, "   ");
        else
            strcat(menutxt, "rf ");
        strcat(menutxt, restStr);
        ln = ln->next;
    }

    /*  for (;;) */
    {
        /*      menu1 (54, 10, "", menutxt, "", &pscr, &mode); */
        mode = 1;
        switch (mode) {
        case 0:
            free(menutxt);
            return;
        default: {
            histRec* hist = histPtr;
            int nBin;
            ln = histTab.strings;
            for (npos = 1; npos < mode; npos++)
                ln = ln->next;
            for (; hist && strcmp(hist->mother, ln->line); hist = hist->next) {
                ;
            }
            if (hist) {
                char xname[200], yname[80], units[80];
                if (hist->nPoints == 0)
                    messanykey(10, 10, "Distibution is empty");
                else
                    /*		  while ((nBin = nBinMenu ())) */
                    nBin = 35;
                {
                    double f[300], df[300], coeff;
                    int i;
                    coeff = nBin / (hist->nPoints * (hist->hMax - hist->hMin));
                    for (i = 0; i < nBin; i++) {
                        int k;
                        f[i] = 0;
                        df[i] = 0;
                        for (k = 0; k < 300 / nBin; k++) {
                            f[i] += hist->f[i * 300 / nBin + k];
                            df[i] += hist->ff[i * 300 / nBin + k];
                        }
                        f[i] *= coeff;
                        df[i] = coeff * sqrt(df[i] - f[i] * f[i] / hist->nPoints);
                    }

                    if (nin_ == 2)
                        strcpy(yname, "Diff. cross section [pb");
                    else
                        strcpy(yname, "Diff. width [GeV");
                    xName(hist->key, hist->plist, xname, units);
                    if (hist->restfr) {
                        ii = 0;
                        while (((hist->restfr[ii]) != 0) && (ii < 10)) {
                            restnumbers[ii] = hist->restfr[ii] + '0';
                            ii++;
                        }
                        if (ii != 0)
                            strcat(xname, ", rest frame ");
                        while (ii < 10) {
                            restnumbers[ii] = ' ';
                            ii++;
                        }
                        trim(restnumbers);
                        strcat(xname, restnumbers);
                    }
                    if (units[0]) {
                        strcat(yname, "/");
                        strcat(yname, units);
                    }
                    strcat(yname, "]");

                    plot_histo(hist->hMin, hist->hMax, nBin, f, df,
                        proces_1.proces, xname, yname);
                }
            }
        }
        }
    }
}

static void
showHist(void)
{
    char histStr[STRSIZ], minStr[STRSIZ], maxStr[STRSIZ], restStr[STRSIZ];
    ;
    linelist ln = histTab.strings;
    char* menutxt;
    void* pscr = NULL;
    int mode = 0;

    int npos = 0;
    int ii;
    int width, width2;
    char restnumbers[10] = { ' ' };

    while (ln) {
        npos++;
        ln = ln->next;
    }
    if (!npos)
        return;

    sscanf(histTab.format, "%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]%*c", histStr,
        minStr, maxStr, restStr);

    width = strlen(histStr);
    width2 = strlen(restStr);

    menutxt = malloc(2 + (width + width2 + 3) * npos);
    menutxt[0] = width + width2 + 3;
    menutxt[1] = 0;

    ln = histTab.strings;
    while (ln) {
        sscanf(ln->line, "%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]%*c", histStr, minStr,
            maxStr, restStr);
        strcat(menutxt, histStr);
        if ((restStr[0] == '|') || (restStr[0] == ' ') || (restStr[0] == '\n'))
            strcat(menutxt, "   ");
        else
            strcat(menutxt, "rf ");
        strcat(menutxt, restStr);
        ln = ln->next;
    }

    for (;;) {
        menu1(54, 10, "", menutxt, "", &pscr, &mode);
        switch (mode) {
        case 0:
            free(menutxt);
            return;
        default: {
            histRec* hist = histPtr;
            int nBin;
            ln = histTab.strings;
            for (npos = 1; npos < mode; npos++)
                ln = ln->next;
            for (; hist && strcmp(hist->mother, ln->line); hist = hist->next) {
                ;
            }
            if (hist) {
                char xname[200], yname[80], units[80];
                if (hist->nPoints == 0)
                    messanykey(10, 10, "Distibution is empty");
                else
                    while ((nBin = nBinMenu())) {
                        double f[300], df[300], coeff;
                        int i;
                        coeff = nBin / (hist->nPoints * (hist->hMax - hist->hMin));
                        for (i = 0; i < nBin; i++) {
                            int k;
                            f[i] = 0;
                            df[i] = 0;
                            for (k = 0; k < 300 / nBin; k++) {
                                f[i] += hist->f[i * 300 / nBin + k];
                                df[i] += hist->ff[i * 300 / nBin + k];
                            }
                            f[i] *= coeff;
                            df[i] = coeff * sqrt(df[i] - f[i] * f[i] / hist->nPoints);
                        }

                        if (nin_ == 2)
                            strcpy(yname, "Diff. cross section [pb");
                        else
                            strcpy(yname, "Diff. width [GeV");
                        xName(hist->key, hist->plist, xname, units);
                        if (hist->restfr) {
                            ii = 0;
                            while (((hist->restfr[ii]) != 0) && (ii < 10)) {
                                restnumbers[ii] = hist->restfr[ii] + '0';
                                ii++;
                            }
                            if (ii != 0)
                                strcat(xname, ", rest frame ");
                            while (ii < 10) {
                                restnumbers[ii] = ' ';
                                ii++;
                            }
                            trim(restnumbers);
                            strcat(xname, restnumbers);
                        }
                        if (units[0]) {
                            strcat(yname, "/");
                            strcat(yname, units);
                        }
                        strcat(yname, "]");

                        plot_histo(hist->hMin, hist->hMax, nBin, f, df,
                            proces_1.proces, xname, yname);
                    }
            }
        }
        }
    }
}

static void
editHist(void)
{
    do
        edittable(1, 4, &histTab, 1, "n_distrib", 0);
    while (correctHistList(1));
}

void manipulateHists(void)
{
    int mode;
    void* pscr = NULL;

    for (;;) {
        char strmen[] = "\034"
                        " Define variables           "
                        " Display 1d distributions   "
                        " Combine 1d distributions   "
                        " Sum 1d distributions       "
                        "----------------------------"
                        " Sum 1d tables              "
                        " Plot 1d table              "
                        " Construct 2d distribution  "
                        "----------------------------"
                        " Sum 2d tables              "
                        " Add constant to table terms"
                        " Multiply 2d tables         "
                        " Multiply 2d table to koeff."
                        " Divide 2d tables           "
                        " Rescale x variable         "
                        " Rescale y variable         "
                        " Plot 3d surface            "
                        " Plot 2d contour            "
                        "----------------------------"
                        " Calculate hi^2             "
                        " Sum hi^2 tables            "
                        " Plot delta hi^2 contour    "
                        " Erase data                 ";

        menu1(50, 7, "", strmen, "n_veg_*", &pscr, &mode);
        switch (mode) {
        case 0:
            return;
        case 1:
            editHist();
            break;
        case 2:
            showHist();
            break;
        case 3:
            combine_root();
            break;
        case 4:
            combine_rootSum();
            break;
        case 5: {
        } break;
        case 6:
            combine_Sum1d();
            break;
        case 7:
            plot_root1d();
            break;
        case 8:
            constract_2d();
            break;
        case 9: {
        } break;
        case 10:
            combine_Sum3d();
            break;
        case 11:
            sum_HistNum();
            break;
        case 12:
            multiply_Sum3d();
            break;
        case 13:
            multiply_HistNum();
            break;
        case 14:
            divide_Sum3d();
            break;
        case 15:
            rescale_xy(1);
            break;
        case 16:
            rescale_xy(2);
            break;
        case 17:
            plot_root3d();
            break;
        case 18:
            plot_rootCountur();
            break;
        case 19: {
        } break;
        case 20:
            calc_hi2();
            break;
        case 21:
            combine_root();
            break;
        case 22:
            plot_deltaHi2();
            break;
        case 23:
            clearHists();
            correctHistList(0);
            messanykey(10, 10, "Data in distributions have been erased");
            break;
        }
    }
}
