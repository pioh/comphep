/*
 * Copyright (C) 2001-2009, CompHEP Collaboration
 * Copyright (C) 1997, Alexander Pukhov
 * ------------------------------------------------------
 */
#include <math.h>
#include <string.h>

#include "chep_crt/include/chep_crt.h"
#include "service2/include/chep_limits.h"
#include "service2/include/files.h"
#include "service2/include/getmem.h"
#include "service2/include/syst.h"
#include "service2/include/unix_utils.h"

#include "diagrams.h"
#include "cweight.h"
#include "model.h"
#include "physics.h"
#include "prepdiag.h"
#include "process.h"
#include "process_core.h"
#include "sos.h"

static void
init_safe(int* exitlevel)
{
    int i;

    strcpy(beam[0].h.name, "p");
    strcpy(beam[1].h.name, "p");
    beam[0].energy = 7000.0;
    beam[1].energy = 7000.0;

    setFinalstatech("m,M,G");

    for (i = 0; i < MAXINOUT; i++) {
        strcpy(hadrons[i].name, "");
        strcpy(hadrons[i].sf_name, "");
        hadrons[i].sf_set = 0;
        hadrons[i].sf_mem = 0;
        hadrons[i].how = 0;
    }
    setExclprtlist("");
    setKeepprtlist("");
    setModelNumberSymb(4);
    *exitlevel = 0;
}

void restoreent(int* exitlevel)
{
    int i;
    int num;
    int d[32];
    shortstr tmps;
    shortstr name;
    int set, mem;
    double e[2], p[2];
    FILE* ff = fopen(scat("%stmp%csafe", pathtouser, f_slash), "r");

    if (ff != NULL) {
        d[0] = fscanf(ff, "#Model %d\n", &num);
        setModelNumberSymb(num);
        d[1] = fscanf(ff, "#nIn %d\n", &num);
        setnin(num);
        d[2] = fscanf(ff, "#nOut %d\n", &num);
        setnout(num);
        d[3] = fscanf(ff, "#beam_1 %[^\n]\n", beam[0].h.name);
        d[4] = fscanf(ff, "#beam_2 %[^\n]\n", beam[1].h.name);
        d[5] = fscanf(ff, "#beam_energy_1 %lf\n", &(beam[0].energy));
        d[6] = fscanf(ff, "#beam_energy_2 %lf\n", &(beam[1].energy));
        for (i = 0; i < 2; ++i) {
            double m = beam[i].h.mass;
            e[i] = beam[i].energy;
            p[i] = sqrt(e[i] * e[i] - m * m);
        }
        if (2 == getnin()) {
            set_sqrts_rap(e[0], p[0], e[1], p[1]);
        } else {
            setsqrtS(e[0]);
            setRapidity(0.0);
        }
        d[7] = fscanf(ff, "#strfun_1 %[^\n]\n", tmps);
        if (strstr(tmps, "LHA") || strstr(tmps, "PDF")) {
            if (get_sf_info(tmps, name, &set, &mem)) {
                strcpy(beam[0].h.sf_name, name);
                beam[0].h.sf_set = set;
                beam[0].h.sf_mem = mem;
            }
        } else {
            strcpy(beam[0].h.sf_name, tmps);
            beam[0].h.sf_set = 0;
            beam[0].h.sf_mem = 0;
        }
        d[8] = fscanf(ff, "#strfun_2 %[^\n]\n", tmps);
        if (strstr(tmps, "LHA") || strstr(tmps, "PDF")) {
            if (get_sf_info(tmps, name, &set, &mem)) {
                strcpy(beam[1].h.sf_name, name);
                beam[1].h.sf_set = set;
                beam[1].h.sf_mem = mem;
            }
        } else {
            strcpy(beam[1].h.sf_name, tmps);
            beam[1].h.sf_set = 0;
            beam[1].h.sf_mem = 0;
        }
        d[9] = fscanf(ff, "#Final_state %[^\n]\n", tmps);
        trim(tmps);
        setFinalstatech(tmps);
        d[10] = fscanf(ff, "#Remove_Virtual%[^\n]\n", tmps);
        trim(tmps);
        setExclprtlist(tmps);
        d[11] = fscanf(ff, "#Keep_Virtual%[^\n]\n", tmps);
        trim(tmps);
        setKeepprtlist(tmps);
        d[12] = fscanf(ff, "#nSubproc(ampl) %d\n", &subproc_f);
        d[13] = fscanf(ff, "#nSubproc(squared) %d\n", &subproc_sq);
        d[14] = fscanf(ff, "#ConservationLaw %d\n", &consLow);
        d[15] = fscanf(ff, "#Nc==inf %d\n", &num);
        setNcinflimit(num);
        d[16] = fscanf(ff, "#ExitCode %d\n", exitlevel);

        if (1 != d[0] || 1 != d[1] || 1 != d[2] || 1 != d[3] || 1 != d[4] || 
            1 != d[5] || 1 != d[6] || 1 != d[7] || 1 != d[8] || 1 != d[9] || 
            1 != d[10] || 1 != d[11] || 1 != d[12] || 1 != d[13] || 1 != d[14] ||
            1 != d[15] || 1 != d[16])
            init_safe(exitlevel);

        fclose(ff);
        if (1 == getnin()) {
            shortstr tmps;
            sprintf(tmps, "%s -> %s", beam[0].h.name, getFinalstatech());
            setProcessch(tmps);
        } else {
            shortstr tmps;
            sprintf(tmps, "%s,%s -> %s", beam[0].h.name, beam[1].h.name, getFinalstatech());
            setProcessch(tmps);
        }
    } else
        init_safe(exitlevel);
}

int restoreent_dump(void)
{
    int subproc = 0;
    int num;
    double num1;
    int d[32];
    shortstr tmps;
    FILE* ff = fopen(scat("%stmp%csafe", pathtouser, f_slash), "r");

    if (ff != NULL) {
        d[0] = fscanf(ff, "#Model %d\n", &num);
        d[1] = fscanf(ff, "#nIn %d\n", &num);
        d[2] = fscanf(ff, "#nOut %d\n", &num);
        d[3] = fscanf(ff, "#beam_1 %[^\n]\n", tmps);
        d[4] = fscanf(ff, "#beam_2 %[^\n]\n", tmps);
        d[5] = fscanf(ff, "#beam_energy_1 %lf\n", &num1);
        d[6] = fscanf(ff, "#beam_energy_2 %lf\n", &num1);
        d[7] = fscanf(ff, "#strfun_1 %[^\n]\n", tmps);
        d[8] = fscanf(ff, "#strfun_2 %[^\n]\n", tmps);
        d[9] = fscanf(ff, "#Final_state %[^\n]\n", tmps);
        d[10] = fscanf(ff, "#Remove_Virtual%[^\n]\n", tmps);
        d[11] = fscanf(ff, "#Keep_Virtual%[^\n]\n", tmps);
        d[12] = fscanf(ff, "#nSubproc(ampl) %d\n", &num);
        d[13] = fscanf(ff, "#nSubproc(squared) %d\n", &subproc);
        d[14] = fscanf(ff, "#ConservationLaw %d\n", &num);
        d[15] = fscanf(ff, "#Nc==inf %d\n", &num);
        d[16] = fscanf(ff, "#ExitCode %d\n", &num);
    }
    return subproc;
}

void saveent(int exitlevel)
{
    shortstr sf1;
    shortstr sf2;
    FILE* ff = fopen(scat("%stmp%csafe", pathtouser, f_slash), "w");

    if (strstr(beam[0].h.sf_name, "LHA") || strstr(beam[0].h.sf_name, "PDF")) {
        sprintf(sf1, "%s:%i:%i", beam[0].h.sf_name, beam[0].h.sf_set, beam[0].h.sf_mem);
    } else {
        sprintf(sf1, "%s", beam[0].h.sf_name);
    }
    if (strstr(beam[1].h.sf_name, "LHA") || strstr(beam[1].h.sf_name, "PDF")) {
        sprintf(sf2, "%s:%i:%i", beam[1].h.sf_name, beam[1].h.sf_set, beam[1].h.sf_mem);
    } else {
        sprintf(sf2, "%s", beam[1].h.sf_name);
    }

    fprintf(ff, "#Model %i\n", getModelNumberSymb());
    fprintf(ff, "#nIn %i\n", getnin());
    fprintf(ff, "#nOut %i\n", getnout());
    fprintf(ff, "#beam_1 %s\n", beam[0].h.name);
    fprintf(ff, "#beam_2 %s\n", beam[1].h.name);
    fprintf(ff, "#beam_energy_1 %f\n", beam[0].energy);
    fprintf(ff, "#beam_energy_2 %f\n", beam[1].energy);
    fprintf(ff, "#strfun_1 %s\n", sf1);
    fprintf(ff, "#strfun_2 %s\n", sf2);
    fprintf(ff, "#Final_state %s\n", getFinalstatech());
    fprintf(ff, "#Remove_Virtual %s\n", getExclprtlist());
    fprintf(ff, "#Keep_Virtual %s\n", getKeepprtlist());
    fprintf(ff, "#nSubproc(ampl) %d\n", subproc_f);
    fprintf(ff, "#nSubproc(squared) %d\n", subproc_sq);
    fprintf(ff, "#ConservationLaw %d\n", consLow);
    fprintf(ff, "#Nc==inf %d\n", getNcinflimit());
    fprintf(ff, "#ExitCode %d\n", exitlevel);
    fclose(ff);
}

void save_sos(int ercode)
{
    unsigned nproc;
    csdiagram cd;
    marktp mark;

    if (ercode == -2) /* User Break */
    {
        saveent(7);
        finish("Restart of the CompHEP program.");
        exit(13); /*  Restart  */
    }

    if (ercode == -1) /* Heap is empty */
    {
        mark.blk_ = NULL;
        mark.pos_ = 0;
        release_(&mark);
    } /* Heap is empty, continue */

    /*  TooLargeNumber  */
    /*  TooManyIdentifiers  */
    /*  RangeCheckError  */
    if ((ercode < 0) || (ercode == 7) || (ercode == 11) || (ercode == 12)) /*  Restart  */
    {
        saveent(9);
        nproc = ftell(diagrq) - sizeof(cd);
        fseek(diagrq, nproc, SEEK_SET);
        FREAD1(cd, diagrq);
        cd.status = -2;
        fseek(diagrq, nproc, SEEK_SET);
        FWRITE1(cd, diagrq);
        finish("*** Error! CompHEP is restarted");
        exit(13); /*  Restart  */
    }

    /*  not disk space  */
    if ((ercode == 40)) {
        /*  freez */
        saveent(6);
        finish("*** Error!");
        exit(0); /*  End of work  */
    }
    if (ercode == 14) {
        saveent(1);
        warnanykey(10, 10, " Check model !");
        finish("*** Error!");
        exit(13);
    }
}

int batch_composite(char* name, shortstr h)
{
    int stop = 0;
    char nm[100], buff[500];
    /*char content[500]; */
    FILE* f = NULL;

    f = fopen(scat("%stmp%ccomposite.dat", pathtouser, f_slash), "r");
    if (f != NULL) {
        while (!stop) {
            if (fscanf(f, "%[^\n]", buff)) {
                int j = 0, i = -1;
                while (buff[i++] != ':')
                    nm[i] = buff[i];
                nm[i - 1] = '\0';
                if (!strcmp(name, nm)) {
                    while (buff[i + j] != '\0') {
                        nm[j] = buff[i + j];
                        j++;
                    }
                    strncpy(h, nm, j);
                    h[j] = 0;
                    fclose(f);
                    return 0;
                }
            } else {
                sprintf(buff, "composite.dat does non contain the particle %s", name);
                fclose(f);
                batch_error(buff, 1);
                stop = 1;
            }
        }
    }
    batch_error("File composite.dat have to be located in tmp", 1);
    return 1;
}
