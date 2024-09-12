/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------
*/
#ifndef __LESHOUCHES__
#define __LESHOUCHES__

#define MAXproc 1024
#define MAXpart 512

/**************************************************************/
/*********      Les Houches process c-structure      **********/
typedef struct processUP
  {
    int IDBeamUP[2];            /* ID of beam particles */
    int PDFgroupUP[2];          /* the author group for beams, according to the PDFlib specification */
    int PDFidUP[2];             /* the PDF set ID for beams, according to the PDFlib specification */
    int switchUP;               /* master switch dictating how the event weights (XWGTUP) are interpreted */
    int NprocRUP;               /* the number of different user subprocesses */
    int listprocRUP[MAXproc];   /* listing of all user process IDs that can appear */
    /* in IDPRUP of HEPEUP for this run */

    double energyUP[2];         /* energy in GeV of beam particles */
    double crossecUP[MAXproc];  /* the cross section for all processes in pb */
    double errorcsUP[MAXproc];  /* the cross section statistical error of processes in pb */
    double maxweightUP[MAXproc];        /*  maximum event weight for all processes */
    char PRnameUP[MAXproc][1000];       /* name of Process */
  }
processUP;

typedef struct process_map
  {
    int file;
    int number;
    char *filename;
    int used;
  }
process_map;

typedef struct process_
  {
    process_map *maps;
    processUP proc_info;
  }
process_;


/**************************************************************/
/*********      Les Houches event c-structure        **********/
typedef struct eventUP
  {
    int NpartUP;                /* number of particle entries in this event */
    int IDprocUP;               /*ID of the process for this event */
    int IDpartUP[MAXpart];      /*particle ID according to Particle Data Group convention */
    int statusUP[MAXpart];      /*status code:
                                   -1 Incoming particle
                                   +1 Outgoing final state particle
                                   -2 Intermediate space-like propagator defining an 
                                   $x$ and $Q^2$ which should be preserved
                                   +2 Intermediate resonance, Mass should be preserved
                                   +3 Intermediate resonance, for documentation only
                                   -9 Incoming beam particles at time $t=-\infty$ */

    int motherUP[2][MAXpart];   /*index of first and last mother */
    int colorflowUP[2][MAXpart];        /*
                                           ICOLUP(1,I) - integer tag for the color flow 
                                           line passing through the color of the particle
                                           ICOLUP(2,I) - integer tag for the color flow 
                                           line passing through the anti-color of the particle */

    double XweightUP;           /* event weight */
    double QscaleUP;            /* scale of the event in GeV, as used for calculation of  PDFs */
    double QEDalphaUP;          /* the QED coupling used for this event */
    double QCDalphaUP;          /* the QCD coupling used for this event */
    double momentumUP[5][MAXpart];      /* lab frame momentum $(P_x, P_y, P_z, E, M)$ of particle in GeV */
    double timelifeUP[MAXpart]; /* invariant lifetime $c\tau$ (distance from 
                                   production to decay) in mm */
    double spinUP[MAXpart];     /*cosine of the angle between the spin-vector of 
                                   particle I and the 3-momentum of the decaying 
                                   particle, specified in the lab frame */

  }
eventUP;

typedef struct eventIN
  {
    int NpartIN;
    int IDprocIN;
    int IDpartIN[MAXpart];
    int statusIN[MAXpart];
    int motherIN[2 * MAXpart];
    int colorflowIN[2 * MAXpart];
    int XweightIN;
    int QscaleIN;
    int QEDalphaIN;
    int QCDalphaIN;
    int momentumIN[5 * MAXpart];
    int timelifeIN[MAXpart];
    int spinUP[MAXpart];
  }
eventIN;

#endif /* LesHouches.h */
