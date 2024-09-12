/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef __PROCESS_CORE__
#define __PROCESS_CORE__

extern int getnin (void);
extern void setnin (int num);
extern int getnout (void);
extern void setnout (int num);
extern int getntot (void);
extern void setntot (int num);
extern int getnx (void);
extern void setnx (int num);
extern double getsqrtS (void);
extern void setsqrtS (double val);
extern double getRapidity (void);
extern void setRapidity (double val);

int getNcinflimit (void);
void setNcinflimit (int num);

extern void set_sqrts_rap (double e1, double p1, double e2, double p2);

extern char * getFinalstatech (void);
extern char * getProcessch (void);
extern char * getExclprtlist (void);
extern char * getKeepprtlist (void);
extern void setFinalstatech (shortstr s);
extern void setProcessch    (shortstr s);
extern void setExclprtlist  (shortstr s);
extern void setKeepprtlist  (shortstr s);
#endif
