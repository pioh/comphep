/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Author: V.Edneral
* ------------------------------------------------------------
*/
#ifndef __SETS__
#define __SETS__

/*
 *   By V.Edneral
 *
 *   int setof(int a,int b,...,_E)
 *      construct and return a set of the specified character values
 *
 *   int inset(int ex,int setrec)
 *      predicate returns true if expression ex is a member of
 *      the set parameter
 *   (It is suppoused enum type setrec takes two bytes here; a,b,ex < 16)
 *
 *   setofbyte setofb(int a,int b,...,_E)
 *   int insetb(unsigned ex, setofbyte setrec)
 *   Here a,b,ex<256; setofbyte=unsigned[16].
 *
 *   setofb_cpy( setofbyte dest, setofbyte source)
 *       for setofbyte copying
 *   setofb_uni( setofbyte a, setofbyte b)
 *       for unification
 *   setofb_aun( setofbyte a, setofbyte b)
 *       for antiunification
 *   setofb_its( setofbyte a, setofbyte b)
 *       for intersection
 *   setofb_eql( setofbyte a, setofbyte b)
 *       for comparision
 *   setofb_eq0( setofbyte a)
 *       for comparision with an empty set
 *   are avaluable
 *
 */

#define _E  (-1)		/* end of set marker */
#define UpTo  (-2)		/* UpTo is analog of .. in Pascal */
#define setrec enum

typedef unsigned setofbyte[16];

/* We suppouse here the length of enum variable is 2 bytes,
   i.e. the length of enum type includes 16 objects as maximum.
   Setofbyte is 32 bytes length for set of byte enum type.
 */


extern int setof (int i,...);
extern int inset (int a, int sp);
extern unsigned *setofb (int i,...);
extern unsigned *setofb_add1 (setofbyte a, int i);
extern int insetb (unsigned a, setofbyte sp);
extern void setofb_zero (setofbyte sp);
extern void setofb_cpy (setofbyte dest, setofbyte source);
extern unsigned *setofb_uni (setofbyte a, setofbyte b);
extern unsigned *setofb_aun (setofbyte a, setofbyte b);
extern unsigned *setofb_its (setofbyte a, setofbyte b);
extern int setofb_eql (setofbyte a, setofbyte b);
extern int setofb_eq0 (setofbyte a);
extern void setofb_dpl (setofbyte a);

#endif
