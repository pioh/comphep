/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __LIMITS_
#define __LIMITS_

#ifndef STRSIZ
#define STRSIZ 2048
#endif

#define PLISTLEN 10

typedef char vshortstr[16];
typedef char shortstr [128];
typedef char midstr   [1024];
typedef char longstr  [2048];

#define VSHORTSTRLEN 16
#define SHORTSTRLEN 128
#define MISSTRLEN 1024
#define LONGSTRLEN 2048

#define MAX_DIM   15
#define MAX_NDMX  50

#define SORTARR(arr,len) {int I=1;long L;  while (I < (len)) \
if (arr[I-1] <= arr[I]) I++; else {L=arr[I-1];arr[I-1]=arr[I];arr[I]=L;if(I>1)I--;}}

#ifndef MAX
#define MAX(a,b)        (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x)          (((x) < 0 ? -(x) : (x)))
#endif

#ifndef FALSE
#define FALSE   0
#endif

#ifndef TRUE
#define TRUE    1
#endif

#define FORMAT " %17.10E"
#define MAXNP 20
#define NAMELEN 20

#define ENERGY(m,p) sqrt((m)*(m) + *(p)**(p) + *(p+1)**(p+1)+ *(p+2)**(p+2))

#define MAXINOUT 10

#ifndef M_E
#define M_E             2.7182818284590452354
#endif
#ifndef M_LOG2E
#define M_LOG2E         1.4426950408889634074
#endif
#ifndef M_LOG10E
#define M_LOG10E        0.43429448190325182765
#endif
#ifndef M_LN2
#define M_LN2           0.69314718055994530942
#endif
#ifndef M_LN10
#define M_LN10          2.30258509299404568402
#endif
#ifndef M_PI
#define M_PI            3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2          1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4          0.78539816339744830962
#endif
#ifndef M_1_PI
#define M_1_PI          0.31830988618379067154
#endif
#ifndef M_2_PI
#define M_2_PI          0.63661977236758134308
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI      1.12837916709551257390
#endif
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2       0.70710678118654752440
#endif

#endif
