/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* ------------------------------------------------------------
*/
#ifndef __PARSER_
#define __PARSER_

typedef void *(*operation) (char *, int, void **);
typedef void *(*rdelement) (char *s);
typedef void (*clear) (void *);

#define braketexpected       1
#define unexpectedcharacter  2
#define operationexpected    3
#define toolongidentifier    4
#define toomanyagruments     5
#define cannotread           6
#define cannotevaluate       7

#define unknownidentifier     8
#define unexpectedoperation   9
#define unknownfunction      10
#define wrongnumberofagr     11
#define typemismatch         12

#define naninoperation       13

#define toolargenumber       14
#define indexuncompatibility 15
#define rangecheckerror      16


#define numbertp 0
#define polytp 1
#define rationtp 2
#define tenstp 3
#define spintp 4
#define vectortp 5
#define indextp 6

extern int rderrcode, rderrpos;


extern void *readExpression (char *source, rdelement rd, operation act, clear del);

#endif
