#ifndef __PARSER_
#define __PARSER_

typedef void* ( * binaction)(char ch, void * m1, void * m2);
typedef void* ( * unaction )(char * ch, void * m);
typedef void* ( * rdelement)(char * s);

#define ok   0
#define recordabsent 1
#define braketexpected 2
#define syntaxerror 3
#define unexpectedcharacter 4
#define operationexpected 5
#define unbalancedbraket 6
#define toolargenumber 7
#define toolongidentifier 8
#define typemismatch 9
#define indexuncompatibility 10
#define toomanyidentifiers 11
#define rangecheckerror 12
#define divisionbyzero 13
#define negativsqrtarg 14
#define unexpectedoperation 15

#define PlusWrongInFacor 18

#define numbertp 0
#define polytp 1
#define rationtp 2
#define tenstp 3
#define spintp 4
#define vectortp 5
#define indextp 6

extern int  rderrcode, rderrpos;


extern void *  readExpression(char * source, binaction  bact,
  unaction uact, rdelement rd);

#endif
