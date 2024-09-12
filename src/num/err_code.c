/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ---------------------------------------------------
*/
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "chep_crt/include/crt_util.h"

#include"err_code.h"

int err_code = 0;

void 
errormessage (void)
{
  switch (err_code)
    {
    case 1:
      warnanykey (10, 10, "Constraints can't be calculated");
      break;
    case 2:
      warnanykey (10, 10, "Zero denominator");
      break;

    case 3:
      warnanykey (10, 10, "Negative particle mass");
      break;

    case 4:
      warnanykey (10, 10, "Energy too small");
      break;

    case 5:
      warnanykey (10, 10, "Can not evaluate constraints");
      break;


    case 10:
      warnanykey (10, 10, "User Break");
      break;

    default:
      warnanykey (10, 10, " Error ?? ");
    }
}
