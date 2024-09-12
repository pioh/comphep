/*
* Copyright (C) 2001-2009, CompHEP Collaboration
* Copyright (C) 1997, Alexander Pukhov 
* ------------------------------------------------------
*/
#include "service2/include/chep_limits.h"
#include "service2/include/getmem.h"
#include "service2/include/tptcmac.h"
#include "service2/include/parser.h"
#include "service2/include/syst.h"

#include "pvars.h"
#include "pre_read.h"

preres pregarbage = NULL;

void 
clearpregarbage (void)
{
  while (pregarbage)
    {
      preres p = pregarbage->next;
      if (pregarbage->varsdeg)
        free (pregarbage->varsdeg);
      free (pregarbage);
      pregarbage = p;
    }

}


static void 
rangecheck (int varrange, preres q)
{
  int i;

  if (q->degp > 255 || q->maxg > 20)
    {
      rderrcode = rangecheckerror;
      return;
    }
  for (i = 0; i < varrange; i++)
    if (q->varsdeg[i] > 127)
      {
        rderrcode = rangecheckerror;
        return;
      }
}


static void 
newrecord (preres * q)
{
  int i;
  *q = pregarbage;
  while (*q != NULL && !(*q)->free)
    *q = (*q)->next;
  if (*q == NULL)
    {
      *q = (preres) m_alloc (sizeof (struct preresrecord));
      (*q)->next = pregarbage;
      pregarbage = *q;
      (*q)->nvar = 0;
      (*q)->varsdeg = NULL;
    }
  else
    for (i = 0; i < (*q)->nvar; i++)
      (*q)->varsdeg[i] = 0;

  (*q)->free = FALSE;
  (*q)->num = 0;
  (*q)->maxp = 0;
  (*q)->degp = 0;
  (*q)->g5 = FALSE;
  (*q)->maxg = 0;
  (*q)->indlist = setof (_E);
}


void * rd_pre (char *s) {
  int i, k;
  preres m;
  long li;

  if (isdigit (s[0]))
    {
      if (1 != sscanf (s, "%ld", &li))
        rderrcode = toolargenumber;
      else
        {
          newrecord (&m);
          m->num = li;
          m->tp = numbertp;
        }
    }
  else
    {
      if (strlen (s) > 6)
        rderrcode = toolongidentifier;
      else
        {
          newrecord (&m);
          m->tp = polytp;

          if (strlen (s) == 2 && isdigit (s[1]) && s[1] != '0')
            {
              switch (s[0])
                {
                case 'p':
                case 'P':
                  m->tp = vectortp;
                  m->maxp = s[1] - '0';
                  m->degp = 1;
                  break;

                case 'm':
                  m->tp = indextp;
                  m->indlist = setof (s[1] - '0', _E);
                  break;
                case 'M':
                  m->tp = indextp;
                  m->indlist = setof (s[1] - '0' + 4, _E);
                }               /*  Case  */

              if (strcmp (s, "G5") == 0)
                {
                  m->tp = spintp;
                  m->g5 = TRUE;
                }
            }

          if (m->tp == polytp)
            {
              i = 0;
              while (i < vardef->nvar && strcmp (s, vardef->vars[i].name))
                i++;
              if (i == vardef->nvar)
                {
                  increaseVars (vardef);
                  strcpy (vardef->vars[i].name, s);
                }
              m->nvar = ALIG (vardef->nvar);
              m->varsdeg = re_alloc (m->varsdeg, m->nvar * sizeof (unsigned));
              for (k = 0; k < m->nvar; k++)
                m->varsdeg[k] = 0;
              m->varsdeg[i] = 1;
            }
        }
    }
  return (void *) m;
}


static void *
uact (char *ch, void *mm)
{
  preres m;

  m = (preres) mm;
  if (strcmp (ch, "G") == 0 || strcmp (ch, "g") == 0)
    {
      if (!inset (m->tp, setof (vectortp, indextp, _E)))
        {
          rderrcode = typemismatch;
          return mm;
        }
      m->tp = spintp;
      m->maxg = 1;
    }
  else if (strcmp (ch, "-") == 0)
    {
      if (m->tp == indextp)
        {
          rderrcode = typemismatch;
          return mm;
        }
      m->num = -m->num;
    }
  else
    rderrcode = unexpectedoperation;
  return mm;
}


static void *
bact (int varrange, char ch, void *mm1, void *mm2)
{
  preres m1, m2, m3;
  int i;
  int d, k;

  m1 = (preres) mm1;
  m2 = (preres) mm2;

  if (m1->nvar < varrange)
    {
      int newsize = ALIG (varrange);
      m1->varsdeg = re_alloc (m1->varsdeg, newsize * sizeof (unsigned));
      for (i = m1->nvar; i < newsize; i++)
        m1->varsdeg[i] = 0;
      m1->nvar = newsize;
    }


  if (m2->nvar < varrange)
    {
      int newsize = ALIG (varrange);
      m2->varsdeg = re_alloc (m2->varsdeg, newsize * sizeof (unsigned));
      for (i = m2->nvar; i < newsize; i++)
        m2->varsdeg[i] = 0;
      m2->nvar = newsize;
    }

  switch (ch)
    {
    case '+':
      if (m1->tp < m2->tp)
        {
          m3 = m1;
          m1 = m2;
          m2 = m3;
        }
      if (m1->tp == indextp ||
          (m1->tp == vectortp && m2->tp != vectortp) ||
          (m2->tp > polytp && m1->tp != m2->tp))
        {
          rderrcode = typemismatch;
          return NULL;
        }

      if (m1->indlist != m2->indlist)
        {
          rderrcode = indexuncompatibility;
          return NULL;
        }
      m1->num += m2->num;
      m1->maxp = MAX (m1->maxp, m2->maxp);
      m1->degp = MAX (m1->degp, m2->degp);
      m1->g5 = m1->g5 || m2->g5;
      m1->maxg = MAX (m1->maxg, m2->maxg);
      if (m1->tp == rationtp)
        for (i = 0; i < varrange; i++)
          m1->varsdeg[i] += m2->varsdeg[i];
      else
        for (i = 0; i < varrange; i++)
          m1->varsdeg[i] = MAX (m1->varsdeg[i], m2->varsdeg[i]);
      break;

    case '*':
      if (m1->tp < m2->tp)
        {
          m3 = m1;
          m1 = m2;
          m2 = m3;
        }
      if (m1->tp == indextp || (m2->tp > polytp && m1->tp != m2->tp))
        {
          rderrcode = typemismatch;
          return NULL;
        }
      if ((m1->indlist & m2->indlist) != setof (_E))
        {
          rderrcode = indexuncompatibility;
          return NULL;
        }
      m1->num *= m2->num;
      m1->maxp = MAX (m1->maxp, m2->maxp);
      m1->degp += m2->degp;
      m1->g5 = m1->g5 || m2->g5;
      m1->maxg += m2->maxg;
      for (i = 0; i < varrange; i++)
        m1->varsdeg[i] += m2->varsdeg[i];
      m1->indlist |= m2->indlist;
      break;

    case '/':
      if (m2->tp > rationtp || m1->tp > rationtp || m1->tp == indextp)
        {
          rderrcode = typemismatch;
          return NULL;
        }

      m1->maxp = MAX (m1->maxp, m2->maxp);
      m1->degp += m2->degp;
      m1->g5 = m1->g5 || m2->g5;
      for (i = 1; i <= varrange; i++)
        m1->varsdeg[i - 1] += m2->varsdeg[i - 1];

      m1->tp = rationtp;
      break;

    case '^':
      d = m2->num;
      if (m2->tp != numbertp || d <= 0 || d > 255)
        {
          rderrcode = rangecheckerror;
          return NULL;
        }

      if (m1->tp > rationtp)
        {
          rderrcode = typemismatch;
          return NULL;
        }

      k = m1->num;
      for (i = 1; i <= d - 1; i++)
        m1->num *= k;
      m1->degp *= d;
      for (i = 1; i <= varrange; i++)
        m1->varsdeg[i - 1] *= d;
      break;

    case '.':
      if (m1->tp < vectortp || m2->tp < vectortp)
        {
          rderrcode = typemismatch;
          return NULL;
        }
      m1->tp = m1->tp == vectortp && m2->tp == vectortp ?
        polytp : tenstp;
      if (m1->indlist != setof (_E) && m1->indlist == m2->indlist)
        {
          rderrcode = indexuncompatibility;
          return NULL;
        }
      else
        m1->indlist |= m2->indlist;
      m1->num *= m2->num;
      m1->maxp = MAX (m1->maxp, m2->maxp);
      m1->degp += m2->degp;
      m1->g5 = m1->g5 || m2->g5;
      m1->maxg += m2->maxg;
      for (i = 0; i < varrange; i++)
        m1->varsdeg[i] += m2->varsdeg[i];
      break;
    default:
      rderrcode = unexpectedoperation;
      m1->free = 1;
      m2->free = 1;
      return NULL;
    }                           /*  Case  */

  m2->free = TRUE;

  rangecheck (varrange, m1);
  return (void *) m1;
}

void *
act_pre (char *name, int n, void **args)
{
  void *args2[2];
  int vrange = vardef->nvar;
  if (n == 2 && !strcmp (name, "-"))
    strcpy (name, "+");
  if (n == 1)
    return uact (name, args[0]);
  if (n == 2)
    return bact (vrange, name[0], args[0], args[1]);
  if (n == 4 && !strcmp (name, "eps"))
    {
      if (!(args2[0] = bact (vrange, '.', args[0], args[1])))
        return NULL;
      if (!(args2[1] = bact (vrange, '.', args[2], args[3])))
        return NULL;
      return bact (vrange, '*', args2[0], args2[1]);
    }
  rderrcode = unexpectedoperation;
  return NULL;
}

void *
act_preF (char *name, int n, void **args)
{
  if (name[0] == '+' || (n == 2 && name[0] == '-') || name[0] == '.')
    {
      rderrcode = unexpectedoperation;
      return NULL;
    }
  return act_pre (name, n, args);
}
