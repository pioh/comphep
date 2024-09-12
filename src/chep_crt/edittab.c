/* 
* Copyright (C) 1997-2009, CompHEP Collaboration
* Author: A.Kryukov
* ---------------------------------------------------------------
*/
/****************************************************************
           Coordinates and etc.  
****************************************************************
               X0   Xt                   Xt
   +----------------------------------->-------------------+
   |           :    :                      Table           |
   |           :    :                                      |
   |           X1   X     X1        Xm     X2   x          |
   |     +------------------------------------->---+       |
   |     |     :    :     :         :      : Screen|       |
   |     |     :    :     :         :      :       |       |
 Y0|...Y1|.....+---------------------------+       |       |
   |     |     |    :     :         :  Win |       |       |
   |     |     |    :     :         :      |       |       |
   |   Y |.....|....CursorF         :      |       |       |
   |     |     |                    :      |       |       |
   |     |     |                    :      |       |       |
   |     |     |                    :      |       |       |
   |   Ym|.....|....................Mouse  |       |       |
   |     |     |                           |       |       |
   |   Y2|.....+---------------------------+       |       |
   |     |                                         |       |
   |     V                                         |       |
   |   y |                                         |       |
   |     +-----------------------------------------+       |
   |                                                       |
   V                                                       |
 Yt|                                                       |
   +-------------------------------------------------------+

****************************************************************/
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#include "service2/include/chep_limits.h"
#include "service2/include/syst.h"
#include "service2/include/getmem.h"

#include "file_scr.h"
#include "crt.h"
#include "crt_util.h"
#include "help.h"
#include "edittab.h"

#define	MPRESSED	1
#define	MRELEASED	2
#define EDITTAB_DEBUG   0

#define SE_HLP          "h_003107"
#define ZE_HLP          "h_003108"
#define F2_HLP          "h_003109"

#define SEPR_CHAR ('|')
#define RESIZE_CHAR '<'
#define KBD_BUFF_SIZE 16
#define CTRL(x) ((x)-'A'+1)

#define X2XT(i)  (c_tab.x0+(i)-X1-1)
#define XT2X(k)  ((k)-c_tab.x0+X1+1)

typedef struct table0
  {
    table *tab;
    linelist top, edit, del;	/* top line */
    int x1, y, x2;		/* x1(x2)-begin(end)  of field ,y - cursos position, */
    int x0, y0;			/* first col. and line that displayed */
    int width;
  }
table0;

static table0 c_tab;		/* current edited table */

static int scrCur;		/* screen cursor */

#define TABCUR  (scrCur - X1+c_tab.x0-2)	/* table cursor */
/*************** Work with boxes ***************/

typedef struct box_struct
  {
    int x1, y1;			/* up-left conner */
    int x2, y2;			/* down-right conner */
    int fg, bg;			/* fore-, back- ground */
    void *buff;			/* buffer to save background */
  }
box_struct;

static box_struct e_box;	/* main window */


/****************** static vars. *******************/

static int X1 = 1;		/* up-left conner of table window */
static int Y1 = 1;
static int X2 = 80;		/* down-right conner of table window */
static int Y2 = 24;

static int maxX2 = 80;
static int maxY2 = 24;
static int x2t = 0;		/* width of table (exclude space char.) */

static int kbd_buff[KBD_BUFF_SIZE] =
{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static int kbd_p1 = 0;		/* out going position in kbd_buff */
static int kbd_p2 = 0;		/* in comming position of kbd_buff */


/***********************************************************************/
/*************** Cursor ON/OFF *********************** Jul. 05, 1996 ***/
/***********************************************************************/

static void 
drawZCur (char c, int x, int y)
/*
   Draw cursor in Zoom mode
   A.K. Feb, 9, 1998.
 */
{
  scrcolor (White, Black);
  goto_xy (x, y);
  print ("%c", c);
  scrcolor (Black, White);
  goto_xy (x, y);
/*      fprintf(stderr,"**** drawZCur(%c,%d,%d)\n",c,x,y); */
/**/
}

static void clearZCur(char c,int x,int y)
/*
   Clear cursor in Zoom mode
   A.K. Feb, 9, 1998.
*/
  {
    scrcolor (Black, White);
    goto_xy (x, y);
    print ("%c", c);
    goto_xy (x, y);
/*      fprintf(stderr,"**** clearZCur(%c,%d,%d)\n",c,x,y); */
/**/
}

/***********************************************************************/
/*************** Horizontal button bar *************** Jan. 31, 1995 ***/
/***********************************************************************/

#define HBB_MAX_NLEN 8

    typedef struct hbb_item
      {
	char name[HBB_MAX_NLEN];	/* name of item */
	int rc;			/* return code */
	int status;		/* 0 - disable, 1 - enable. */
      }
    hbb_item;

    typedef struct hbb_struct
      {
	hbb_item *hbbar;	/* pointer to array of hbbar_item */
	int nb;			/* number of buttons */
	char bg;		/* background char */
	int x, y;		/* position */
	int status;		/* 0 - disable, 1 - enable. */
	int onscreen;		/*0 - offscreen, 1 - onscreen */
	int fgc, bgc;		/* foreground and background colours */
	int fgc_off, bgc_off;	/* fg and bg colours for off screen */
      }
    hbb_struct;

    static int hbb_new (hbb_struct * hbb, hbb_item * hbbar, char bg)
    {
      int nb;
      for (nb = 0; strcmp (hbbar[nb].name, ""); nb++)
	{
	}
      hbb->hbbar = hbbar;
      hbb->nb = nb;
      hbb->bg = bg;
      hbb->x = 1;		/* default value, see alsow hbb_place */
      hbb->y = 1;
      hbb->status = 1;
      hbb->onscreen = 0;
      hbb->fgc = White;
      hbb->bgc = Black;
      hbb->fgc_off = Black;
      hbb->bgc_off = White;
      return 0;
    }

    static int hbb_place (hbb_struct * hbb, int x, int y)
    {
      hbb->x = x;
      hbb->y = y;
#if EDITTAB_DEBUG
      printf ("hbb_place: hbb->x=%d, x=%d\n", hbb->x, y);
#endif
      return 0;
    }

    static int hbb_open (hbb_struct * hbb)
    {
      int x0 = where_x (), y0 = where_y ();
      int i;
      if (hbb->status == 0 || hbb->onscreen)
	  return 0;
        goto_xy (hbb->x, hbb->y);
        scrcolor (hbb->fgc, hbb->bgc);
      for (i = 0; i < hbb->nb; i++)
	{
	  if (where_x () + strlen (hbb->hbbar[i].name) + 1 >= X2)
	    {
	      hbb->nb = i;
	      break;
	    }
	  print ("%s", hbb->hbbar[i].name);
	  goto_xy (where_x () + 1, hbb->y);
	}

      scrcolor (hbb->fgc_off, hbb->bgc_off);
      goto_xy (x0, y0);
      hbb->onscreen = 1;
      return 0;
    }

    static int hbb_close (hbb_struct * hbb)
    {
      int x0 = where_x (), y0 = where_y ();
      int i, j;
      if (hbb->status == 0 || !hbb->onscreen)
	  return 0;
        goto_xy (hbb->x, hbb->y);
        scrcolor (hbb->fgc_off, hbb->bgc_off);
      for (i = 0; i < hbb->nb; i++)
	for (j = 0; j <= strlen (hbb->hbbar[i].name); j++)
	    print ("%c", hbb->bg);
        goto_xy (x0, y0);
        hbb->onscreen = 0;
        return 0;
    }

    static int hbb_select (hbb_struct * hbb)
    {
      int i, x;
      if (!hbb->status || !hbb->onscreen)
	  return 0;
      if (mouse_info.row != hbb->y)
	  return 0;		/* out of bar */
        x = mouse_info.col - hbb->x;
      if (x < 0)
	  return 0;
      for (i = 0; hbb->hbbar[i].status; i++)
	{
	  x -= strlen (hbb->hbbar[i].name);
	  if (x < 0)
	    return hbb->hbbar[i].rc;
	  if (--x < 0)
	    return 0;
	}

      return 0;
    }

    static void hbb_colors (hbb_struct * hbb, int fgc, int bgc, int fgc_off, int bgc_off)
    {
      hbb->fgc = fgc;
      hbb->bgc = bgc;
      hbb->fgc_off = fgc_off;
      hbb->bgc_off = bgc_off;
    }

/***********************************************************************/
/********* List of horizontal button bars *********** Feb. 04, 1995 ****/
/***********************************************************************/

    typedef struct hbl_struct
    {
      struct hbl_struct *next;
      hbb_struct *hbb;
    }
    hbl_struct;

    static hbl_struct *hbl_add (hbl_struct * hbl, hbb_struct * hbb)
    {
      hbl_struct *new;
        new = (hbl_struct *) m_alloc (sizeof (hbl_struct));
      if (!new)
	  exit (99);
        new->hbb = hbb;
        new->next = hbl;
        return new;
    }

    static int hbl_delete (hbl_struct * hbl)
    {
      hbl_struct *w;
      for (; hbl;)
	{
	  w = hbl;
	  hbl = hbl->next;
	  free (w);
	};
      return 0;
    }

    static int hbl_select (hbl_struct * hbl)
    {
      hbl_struct *w;
      int rc = 0;
      for (w = hbl; w && !rc; w = w->next)
	  rc = hbb_select (w->hbb);
#if EDITTAB_DEBUG
        printf ("hbl_select: rc=%d\n", rc);
#endif
        return rc;
    }

    static int hbl_open (hbl_struct * hbl)
    {
      hbl_struct *w;
      for (w = hbl; w; w = w->next)
	  hbb_open (w->hbb);
        return 0;
    }

    static int hbl_close (hbl_struct * hbl)
    {
      hbl_struct *w;
      for (w = hbl; w; w = w->next)
	  hbb_close (w->hbb);
        return 0;
    }

    static void hbl_colors (hbl_struct * hbl, int fgc, int bgc, int fgc_off, int bgc_off)
    {
      hbl_struct *w;
      for (w = hbl; w; w = w->next)
	{
	  hbb_colors (w->hbb, fgc, bgc, fgc_off, bgc_off);
	};
    }

/*********************** statics *********************/

    static hbb_item up_items[6] =
    {
      {"Clr", CTRL ('C'), 1},
      {"Rest", CTRL ('R'), 1},
      {"Del", CTRL ('D'), 1},
      {"Size", CTRL ('S'), 1},
      {"", 0, 0}
    };

    static hbb_item down_items[9] =
    {
      {"F1", KB_F1, 1},
      {"F2", KB_F2, 1},
      {"Top", 20, 1},
      {"Bottom", 2, 1},
      {"GoTo", CTRL ('G'), 1},
      {"Find", CTRL ('F'), 1},
      {"Zoom", CTRL ('Z'), 1},
      {"ErrMes", CTRL ('E'), 1},
      {"", 0, 0}
    };

    static hbb_struct up_bar, down_bar;
    static hbl_struct *hbl;

/****************************************************************/
/************* Mouse functions *************** Feb. 04, 1995 ****/
/****************************************************************/

    static int m_inbox (int x1, int y1, int x2, int y2)
    {
      return (x1 < mouse_info.col && mouse_info.col < x2
	      && y1 < mouse_info.row && mouse_info.row < y2);
    }

/**************** Stack command **************/

    static void push_key (int kk)
    {
      kbd_buff[kbd_p2] = kk;
      kbd_p2++;
      if (kbd_p2 == KBD_BUFF_SIZE)
	kbd_p2 = 0;
      if (kbd_p2 == kbd_p1)
	exit (99);
    }

    static int pop_key (void)
    {
      int kk;
      if (kbd_p1 != kbd_p2)
	{
	  kk = kbd_buff[kbd_p1];
	  kbd_p1++;
	  if (kbd_p1 == KBD_BUFF_SIZE)
	    kbd_p1 = 0;
#if EDITTAB_DEBUG
	  printf ("kbd_buff: p1=%d, kk=%d, %d\n"
		  ,kbd_p1, kk
		  ,kbd_buff[kbd_p1]);	/* debug */
#endif
	}
      else
	{
	  kk = inkey ();
	  if (kbd_p2 == KBD_BUFF_SIZE)
	    kbd_p2 = 0;
	  if (kbd_p1 == KBD_BUFF_SIZE)
	    kbd_p1 = 0;
	  switch (kk)
	    {
	    case KB_UP:
	    case KB_DOWN:
	      push_key (1);
	      break;

	    case KB_LEFT:
	    case KB_RIGHT:
	      push_key (1);
	      break;
	    };
	}
      return kk;
    }

/*********************************************************************/
/***************** Draw table ********************** Feb. 04, 1995 ***/
/*********************************************************************/

    static void redrawTable (table * tab)
    {
      int i;
      linelist lastline = c_tab.top;
      char buffer[STRSIZ];

        strcpy (buffer, tab->format + c_tab.x0 - 1);
        buffer[c_tab.width] = 0;

        scrcolor (Black, White);
        goto_xy (X1 + 1, Y1 + 1);
        print ("%s", buffer);

      for (i = Y1 + 2; i < Y2; i++)
	{
	  goto_xy (X1 + 1, i);
	  if (lastline != NULL)
	    {
	      strcpy (buffer, lastline->line + c_tab.x0 - 1);
	      buffer[c_tab.width] = 0;
	      print (buffer);
	      lastline = lastline->next;
	    }
	  print ("%*.*s", X2 - where_x (), X2 - where_x (), "");
	}

    }				/* redrawTable */

/*************** Shift cursor ***********************/

    static int shiftY1 (int y)
    {
      int res = 0;
      switch (y)
	{
	case 1:
	  if (c_tab.edit->next)
	    {
	      c_tab.edit = c_tab.edit->next;
	      c_tab.y++;
	    }
	  break;
	case -1:
	  if (c_tab.edit->pred)
	    {
	      c_tab.edit = c_tab.edit->pred;
	      c_tab.y--;
	    }
	  break;
	}
      return res;
    }				/* shiftY */

    static int shiftY (int y)
    {
      int i, j = 1;
      if (y < 0)
	{
	  j = -1;
	  y = -y;
	};
      for (i = 0; i < y; i++)
	{
	  shiftY1 (j);
	};
      return 0;
    }
#define MINWIDTH 15		/* desirable minimum width of table */

    static void setFieldFrame (int xPos, char *format)
    {
      int z;
      for (z = xPos - X1 + c_tab.x0 - 1; z < x2t && format[z] != SEPR_CHAR; z++);
        c_tab.x2 = z + X1 - c_tab.x0 + 1;

      for (z = xPos - X1 + c_tab.x0 - 3; z >= 0 && format[z] != SEPR_CHAR; z--);
        c_tab.x1 = z + X1 - c_tab.x0 + 3;
    }



    static void initEdit (table * tab, int xx1, int yy1, int ntop)
    {
      int i;
      /* Coord. of box on the screen */
        X1 = xx1;
        Y1 = yy1 + 1;
        maxX2 = maxCol ();
        maxY2 = maxRow ();
        X2 = maxX2;
        Y2 = maxY2;
      /* Initialization of static */
        c_tab.tab = tab;

      for (x2t = strlen (tab->format) - 1; x2t && tab->format[x2t] != SEPR_CHAR; x2t--)
	{
	};
      tab->format[x2t + 1] = 0;
      X2 = X1 + strlen (tab->format);
      X2 = MAX (X2, X1 + MINWIDTH);
      if (X2 > maxX2)
	X2 = maxX2;
      c_tab.top = tab->strings;
      if (!c_tab.top)
	ntop--;
      for (i = 1; (i < ntop) && (c_tab.top->next); i++)
	c_tab.top = c_tab.top->next;
      c_tab.edit = c_tab.top;
      c_tab.del = NULL;
      c_tab.del = (linelist) m_alloc (sizeof (linerec) + 2);
      c_tab.del->pred = NULL;
      c_tab.del->next = NULL;
      strcpy (c_tab.del->line, tab->format);
      for (i = 0; c_tab.del->line[i]; i++)
	{
	  if (c_tab.del->line[i] != SEPR_CHAR)
	    c_tab.del->line[i] = ' ';
	}
      c_tab.x0 = 1;
      c_tab.y0 = ntop;
      c_tab.width = X2 - X1 - 1;
      c_tab.x1 = X1 + 1;
      c_tab.y = Y1 + 2;
      c_tab.x2 = 0;
      scrCur = X1 + 1;
      setFieldFrame (scrCur, tab->format);
    }				/* initEdit */

#define TEMPORAL (-100)
    static int s_edit (int x1, int x2, int x0, int y1, char *s, int show)
    {
      /* x1,x2,y1 - coord of string window (begin,end,Y)
         x0 - the first char. from s which display 
         Return o if exit by [esc]    A.K.   Feb. 9, 1998
       */
      int i;
      int rc = 0, rc1 = 0;
      int key = 0;
      int len, z, slen0 = strlen (s);
      int dx, dy;

        len = x2 - x1 + 1;

      while (key != KB_ESC && key != TEMPORAL)
	{

	  if (show && (((key >= 32) && (key <= 127))
		       || key == KB_BACKSP || key == KB_DC || key == KB_ENTER
		       || key == CTRL ('C') || key == CTRL ('D')
		       || key == CTRL ('S') || key == CTRL ('R')))
	    {
	      be_be ();
	      key = 0;
	    }

	  switch (key)
	    {

	    case KB_BACKSP:
	      if (scrCur > x1)
		scrCur--;
	      else
		break;
	    case KB_DC:
	      for (i = x0 + scrCur - x1; s[i]; i++)
		s[i] = s[i + 1];
	      s[i - 1] = ' ';
	      rc1 = 1;
	      break;

	    case KB_LEFT:
	      dx = pop_key ();
	      scrCur--;
	      break;

	    case KB_RIGHT:
	      dx = pop_key ();
	      scrCur++;
	      break;

	    case KB_HOME:
	      scrCur = x1;
	      x0 = 0;
	      break;

	    case KB_END:
	      z = strlen (s) - x0 + x1 - x2 - 1;
	      if (z > 0)
		{
		  scrCur = x2;
		  x0 += z;
		}
	      else
		scrCur = z + x2;
	      break;

	    case KB_ENTER:
	      push_key (TEMPORAL);
	      push_key (KB_ENTER);
	      break;

	    case KB_UP:
	    case KB_DOWN:
	      dy = pop_key ();
	      push_key (TEMPORAL);
	      push_key (key);
	      push_key (dy);
	      break;

	    case KB_TAB:;
	      push_key (TEMPORAL);
	      push_key (key);
	      break;

	    case CTRL ('C'):
	      for (i = x0 + scrCur - x1; i < strlen (s); i++)
		s[i] = ' ';
	      rc1 = 1;
	      break;

	    case KB_PAGEU:
	    case KB_PAGED:
	    case CTRL ('B'):
	    case CTRL ('D'):
	    case CTRL ('T'):
	    case CTRL ('E'):
	    case CTRL ('F'):
	    case CTRL ('G'):
	    case CTRL ('S'):
	    case CTRL ('Z'):
	    case KB_F1:
	    case KB_F2:
	      push_key (TEMPORAL);
	      push_key (key);
	      break;


	    case CTRL ('R'):
	      rc1 = 0;
	      push_key (TEMPORAL);
	      push_key (key);
	      break;

	    case KB_MOUSE:
	      if (mouse_info.but1 == MRELEASED)
		{
#if EDITTAB_DEBUG
		  printf (" mouse: but1=%d, x=%d, y=%d\n",
			  mouse_info.but1, mouse_info.row,
			  mouse_info.col);
		  printf (" X1=%d, Y1=%d\n", X1, Y1);
		  printf (" s-cursor: x=%d\n", x);
#endif
		  if (mouse_info.col == X1 && mouse_info.row == Y1 - 1)
		    push_key (KB_ESC);
		  else if ((rc = hbl_select (hbl /*->next*/ )))
		    push_key (rc);
		  else if (m_inbox (x1 - 1, y1 - 1, x2 + 1, y1 + 1))
		    scrCur = mouse_info.col;
		  else if (m_inbox (e_box.x1, e_box.y1, e_box.x2, e_box.y2))
		    {
		      push_key (TEMPORAL);
		      push_key (KB_MOUSE);
		    };
		};
	      key = 0;
	      break;

	    default:
	      if ((key >= 32) && (key <= 127))
		{
		  if (s[strlen (s) - 1] == ' ' || strlen (s) < slen0)
		    {
		      for (i = strlen (s) - 2; i >= x0 + scrCur - x1; i--)
			s[i + 1] = s[i];
		      s[x0 + scrCur - x1] = (char) key;
		      /*if (scrCur<x2) */ scrCur++;
		      /*else {if (s[x0+scrCur-x1+1]) x0++;} */
		    }
		  else
		    be_be ();
		  rc1 = 1;
		}
	      key = 0;
	      break;
	    }			/* switch */

	  if (scrCur > x2)
	    {
	      push_key (TEMPORAL);
	      push_key (KB_RIGHT);
	      push_key (1);
	    }
	  else if (scrCur < x1)
	    {
	      push_key (TEMPORAL);
	      push_key (KB_LEFT);
	      push_key (1);
	    }
	  else
	    {
	      scrcolor (White, Blue);
	      goto_xy (x1, y1);
	      print ("%.*s", len, s + x0);
	      scrcolor (Blue, Yellow);
	      goto_xy (scrCur, y1);
	      print ("%c", s[x0 + scrCur - x1]);
	    }
	  goto_xy (scrCur, y1);

	  key = pop_key ();
	}			/* while */

/*  if (key == KB_ENTER) rc1 = 1;   A.K.   Feb. 9, 1998 */
      if (key == KB_ESC)
	{
	  push_key (KB_ESC);	/* rc1 = 0; */
	}			/*  A.K.   Feb. 9, 1998 */

      scrcolor (Black, White);
      goto_xy (x1, y1);
      print ("%.*s", len, s + x0);
      return rc1;
    }				/* s_edit */

/****************** Zoom *******************/

    static int wZoom;		/* Width of zoom window without border */
    static int hZoom;		/* High of zoom window without border */
    static int sLen;		/* The length of the zoomed string */
    static int sEnd;		/* The real end of the zoomed string */

    static int xy2z (int x, int y, int *z)
    {
/*                                        May 12, 1994
   xy2z transorm 2-dim coordinate to linear coordinate along
   the working string.

   x,y - 2-dim. coordinate. Start point is (0,0). Last on is 
   (hZoom-1,wZoom-1).
   *z   - linear coordinate. Start point is 0. Last one is sLen-1.

   return - 0 if ok and 1 if tranformation impossible.
 */
      int z1;			/* Work vars. */

        z1 = wZoom * y + x;
      if (z1 < 0)
	  return -1;
      else if (z1 >= sLen)
	  return 1;
      else
	{
	  *z = z1;
	  return 0;
	};
    }

    static int z2xy (int z, int *x, int *y)
    {
/*                                        May 12, 1994
   z2xy transorm linear coordinate along the working string to 
   2-dim. coordinate.

   *x,*y - 2-dim. coordinate. Start point is (0,0). Last on is 
   (hZoom-1,wZoom-1).
   z   - linear coordinate. Start point is 0. Last one is sLen-1.

   return - 0 if ok and 1 if tranformation impossible.
 */
      int x1, y1;		/* Work vars. */

        y1 = z / wZoom;
        x1 = z - (wZoom * y1);
      /*printf("*** z->(x,y)=%d->(%d,%d)",z,*x,*y); inkey(); */
      if ((x1 < 0) || (x1 >= wZoom) || (y1 < 0) || (y1 >= hZoom))
	  return 1;
      else
	{
	  *x = x1;
	  *y = y1;
	  return 0;
	};
    }


    static void drawLStr (int x, int y, char *s, int z)
    {
/*
   drawLStr draw logn string in the window prepare before.
   x,y - the start point on the screen.
   s   - drawable string.
   z   - number of intial char. from s.
 */
      int i;			/* cycle counter */
      int x1 = 0, y1 = 0;	/* local coordinate */
      char ws[STRSIZ];		/* work string */

      if (z >= sEnd)
	  sEnd = z + 1;
      if (z2xy (z, &x1, &y1))
	  exit (99);
        z = y1 * wZoom;
      for (; y1 < hZoom && z < sEnd + 1; y1++)
	{
	  for (i = x1; (i < wZoom) && s[z + i] && z + i < sEnd + 1; i++)
	    ws[i] = s[y1 * wZoom + i];
	  ws[i] = 0;
	  goto_xy (x + x1, y + y1);
	  print (&(ws[x1]));
	  x1 = 0;
	  z += wZoom;
	};
    }


    static void clZWin (int x, int y)
    {
/*
   clZein clear window for zoom.
   x,y - the start point on the screen.
   s   - drawable string.
   z   - number of intial char. from s.
 */
      int i;			/* cycle counter */
      int y1 = 0;		/* local coordinate */
      char ws[STRSIZ];		/* work string */

      for (; y1 < hZoom; y1++)
	{
	  for (i = 0; (i < wZoom) && (y1 * wZoom + i) < sLen; i++)
	    ws[i] = ' ';
	  if (y1 == (hZoom - 1))
	    for (; i < wZoom; i++)
	      ws[i] = '<';
	  ws[i] = 0;
	  goto_xy (x, y + y1);
	  print (ws);
	};
    }

    static int insChar (char c, char *s, int z)
    {
      int i;			/* Cycle var. */
      int l;			/* length of str. */

        l = strlen (s);
      if (s[l - 1] != ' ')
	  return 1;
      for (i = l - 2; i >= z; i--)
	  s[i + 1] = s[i];
        s[z] = c;
        sEnd++;
        return 0;
    }

    static int delChar (char *s, int z)
    {
      int i;			/* Cycle var. */
      int l;			/* length of str. */

        l = strlen (s);
      for (i = z + 1; i < l; i++)
	  s[i - 1] = s[i];
        s[l - 1] = ' ';
        sEnd--;
        return 0;
    }

    static int zoomStr (char *ws, int x, int y, int b, int show)
    {
/*                                      May 12, 1994 
   zoom shows long fields in separate window.

   ws - zoomed string
   x,y - the coordinates of left-top angle
   b - the width of border;
   show - showed flag.
   return - 0

   Modifications:
   A.K. Feb, 9, 1998.      clear and draw cursor.
   A.K. Feb. 9, 1998       restore old string if press [esc]
 */
      void *pbuff;		/* pointer to saving buffer */
      int key;			/* code of pressed key */
      int x1, y1;		/* shifted screen coordinate. the first position 
				   of the text */
      int edit = 0;		/* Is it edited? Not implemented */
      int zc, zc1;		/* Linear curses coordinate. zc1-previous pos. */
      int xc, yc, xc1, yc1;	/* Local curses coordinate. xc1,yc1 - previous pos. */
      int dx, dy;
      int rc;

        wZoom = 75;
      if (strlen (ws) >= STRSIZ)
	  exit (99);
      if (show)
	  trim (ws);
        sLen = strlen (ws);
      for (sEnd = sLen - 1; sEnd && ws[sEnd] == ' '; sEnd--)
	{
	};
      sEnd++;
      hZoom = STRSIZ / wZoom + 1;
      z2xy (sLen - 1, &xc, &yc);
      hZoom = yc;
      if (!hZoom)
	wZoom = xc + 1;
      /*printf("*** (x,y)=(%d,%d)",x,y); inkey(); */
      hZoom++;
      x1 = x + b;
      y1 = y + b;
      zc = 0;
      z2xy (zc, &xc, &yc);

/*get_text(x,y,x+wZoom-1+2*b,y+hZoom-1+2*b,&pbuff);     ???? */
      get_text (x, y, x + wZoom - 1 + 2 * b, y + hZoom + 2 * b, &pbuff);
/*printf(" x=%d,y=%d,wZoom=%d,hZoom=%d,b=%d\n",x,y,wZoom,hZoom,b); */
      scrcolor (Black, White);
      if (b > 0)
	chepbox (x, y, x + wZoom - 1 + 2 * b, y + hZoom - 1 + 2 * b);

      key = 0;
      clZWin (x1, y1);

      for (; key != KB_ESC && key != KB_ENTER;)
	{
	  zc1 = zc;
	  switch (key)
	    {
	    case 0:
	      drawLStr (x1, y1, ws, 0);
	      drawZCur (ws[0], x1 + xc, y1 + yc);
	      goto_xy (x1 + xc, y1 + yc);
	      break;

	    case KB_LEFT:
	      dx = pop_key ();
	      for (; dx > 0 && zc > 0; dx--)
		{
		  zc--;
		  z2xy (zc, &xc, &yc);
		  goto_xy (x1 + xc, y1 + yc);
		};
	      break;

	    case KB_RIGHT:
	      dx = pop_key ();
	      for (; dx > 0 /*&& zc<sLen-1 */ ; dx--)
		{
		  zc++;
		  if (zc == sLen - 1)
		    zc = 0;
		  z2xy (zc, &xc, &yc);
		  goto_xy (x1 + xc, y1 + yc);
		};
	      break;

	    case KB_UP:
	      dy = pop_key ();
	      for (; dy > 0 && yc > 0; dy--)
		{
		  yc--;
		  xy2z (xc, yc, &zc);
		  goto_xy (x1 + xc, y1 + yc);
		};
	      break;

	    case KB_DOWN:
	      dy = pop_key ();
	      for (; dy > 0 && yc < hZoom - 1; dy--)
		{
		  yc++;
		  if (xy2z (xc, yc, &zc))
		    yc--;
		  goto_xy (x1 + xc, y1 + yc);
		};
	      break;

	    case KB_PAGEU:
	    case KB_HOME:
	      zc = 0;
	      z2xy (zc, &xc, &yc);
	      goto_xy (x1 + xc, y1 + yc);
	      break;

	    case KB_PAGED:
	    case KB_END:
	      zc = sLen - 1;
	      z2xy (zc, &xc, &yc);
	      goto_xy (x1 + xc, y1 + yc);
	      break;

	    case KB_DC:	/* Del. char. */
	      if (show)
		break;
	      if (!delChar (ws, zc))
		edit = 1;
	      drawLStr (x1, y1, ws, zc);
	      goto_xy (x1 + xc, y1 + yc);
	      break;

	    case KB_BACKSP:	/* Del. char. before */
	      if (show)
		break;
	      if (zc <= 0)
		break;
	      zc--;
	      z2xy (zc, &xc, &yc);
	      goto_xy (x1 + xc, y1 + yc);
	      if (!delChar (ws, zc))
		edit = 1;
	      drawLStr (x1, y1, ws, zc);
	      goto_xy (x1 + xc, y1 + yc);
	      break;

	    case KB_MOUSE:
	      if (mouse_info.but1 == MRELEASED)
		{
#if EDITTAB_DEBUG
		  printf (" mouse: but1=%d, x=%d, y=%d\n",
			  mouse_info.but1, mouse_info.row,
			  mouse_info.col);
		  printf (" Z-cursor: xc=%d, yc=%d, zc=%d\n", xc, yc, zc);
#endif
		  if ((rc = hbl_select (hbl->next)))
		    {
		      push_key (rc);
		      if (rc == KB_UP)
			push_key (1);
		      else if (rc == KB_DOWN)
			push_key (1);
		      else if (rc == KB_LEFT)
			push_key (1);
		      else if (rc == KB_RIGHT)
			push_key (1);
		    }
		  else if (m_inbox (x, y, x + wZoom - 1 + 2 * b, y + hZoom - 1 + 2 * b))
		    {
		      dx = mouse_info.col - x1 - xc;
		      if (dx > 0)
			{
			  push_key (KB_RIGHT);
			  push_key (dx);
			}
		      else if (dx < 0)
			{
			  push_key (KB_LEFT);
			  push_key (-dx);
			};
		      dy = mouse_info.row - y1 - yc;
		      if (dy > 0)
			{
			  push_key (KB_DOWN);
			  push_key (dy);
			}
		      else if (dy < 0)
			{
			  push_key (KB_UP);
			  push_key (-dy);
			};
		    };
		};
	      break;

	    case KB_F1:
	      show_help (ZE_HLP);
	      scrcolor (Black, White);
	      break;

	    case KB_F2:
	      show_help ("zoomed");
	      scrcolor (Black, White);
	      break;

	    default:
	      if (show)
		break;
	      if ((key >= 32) && (key <= 127))
		{
		  if (insChar ((char) key, ws, zc))
		    break;
		  edit = 1;
		  drawLStr (x1, y1, ws, zc);
		  if (zc < sLen - 1)
		    zc++;
		  z2xy (zc, &xc, &yc);
		  goto_xy (x1 + xc, y1 + yc);
		};
	      break;
	    };			/* switch (key) */

	  if (!show)
	    {
	      goto_xy (x + 1, y);
	      scrcolor (White, Black);
	      print ("%d%c%c%c", zc + 1, boxFrame[5], boxFrame[5], boxFrame[5]);
	      goto_xy (x1 + xc, y1 + yc);
	      drawZCur (ws[zc], x1 + xc, y1 + yc);
	      if (zc != zc1)
		{
		  z2xy (zc1, &xc1, &yc1);
		  clearZCur (ws[zc1], x1 + xc1, y1 + yc1);
		};
	      scrcolor (Black, White);
	    }
	  key = pop_key ();
	};			/* while */

      if (key == KB_ESC)
	edit = 0;		/* A.K. Feb. 9, 1998 */

      scrcolor (White, Black);
      put_text (&pbuff);
      return edit;
    }
/*************** Resize Table *************/

    static void squeezeField (int k)
    {
      /* Squeeze the current field in table */
      /* k - how much */
      linelist llist;
      int i;
      char *s = c_tab.tab->format;

      for (i = c_tab.x2 + 1; s[X2XT (i) - 1]; i++)
	  s[X2XT (i - k) - 1] = s[X2XT (i) - 1];
        s[X2XT (i - k) - 1] = 0;
        s[X2XT (c_tab.x2 - k) - 1] = RESIZE_CHAR;
        s[X2XT (i - k)] = 0;
        llist = c_tab.tab->strings;
      while (llist)
	{
	  for (i = c_tab.x2; llist->line[X2XT (i)]; i++)
	    llist->line[X2XT (i - k)] = llist->line[X2XT (i)];
	  llist->line[X2XT (i - k)] = 0;
	  llist = llist->next;
	}			/* while (llist) */
      if (c_tab.del)
	{
	  llist = c_tab.del;
	  for (i = c_tab.x2; llist->line[X2XT (i)]; i++)
	    llist->line[X2XT (i - k)] = llist->line[X2XT (i)];
	  llist->line[X2XT (i - k)] = 0;
	}			/* if (c_tab.del) */
      c_tab.x2 -= k;
      x2t -= k;
    }				/* squeezeField */

    static void spreadLine (char *s, int n, int k)
    {
      /* return the s where after s[n] insert k spaces */
      int i = strlen (s);
      for (; i > n; i--)
	  s[i + k] = s[i];
      for (i = k; i; i--)
	  s[n + i] = ' ';
    }				/* spreadLine */

    static void spreadField (int k)
    {
      /* sread the currebt field of table */
      /* k - how much */
      linelist llist, w, w0, w1;
      int len, i;
      char *s;
        s = c_tab.tab->format;
        len = strlen (s);
      for (i = len; s[i] != SEPR_CHAR; i--)
	{
	};
      len = i;
      s[len + 1] = 0;
      s[X2XT (c_tab.x2) - 1] = ' ';
      spreadLine (s, X2XT (c_tab.x2) - 1, k);
      s[X2XT (c_tab.x2 + k) - 1] = RESIZE_CHAR;

      w0 = 0;
      llist = c_tab.tab->strings;
      while (llist)
	{
	  len = strlen (llist->line);
	  w1 = (linelist) m_alloc (len + k + 1 + 2 * sizeof (linelist));
	  w1->pred = 0;
	  w1->next = 0;
	  strcpy (w1->line, llist->line);
	  spreadLine (w1->line, X2XT (c_tab.x2) - 1, k);
	  if (!w0)
	    w0 = w1;
	  else
	    {
	      w1->pred = w;
	      w->next = w1;
	    }
	  w = w1;
	  if (c_tab.top == llist)
	    c_tab.top = w;
	  if (c_tab.edit == llist)
	    c_tab.edit = w;
	  w1 = llist;
	  llist = llist->next;
	  free (w1);
	}			/* while (llist) */
      c_tab.tab->strings = w0;

      if (c_tab.del)
	{
	  len = strlen (c_tab.del->line);
	  w1 = (linelist) m_alloc (len + k + 1 + 2 * sizeof (linelist));
	  w1->pred = 0;
	  w1->next = 0;
	  strcpy (w1->line, c_tab.del->line);
	  spreadLine (w1->line, X2XT (c_tab.x2) - 1, k);
	  free (c_tab.del);
	  c_tab.del = w1;
	}			/* if (c_tab.del) */
      c_tab.x2 += k;
      x2t += k;
    }				/* spreadField */

    static void minmaxFSize (int *min, int *max)
    {
      /* MIN,MAX - valid size of the current field */
      linelist llist = c_tab.tab->strings;
      char *s = c_tab.tab->format;
      int i;
       *max = c_tab.x2 - c_tab.x1 + STRSIZ - strlen (s) - 1;

       *min = 0;
      for (i = c_tab.x2 - 1; i >= c_tab.x1 && s[c_tab.x0 + i - X1 - 2] == ' '; i--)
	{
	};
      *min = MAX (*min, i - c_tab.x1 + 2);

      while (llist)
	{
	  s = llist->line;
	  for (i = c_tab.x2; i >= c_tab.x1 && s[c_tab.x0 + i - X1 - 2] == ' '; i--)
	    {
	    };
	  *min = MAX (*min, i - c_tab.x1 + 1);
	  llist = llist->next;
	}			/* while (llist) */
    }				/* minmaxFSize */

    static int newFSize (int kk, int min, int max)
    {
      /* return how much to resize field */
      /* kk - current size */
      int key;
      char ss[100];
        key = kk;
        sprintf (ss, " MIN=%d, MAX=%d, current size: ", min, max);
      if (correctInt (X1 + 1, Y1, ss, &kk, 1) && kk >= min && kk <= max)
	  return key = kk;
        return key;
    }				/* newFSize */

    static int resizeField (void)
    {
      int k, kk = c_tab.x2 - c_tab.x1 + 1;
      int min, max, rc;
        minmaxFSize (&min, &max);
        k = newFSize (kk, min, max);
        rc = 1;
      if (k < kk)
	  squeezeField (kk - k);
      else if (k > kk)
	  spreadField (k - kk);
      else
	  rc = 0;
        return rc;
    }				/* resizeField */

/******************* end resize ***********/

    static void shiftL (int dx, table * tab)
    {
      if (!dx)
	return;
      dx = MIN (dx, c_tab.x0 - 1);
      dx = MAX (dx, -x2t + (X2 - X1) + c_tab.x0 - 3);

      if (dx)
	{
	  c_tab.x0 -= dx;
	  c_tab.x1 += dx;
	  c_tab.x2 += dx;
	  scrCur += dx;
	  redrawTable (tab);
	}
    }


    static char fLine[16] = "               ";	/* fine line */
    int edittable (int xx, int yy, table * tab, int ntop, char *hlpf, int show)
    {
      int key;
      void *pscr;
      int i, l;
      char ws[STRSIZ];
      int edit = FALSE;
      int dx, dy;		/*  Mouse */
      int rc;

        clearTypeAhead ();
        initEdit (tab, xx, yy, ntop);
        e_box.x1 = X1;
        e_box.y1 = Y1;
        e_box.x2 = X2;
        e_box.y2 = Y2;
        get_text (X1, Y1 - 1, X2, Y2, &pscr);

        scrcolor (Black, White);
        chepbox (X1, Y1, X2, Y2);
        scrcolor (White, Blue);
        goto_xy (X1, Y1 - 1);
        print ("%*s", X2 - X1 + 1, "");
        l = strlen (tab->headln);
      {
	int z = X1 + (X2 - X1 - 1 - l) / 2;
	if (z < X1 + 18)
	    z = X1 + 18;
	if (z < X2 - 6)
	  {
	    goto_xy (z, Y1 - 1);
	    print ("%.*s", X2 - 6 - z, tab->headln);
	  }
      }
      scrcolor (Red, Blue);
      goto_xy (X1, Y1 - 1);
      print ("%s", "*");
      /* up_bottom bar */
      hbb_new (&up_bar, &up_items[0], boxFrame[1]);
      hbb_place (&up_bar, 1 + X1, Y1);

      /* down_bottom bar */
      hbb_new (&down_bar, &down_items[0], boxFrame[5]);
      hbb_place (&down_bar, 1 + X1, Y2);

      /* init list of hbb's */
      hbl = hbl_add (0, &down_bar);
      if (!show)
	hbl = hbl_add (hbl, &up_bar);
      hbl_colors (hbl, Black, BGmain, Black, White);
      hbl_open (hbl);

      redrawTable (tab);
      key = 0;
      while (key != KB_ESC)
	{
	  if (c_tab.top)
	    {
	      strcpy (ws, c_tab.edit->line + c_tab.x0 + c_tab.x1 - X1 - 2);
	      ws[c_tab.x2 - c_tab.x1 + 1] = 0;
	      if (s_edit (MAX (X1 + 1, c_tab.x1), MIN (X2 - 1, c_tab.x2),
			  MAX (0, X1 - c_tab.x1 + 1),
			  c_tab.y, ws, show))
		{
		  edit = TRUE;
		  for (i = 0; ws[i]; i++)
		    c_tab.edit->line[c_tab.x0 + c_tab.x1 + i - X1 - 2] = ws[i];
		}

	    }
	  key = pop_key ();

	  if (key == KB_MOUSE && mouse_info.but1 == MRELEASED)
	    {
#if EDITTAB_DEBUG
	      printf ("mouse: but1=%d, x=%d, y=%d\n", mouse_info.but1, mouse_info.row, mouse_info.col);
	      printf (" cursor: x=%d, x1=%d, y=%d\n", c_tab.x1, c_tab.x2, c_tab.y);
#endif

	      if (mouse_info.col == X1 && mouse_info.row == Y1 - 1)
		push_key (KB_ESC);
	      else if ((rc = hbl_select (hbl)))
		push_key (rc);
	      else if (c_tab.top && m_inbox (X1, Y1 + 1, X2, Y2)
		       && tab->format[mouse_info.col - X1 + c_tab.x0 - 2] != SEPR_CHAR)
		{
		  linelist ll = c_tab.edit;
		  int dy = mouse_info.row - c_tab.y;

		  if (dy > 0)
		    for (i = dy; i && ll; i--)
		      ll = ll->next;
		  if (dy < 0)
		    for (i = -dy; i && ll; i--)
		      ll = ll->pred;
		  if (ll)
		    {
		      c_tab.edit = ll;
		      c_tab.y += dy;
		      scrCur = mouse_info.col;
		      setFieldFrame (scrCur, tab->format);
		    }
		}
	    }

	  if (c_tab.top == NULL && key != KB_ENTER && key != KB_F1 && key != KB_F2
	      && key != KB_DOWN && key != KB_ESC)
	    key = 0;
	  switch (key)
	    {
	    case KB_PAGEU:	/*  PgUp    */
	      for (i = Y1 + 2; i < Y2; i++)
		{
		  if (c_tab.top->pred != NULL)
		    {
		      c_tab.top = c_tab.top->pred;
		      c_tab.y0--;
		      c_tab.edit = c_tab.edit->pred;
		    }
		}
	      redrawTable (tab);
	      break;

	    case KB_PAGED:	/*  PgDn    */
	      for (i = Y1 + 2; i < Y2; i++)
		{
		  if (c_tab.edit->next != NULL)
		    {
		      c_tab.top = c_tab.top->next;
		      c_tab.y0++;
		      c_tab.edit = c_tab.edit->next;
		    }
		}
	      redrawTable (tab);
	      break;

	    case KB_UP:
	      dy = pop_key ();
	      if (c_tab.y <= (Y1 + 1 + dy))
		{
		  for (; dy && c_tab.top->pred != NULL; dy--)
		    {
		      c_tab.top = c_tab.top->pred;
		      c_tab.y0--;
		      c_tab.edit = c_tab.edit->pred;
		    }
		  redrawTable (tab);
		}
	      if (c_tab.y > (Y1 + 1 + dy))
		shiftY (-dy);
	      break;

	    case KB_ENTER:
	      if (show)
		break;
	      if (!c_tab.del)
		break;
	      {
		linelist del2 = (linelist) m_alloc (sizeof (linerec) + 2);
		strcpy (del2->line, c_tab.del->line);
		if (!c_tab.edit)
		  {
		    c_tab.top = c_tab.del;
		    c_tab.edit = c_tab.del;
		    c_tab.tab->strings = c_tab.top;
		    c_tab.y0 = 1;
		  }
		else
		  {
		    c_tab.del->next = c_tab.edit->next;
		    if (c_tab.edit->next)
		      c_tab.edit->next->pred = c_tab.del;
		    c_tab.del->pred = c_tab.edit;
		    c_tab.edit->next = c_tab.del;
		  }		/* if */
		c_tab.del = del2;
	      }
	      redrawTable (tab);
	      edit = TRUE;
	      push_key (KB_DOWN);
	      push_key (1);
	      break;

	    case KB_DOWN:	/*  LineDn  */
	      dy = pop_key ();
	      if (!c_tab.top && !show)
		push_key (KB_ENTER);
	      else
		{
		  if (c_tab.y >= (Y2 - dy))
		    {
		      for (; dy && c_tab.edit->next != NULL; dy--)
			{
			  c_tab.top = c_tab.top->next;
			  c_tab.y0++;
			  c_tab.edit = c_tab.edit->next;
			}
		      redrawTable (tab);
		    }
		  if (c_tab.y < (Y2 - dy))
		    shiftY (dy);
		  key = 0;
		}
	      break;

	    case KB_TAB:
	      for (; TABCUR < x2t && tab->format[TABCUR] != SEPR_CHAR; scrCur++);
	      scrCur++;
	      if (TABCUR >= x2t)
		scrCur -= TABCUR;
	      setFieldFrame (scrCur, tab->format);
	      if (scrCur <= X1 || scrCur >= X2)
		shiftL (X1 - scrCur + 1, tab);
	      break;

	    case KB_RIGHT:
	      dx = pop_key ();	/* suppose  dx==1 */
	      if (TABCUR >= x2t)
		scrCur--;
	      else
		{
		  if (tab->format[TABCUR] == SEPR_CHAR)
		    {
		      scrCur++;
		      setFieldFrame (scrCur, tab->format);
		    }
		  if (scrCur >= X2)
		    shiftL (-20, tab);
		}
	      break;

	    case KB_LEFT:
	      dx = pop_key ();	/* suppose  dx==1 */
	      if (TABCUR < 0)
		scrCur++;
	      else
		{
		  if (tab->format[TABCUR] == SEPR_CHAR)
		    {
		      scrCur--;
		      setFieldFrame (scrCur, tab->format);
		    }
		  if (scrCur <= X1)
		    shiftL (20, tab);
		}
	      break;

	    case CTRL ('B'):	/* Bottom    */
	      for (i = Y1 + 2; c_tab.edit->next != NULL; i++)
		{
		  c_tab.top = c_tab.top->next;
		  c_tab.y0++;
		  c_tab.edit = c_tab.edit->next;
		}
	      redrawTable (tab);
	      break;

	    case CTRL ('D'):
	      if (show)
		{
		  key = 0;
		  break;
		}
	      if (c_tab.del)
		free (c_tab.del);
	      c_tab.del = c_tab.edit;
	      if (c_tab.top == c_tab.edit)
		{
		  if (c_tab.top->next)
		    {
		      c_tab.top = c_tab.top->next;
		      c_tab.top->pred = c_tab.del->pred;
		      if (c_tab.del->pred)
			c_tab.del->pred->next = c_tab.top;
		      else
			c_tab.tab->strings = c_tab.top;
		    }
		  else
		    {
		      c_tab.top = c_tab.top->pred;
		      if (c_tab.top)
			{
			  c_tab.top->next = 0;
			  c_tab.y0--;
			}
		      else
			{
			  c_tab.tab->strings = c_tab.top;
			  c_tab.y0 = 0;
			}
		    }
		  c_tab.edit = c_tab.top;
		}
	      else
		{
		  if (c_tab.edit->next)
		    {
		      c_tab.edit = c_tab.edit->next;
		      c_tab.edit->pred = c_tab.del->pred;
		      c_tab.del->pred->next = c_tab.edit;
		    }
		  else
		    {
		      c_tab.edit = c_tab.edit->pred;
		      c_tab.edit->next = 0;
		      c_tab.y--;
		    }
		}
	      c_tab.del->next = NULL;
	      c_tab.del->pred = NULL;
	      redrawTable (tab);
	      edit = TRUE;
	      break;

	    case CTRL ('F'):
	      {
		static linelist lastLine = NULL;
		static int lastPos = -1;
		linelist nn = c_tab.edit;
		if ((lastLine == c_tab.edit && lastPos == TABCUR)
		    || correctStr (X1 + 1, Y2, " Find:", fLine, 15, 1))
		  {
		    int dy = 0;
		    char *ip = strstr (nn->line + TABCUR, fLine);
		    if (!ip)
		      {
			nn = nn->next;
			dy++;
			for (; nn && !(ip = strstr (nn->line, fLine)); nn = nn->next)
			  dy++;
		      }
		    if (ip)
		      {
			scrCur += ((ip - nn->line) - TABCUR + strlen (fLine));
			if (TABCUR >= x2t)
			  scrCur--;
			setFieldFrame (scrCur, tab->format);
			if (scrCur <= X1 || scrCur >= X2)
			  shiftL (X1 + strlen (fLine) + 1 - scrCur, tab);
			if (dy)
			  {
			    push_key (KB_DOWN);
			    push_key (dy);
			  }
			lastLine = nn;
			lastPos = TABCUR;
		      }
		    else
		      {
			lastLine = NULL;
			lastPos = -1;
			warnanykey (10, 5, "String not found");
		      }
		  }
	      }
	      break;

	    case CTRL ('G'):
	      {
		int nn = 1;
		if (correctInt (X1 + 1, Y2, "  Go to: ", &nn, 1))
		  {
		    dy = nn - (c_tab.y0 + (c_tab.y - Y1 - 2));
		    if (dy > 0)
		      {
			push_key (KB_DOWN);
			push_key (dy);
		      }
		    else if (dy < 0)
		      {
			push_key (KB_UP);
			push_key (-dy);
		      }
		  }
	      }
	      break;

	    case CTRL ('E'):
	      if (errorText[0])
		warnanykey (10, 10, errorText);
	      else
		be_be ();
	      break;

	    case CTRL ('T'):	/*  Top    */
	      push_key (KB_UP);
	      push_key (c_tab.y0 + c_tab.y - Y1 - 3);
	      break;

	    case CTRL ('S'):
	      if (show)
		break;
	      if (c_tab.tab->format[X2XT (c_tab.x2) - 1] == RESIZE_CHAR
		  && resizeField ())
		{
		  edit = 1;
		  redrawTable (tab);
		}
	      break;

	    case CTRL ('Z'):
	      strcpy (ws, &c_tab.edit->line[c_tab.x0 + c_tab.x1 - X1 - 2]);
	      ws[c_tab.x2 - c_tab.x1 + 1] = 0;
	      hbl_close (hbl);
	      if (zoomStr (ws, 3, 3, 1, show))
		{
		  edit = 1;
		  for (i = 0; ws[i]; i++)
		    c_tab.edit->line[c_tab.x0 + c_tab.x1 + i - X1 - 2] = ws[i];

		};
	      hbl_open (hbl);
	      break;

	    case KB_F1:
	      show_help (hlpf);
	      break;
	    case KB_F2:
	      show_help ("edittab");
	      break;

	    }

	  scrcolor (White, Black);
	  goto_xy (X2 - 6, Y1 - 1);
	  print ("%-4d", c_tab.y0 + (c_tab.y - Y1 - 2));
	  goto_xy (X2, Y1);
#if  EDITTAB_DEBUG
	  printf ("Edittab: key=%d\n", key);
#endif
	}			/* while */

      if (c_tab.del)
	free (c_tab.del);
      if (hbl)
	{
	  hbl_close (hbl);
	  hbl_delete (hbl);
	};
      scrcolor (White, Black);
      put_text (&pscr);
      return edit;
    }
