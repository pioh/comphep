/*
* Copyright (C) 2002-2009, CompHEP Collaboration
* ------------------------------------------------------
*/
#ifndef _TAGPARSER_H_
#define _TAGPARSER_H_

/***********************************************************/
/* compute int lenth                                       */
extern int intlen (int src);

/***********************************************************/
/* add extra comand COM to end of tag T,                   */
/* returns 1 if all is normal, 0 otherwise.                */
extern int add_com (string_comnd com, elementary_tag * t);

/***********************************************************/
/* remove comand with name CNAME from tag T,               */
/* return 1 if all is normal,                              */
/* -1 if the tag does not contain the command.             */
extern int remove_com (char *cname, elementary_tag * t);

/***********************************************************/
/* add extra tag T to end of tags base BASE,               */
/* return 1 if all is normal, 0 otherwise.                 */
extern int add_tag (elementary_tag * t, tags * base);

/***********************************************************/
/* remove tag with numner NUM from tags base BASE,         */
/*  return 1 if all is fine, -1 - otherwise                */
extern int remove_tag (char *tname, tags * base);

/***********************************************************/
/* remove tag with name TNAME from tags base BASE,         */
/*  return 1 if all is normal,                             */
/* -1 if the base does not contain the tag.                */
extern int remove_tag_num (int num, tags * base);

/***********************************************************/
/* remove tags with null names in tags base BASE,          */
/*  return 1 if all is fine, -1 oterwise                   */
extern int remove_null_tags (tags * base);

/***********************************************************/
/* get number of tag from tags base BASE with tagname NAME */
/*  return the tag number in the base                      */
/* -1 if the base does not contain the tag.                */
extern int get_tag (int first, tags * base, char *name);

/***********************************************************/
/* get number of tag from tags base BASE with tagname NAME */
/* and which contains command COM with name COM.NAME       */
/*  return the tag number in the base                      */
/*  COM.VALUE contains value form the tag                  */
/* -1 if the base does not contain the tag.                */
extern int
get_tag_with1com (int first, tags * base, char tagname[2048], string_comnd * com);

/***********************************************************/
/* get number of tag from tags base BASE with tagname NAME */
/* and which contains commands COM1 (with name COM1.NAME)  */
/* and  COM2 (with name COM2.NAME)                         */
/*  return the tag number in the base                      */
/*  COM1.VALUE and COM2.VALUE contain values form the tag  */
/* -1 if the base does not contain the tag.                */
extern int get_tag_with2com (int first, tags * base, char *tname,
			     string_comnd * com1, string_comnd * com2);

/***********************************************************/
/* get number of tag from tags base BASE with tagname NAME */
/* and which contains command COM                          */
/*  return the tag number in the base                      */
/* -1 if the base does not contain the tag.                */
extern int get_tag_with_exactcom (int first, tags * base,
				  char *tagname, string_comnd com);

/***********************************************************/
/* get number of tag from tags base BASE with tagname NAME */
/* and which contains commands COM1 and COM2               */
/*  return the tag number in the base                      */
/* -1 if the base does not contain the tag.                */
extern int get_tag_with_exact2com (int first, tags * base,
		   char *tagname, string_comnd * com1, string_comnd * com2);

/***********************************************************/
/* get int value of command with name CNAME from tag T     */
/*  return the the value                                   */
/* -1 if the tag does not contain the command.             */
extern int get_ival (int n, char *cname, elementary_tag * t);

/***********************************************************/
/* get double value of command with name CNAME from tag T  */
/* return the the value,                                   */
/* -1.0 if the tag does not contain the command.           */
extern double get_fval (int n, char *cname, elementary_tag * t);

/***********************************************************/
/* get string value of command with name CNAME from tag T  */
/* return the the value,                                   */
/* empty string if the tag does not contain the command.   */
extern char *get_cval (int n, char *cname, elementary_tag * t);

/***********************************************************/
/* replace command with name COM.NAME in tag T.            */
/* return the tag number in the base,                      */
/* -1 if the tag does not contain the command COM          */
/* (with name COM.NAME).                                   */
extern int replace_com (string_comnd com, elementary_tag * t);

/***********************************************************/
/* change int value of command COM from tag T to NEWVAL    */
/* SHIFT is number of backspaces from COM.NAME= and        */
/* COM.VALUE                                               */
/* return 1 if all is normal,                              */
/* -1 if the tag does not contain the command COM          */
/* (with name COM.NAME).                                   */
extern int change_ival (int newval, int shift, string_comnd com,
			elementary_tag * t);

/***********************************************************/
/* change double value of command COM from tag T to NEWVAL */
/* return 1 if all is normal,                              */
/* -1 if the tag does not contain the command COM          */
/* (with name COM.NAME).                                   */
extern int change_fval (double newval, string_comnd com,
			elementary_tag * t);

/***********************************************************/
/* change string value of command COM from tag T to NEWVAL */
/* return 1 if all is normal,                              */
/* -1 if the tag does not contain the command COM          */
/* (with name COM.NAME).                                   */
extern int change_cval (char *newval, string_comnd com,
			elementary_tag * t);

extern int tag_contain_com (string_comnd * com, elementary_tag * t);
extern int tag_contain_exactcom (string_comnd com, elementary_tag * t);

int event_parser (char *in, char **eventword);

#endif /* parser.h  */
