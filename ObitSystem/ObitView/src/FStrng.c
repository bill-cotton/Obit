/* $Id: FStrng.c,v 1.1 2005/07/23 00:20:08 bcotton Exp $ */
/*  Implementation for class FStrng string  */ 
/* an FStrng is an implementation of character strings with associated
   functions */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1998,2002-2008
*  Associated Universities, Inc. Washington DC, USA.
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*-----------------------------------------------------------------------*/
#include <stdlib.h> 
#include <stdio.h> 
#include "FStrng.h" 
#include "messagebox.h" 
    
/**
 * Constructors from gchar*
 * \param str  C string to insert
 * \return new FStrng  KillFStrng when done
 */
FStrng*  MakeFStrng(gchar *str) 
{ 
  FStrng *me; 
  if (!str) return NULL; 
  me = (FStrng *)g_malloc(sizeof(FStrng)); 
  me->length = strlen(str); 
  me->sp = (gchar*) g_malloc0((size_t)me->length+1);
  strcpy (me->sp, str); 
  return me; 
} /* end MakeFStrng */ 
  
/**
 * Constructors from size 
 * \param   length number of characters
 * \return  new FStrng  KillFStrng when done
 */
FStrng*  MakeFStrngSize(olong length) 
{ 
  FStrng *me; 
  me = (FStrng *) g_malloc(sizeof(FStrng)); 
  me->length = length; 
  me->sp = (gchar*) g_malloc0((size_t)me->length+1);
  return me; 
} /* end MakeFStrngSize */ 
  
/**
 * Destructor
 * \param  str FStrng to delete
 */
void  KillFStrng(FStrng* str) 
{ 
  if (!str) return;  /* nobody home */ 
  if (str->sp) g_free(str->sp); str->sp=NULL;/* delete actual string */ 
  if (str) g_free(str); str = NULL;
}  /* end KillFStri\ng */ 
  
/**
 * Fill with gchar array 
 * \param out  FStrng to fill
 * \param in   source string
 */
void FStrngFill(FStrng *out, gchar * in) 
{
  gboolean remake; 
  olong i, length; 

  if (!out) return; 
  if (!in) return; 
  length = strlen(in); 
  remake = (out->sp == 0); 
  remake = remake || (out->length!=length); 
  if (remake)   /* If strings aren't same size remake output  */ 
    {if (out->sp!=0) {g_free(out->sp); out->sp=NULL;} 
    out->sp = (gchar*) g_malloc(length+1); 
    out->length = length;
    } 
  for (i=0; i<length; i++) out->sp[i] = in[i]; out->sp[i] = 0; 
} /* End StrngFill */ 
  
  
  
/* Substrings  */ 
/**
 * Return substring given by 0-rel indices first and last
 * \param  me     FStrng to extract
 * \param  first  first character 
 * \param  last   last character
 * \return  new FStrng  KillFStrng when done
 */
FStrng* subFstrng(FStrng* me, olong first, olong last) 
{ 
  olong j, i; 
  gchar *tout, szErrMess[100]; 
  FStrng *out; 
  if (!me) return NULL;
  if (((first<0) || (last<0))  || ((first>me->length) 
				   || (last>me->length))) {
    sprintf (szErrMess, "substring: Error specifying substring %d  %d %d ", 
	     first, last, me->length); 
    MessageShow(szErrMess); 
    return NULL;
  } 

  tout = (gchar*)g_malloc(last-first+1); 
  j=0;
  for (i=first; i<=last; i++) tout[j++] = me->sp[i]; 
  out = MakeFStrng(tout); 
  if (tout) g_free(tout); tout = NULL;
  return out;
}  /*  End of subFStrng  */ 
  
/**
 * Replace characters first to last with string s.
 * \param  me     FStrng to modify
 * \param  first  first character 
 * \param  last   last character
 * \param  s      Source FStrng
 */
void repFStrng(FStrng* me, olong first, olong last, const FStrng *s) 
{ 
  olong i, j; 
  gchar szErrMess[100];
  if (!me) return;
  if (!s) return;
  
  if ((first<0) || (last<0) || (first>me->length) || (last>me->length) 
      || ((last-first+1)<s->length)) 
    {sprintf (szErrMess, "repFstring: Error specifying substring, %d  %d %d ", 
	      first, last, s->length); 
    MessageShow (szErrMess); 
    return;} 
  j=s->length-1; 
  for (i=first; i<last; i++) me->sp[i] = ' ';  /* blank fill  */ 
  i=last; 
  while (j>=0) me->sp[i--] = s->sp[j--]; /* right justify copy  */ 
} /* End of repFStrng  */ 
  
/**
 * Search for character in a string; return index if found; -1 if not.
 * \param  me     FStrng to search
 * \param  ch     Character to search for
 * \return 0-rel index if found, else -1
 */
olong FStrngFind_char(FStrng* me, gchar ch) 
{ 
  olong i; 
  if (!me) return -1;
  for (i=0; i<me->length; i++) if (me->sp[i]==ch) return i; 
  return -1;
}  /*End of FStrngFind_gchar  */ 
  
/**
 * Copy contents of FStrng in to FStrng out
 * \param  out  destination FStrng
 * \param  in   source FStrng
 */
void FStrngCopy(FStrng* out, FStrng *in) 
{
  gboolean remake; 
  olong i; 
  if (!out) return;
  if (!in) return;

  remake = (out==0) || (out->sp == 0); 
  remake = remake || (out->length!=in->length); 
  if (remake)   /* If strings aren't same size remake output  */ 
    {if (out->sp!=0) g_free(out->sp); out->sp=NULL;
    out->sp = (gchar*) g_malloc(in->length+1); 
    out->length = in->length;} 
  for (i=0; i<in->length; i++) out->sp[i] = in->sp[i]; out->sp[i] = 0; 
}  /* End of FStrngCopy  */ 
  
/**
 * Compare contents of FStrng in to FStrng out
 * \param  out  first FStrng
 * \param  in   second FStrng
 * \return TRUE if the contents are the same else FALSE
 */
/* Comparison  */ 
gboolean FStrngComp (FStrng *s1, FStrng *s2) 
{
  gboolean test; 
  olong i; 
  if (!s1) return FALSE;
  if (!s2) return FALSE;
  test = (s1->length == s2->length); 
  if (!test) return FALSE; 
  for (i=0; i<s1->length; i++) test = (test && (s1->sp[i]==s2->sp[i])); 
  return test; 
}  /*  End FStrngComp  */ 
  
/**
 * Concatenate out = s1+s2
 * \param  s1 first FStrng
 * \param  s2 second FStrng
 * \return  new FStrng  KillFStrng when done
 */
FStrng* FStrngConcat (FStrng *s1, FStrng *s2) 
{ 
  olong i, j; 
  gchar *tout; 
  FStrng *out; 
  if (!s1) return NULL;
  if (!s2) return NULL;

  tout = (gchar*)g_malloc(s1->length+s2->length+1); 
  j = 0; 
  for (i=0; i<s1->length; i++) tout[j++] = s1->sp[i]; 
  for (i=0; i<s2->length; i++) tout[j++] = s2->sp[i]; 
  tout[j] = 0; 
  out = MakeFStrng(tout); 
  if (tout) g_free(tout); tout = NULL;
  return out; 
}  /*  End FStrngConcat  */ 
  
  
  
