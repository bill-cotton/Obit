/* $Id: FStrng.h,v 1.1 2005/07/23 00:20:08 bcotton Exp $ */
/*  header file for FStrng string class  */ 
/* an FStrng is an implementation of character strings with associated
   functions */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996, 2008
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

#include <string.h> 
#include <glib.h> 
#include "ObitTypes.h"
#ifndef FSTRNG_H 
#define FSTRNG_H 
  
typedef struct FSTRNG { 
    olong length;    /*  Length of string  */ 
    gchar *sp;          /*  String pointer    */ 
  } FStrng; 
  
  
/** Constructors  */ 
FStrng*  MakeFStrng(gchar *str) ; 
FStrng*  MakeFStrngSize(olong length) ; 

/*  Destructor  */ 
void  KillFStrng(FStrng* str); 
  
/* fill with char array */ 
void FStrngFill(FStrng *out, gchar *in); 
  
/* Substring  */ 
/* Return substring given by 0-rel indices first and last  */ 
FStrng* subFStrng(FStrng* this, olong first, olong last); 
  
/* Replace characters first to last with string s.  */ 
void repFstrng(FStrng* this, olong first, olong last, const FStrng *s); 
  
/* Search for character in a string; return index if found; -1 if not.  */ 
olong FStrngFind_char(FStrng* this, gchar ch); 
  
/* Copy FStrng in to FStrng out  */ 
void FStrngCopy(FStrng* out, FStrng *in); 
  
/* Comparison  */ 
gboolean FStrngComp (FStrng *s1, FStrng *s2); 
  
/* Concatenate  */ 
FStrng* FStrngConcat(FStrng *s1, FStrng *s2); 
  
  
#endif /*  FSTRNG_H  */ 
