/* $Id$ */
/* header file for ScrollText (Unix/Xwindow) routines */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996
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
#include "textfile.h"

#ifndef SCROLLTEXT_H
#define SCROLLTEXT_H

#define MAX_LINE 1024  /* maximum number of lines in ScrollText */

typedef struct {
Widget Parent;           /* parent window */
char*  Title;            /* title of window */
Widget ScrollBox;        /* scrolling text window */
Widget TextScrollBar;    /* scroll bar */
Widget TextDraw;         /* text drawing area */
TextFileProc DismissProc;/* procedure to call when box destructs */
int    TextDraw_wid;     /* width of TextDraw */
int    TextDraw_hei;     /* height of TextDraw */
GC     gc;               /* graphics context for box */
char*  lines[MAX_LINE];  /* pointers to the lines of text to be displayed */
int    num_lines;        /* total number of lines */
int    first;            /* first line in display */
int    number;           /* number of lines in display */
int    max_lines;        /* maximum number of lines that will fit in display */
int    max_scroll;       /* maximum scroll value */
int    baseskip;         /* baseline skip (distance between rows) */
} ScrollTextInfo, *ScrollTextPtr;

/* public function prototypes */
/* create ScrollText and fill it with the contents of a TextFile */
/* argument actually a TextFilePtr - which is destroyed*/
void ScrollTextCopy (XPointer TFilePtr);

/* create/initialize ScrollText structures */
ScrollTextPtr ScrollTextMake (Widget Parent, char* Title);

/* destroy ScrollText structures and associated TextFile */
/* returns 1 if OK */
int ScrollTextKill( ScrollTextPtr STextPtr);

/* copy text from TextFile to ScrollText */
/* returns 1 if OK */
int ScrollTextFill (ScrollTextPtr STextPtr, TextFilePtr TFilePtr);

/* expose event */
/* clientData = ScrollText */
void STextExposeCB (Widget w, XtPointer clientData, XtPointer callData);

/* Initialization of ScrollText after text strings loaded */
void ScrollTextInit (ScrollTextPtr STextPtr);

/* Move scroll box to the bottom */
void ScrollTextBottom (ScrollTextPtr STextPtr);

#endif /* end SCROLLTEXT_H */
