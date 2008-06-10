/* $Id$  */
/* ScrollText routines for ObitView */
/* Scrollable boxes displaying text */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,2002-2008
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
#include <Xm/Xm.h> 
#include <Xm/DialogS.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h>
#include <Xm/ScrollBar.h>
#include <Xm/Form.h>
#include <Xm/ScrolledW.h>
#include <Xm/PushB.h>
#include <X11/Intrinsic.h>
#include "scrolltext.h"
#include "messagebox.h"

/**
 *  \file scrolltext.c
 * displays scrolling text dialog.
 */

/*--------------- file global data ----------------*/
#define SCROLLBOX_WIDTH  500  /* Width of scroll box */
#define SCROLLBOX_HEIGHT 420  /* Height of scroll box */
#define MAXCHAR_LINE     132  /* maximum number of characters in a line */

/*---------------Private function prototypes----------------*/
/* resize event */
void STextResizeCB (Widget w, XtPointer clientData, XtPointer callData);
/* scrollbar changed */
void STextScrollCB (Widget w, XtPointer clientData, XtPointer callData);
/* Dismiss button hit */
void STextDismissButCB (Widget w, XtPointer clientData, XtPointer callData);

/*-----------------public functions--------------*/

/**
 * Create ScrollText and fill it with the contents of a TextFile
 * \param TFilePtr Text file to display a TextFilePtr , is destroyed.
 */
void ScrollTextCopy (XPointer TFilePtr)
{
  ScrollTextPtr STextPtr;
  TextFilePtr   TFile;
  int rcode;
  
  /* make it */
  TFile = (TextFilePtr)TFilePtr;
  STextPtr = ScrollTextMake (TFile->w, TFile->FileName);
  
  /* copy text */
  rcode = 0;
  if (STextPtr) rcode = ScrollTextFill(STextPtr, TFile);
  
  /* Error Message */
  if (!rcode) MessageShow ("Error loading text file to scrolling window");
  
  /* delete TextFile */
  TextFileKill(TFile);
  
  /* delete ScrollText if something went wrong */
  if (!rcode) ScrollTextKill(STextPtr);
} /* end ScrollTextCopy */

/**
 * Create/initialize ScrollText structures
 * \param parent     parent widget
 * \param Title      Title for new dialog
 */
ScrollTextPtr ScrollTextMake (Widget Parent, char* Title)
{
  ScrollTextPtr STextPtr;
  int loop;
  Dimension text_width, text_height, scrollbar_width, butt_height;
  Widget form, DismissButton;
  XFontStruct *XFont;
  int value, increment, slider_size,page_increment;
  
  /* allocate */
  STextPtr = (ScrollTextPtr)g_malloc (sizeof(ScrollTextInfo));
  if (!STextPtr) return NULL;
  
  /* initialize */
  STextPtr->Parent = Parent;
  STextPtr->DismissProc = NULL;
  STextPtr->num_lines = 0;
  STextPtr->first = 0;
  STextPtr->number = 0;
  STextPtr->max_lines = 0;
  STextPtr->Title = (char*)g_malloc(strlen(Title)+1);
  strcpy (STextPtr->Title, Title);
  for (loop=0; loop<MAX_LINE; loop++) STextPtr->lines[loop] = NULL;
  
  /* create main widget */
  STextPtr->ScrollBox = 
    XtVaCreatePopupShell (STextPtr->Title,
			  xmDialogShellWidgetClass, 
			  STextPtr->Parent,
			  XmNautoUnmanage, False,
			  XmNwidth,  (Dimension)SCROLLBOX_WIDTH,
			  XmNheight, (Dimension)SCROLLBOX_HEIGHT,
			  XmNdeleteResponse, XmDESTROY,
			  NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("ScrollTextForm", xmFormWidgetClass,
				  STextPtr->ScrollBox,
				  XmNautoUnmanage, False,
				  XmNwidth,  (Dimension)SCROLLBOX_WIDTH,
				  XmNheight, (Dimension)SCROLLBOX_HEIGHT,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* Play button */
  DismissButton = 
    XtVaCreateManagedWidget (" Dismiss ", 
			     xmPushButtonWidgetClass, 
			     form, 
			     XmNbottomAttachment, XmATTACH_FORM,
			     XmNrightAttachment, XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     NULL);
  XtAddCallback (DismissButton, XmNactivateCallback, STextDismissButCB, 
		 (XtPointer)STextPtr);
  
  /* drawing area */
  scrollbar_width = 15;
  text_width = SCROLLBOX_WIDTH - scrollbar_width;
  XtVaGetValues (DismissButton, XmNheight, &butt_height, NULL);
  text_height = SCROLLBOX_HEIGHT - butt_height;
  STextPtr->TextDraw_wid = text_width;
  STextPtr->TextDraw_hei = text_height;
  value = 1;
  slider_size = 20;
  increment = 1;
  page_increment = 20;
  /* plane scroll */
  STextPtr->TextScrollBar = 
    XtVaCreateManagedWidget ("TextScrollBar", 
			     xmScrollBarWidgetClass, 
			     form,
			     XmNheight,       text_height,
			     XmNwidth,        scrollbar_width,
			     XmNmaximum,         200,
			     XmNminimum,           1,
			     XmNvalue,             1,
			     XmNshowValue,       True,
			     XmNorientation,   XmVERTICAL,
			     XmNprocessingDirection, XmMAX_ON_BOTTOM,
			     XmNrightAttachment,  XmATTACH_FORM,
			     XmNtopAttachment, XmATTACH_FORM,
			     XmNbottomAttachment,  XmATTACH_WIDGET,
			     XmNbottomWidget,      DismissButton,
			     NULL);
  XmScrollBarSetValues (STextPtr->TextScrollBar, value, slider_size, 
			increment, page_increment, False);
  XtAddCallback(STextPtr->TextScrollBar, XmNvalueChangedCallback, 
                STextScrollCB, (XtPointer)STextPtr);
  STextPtr->max_scroll = 200; /* maximum scroll value */
  
  STextPtr->TextDraw = 
    XtVaCreateManagedWidget ("textdraw", xmDrawingAreaWidgetClass, 
			     form, 
			     XmNwidth,              text_width,
			     XmNheight,             text_height,
			     XmNtopAttachment,   XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     XmNrightAttachment,  XmATTACH_WIDGET,
			     XmNrightWidget,  STextPtr->TextScrollBar,
			     XmNbottomAttachment,  XmATTACH_WIDGET,
			     XmNbottomWidget,      DismissButton,
			     NULL); 
  /*   Add callbacks to handle exposures,resize.  */
  XtAddCallback (STextPtr->TextDraw, XmNexposeCallback, STextExposeCB, 
		 (XtPointer)STextPtr);
  XtAddCallback (STextPtr->TextDraw, XmNresizeCallback, STextResizeCB,
		 (XtPointer)STextPtr);
  
  /* set it up */
  XtManageChild (STextPtr->ScrollBox);
  
  /* create graphics context for box */
  STextPtr->gc = XCreateGC (XtDisplay (STextPtr->ScrollBox), 
			    DefaultRootWindow (XtDisplay(STextPtr->ScrollBox)),
			    0, NULL );
  /* how tall are the characters */
  XFont = XQueryFont(XtDisplay (STextPtr->ScrollBox), 
		     XGContextFromGC(STextPtr->gc));
  STextPtr->baseskip = XFont->ascent + XFont->descent + 2;
  
  /* figure out how many lines will fit */
  STextPtr->max_lines = ((int)text_height) / STextPtr->baseskip;
  
  
  return STextPtr; /* return structure */
} /* end ScrollTextMake */

/**
 * Destroy ScrollText structures
 * \param STextPtr  Scrolling text dialog to delete
 * \return 1 if OK
 */
int ScrollTextKill( ScrollTextPtr STextPtr)
{
  int loop;
  
  if (!STextPtr) return 0; /* anybody home? */
  
  /* free up text strings */
  for (loop=0; loop<MAX_LINE; loop++) 
    {if (STextPtr->lines[loop]) {g_free (STextPtr->lines[loop]);
    STextPtr->lines[loop]=NULL;}}
  if (STextPtr->Title) g_free (STextPtr->Title); STextPtr->Title=NULL;
  
  /* release graphics context */
  if (STextPtr->gc) XtReleaseGC (STextPtr->ScrollBox, STextPtr->gc);
  
  /* kill da wabbit */
  XtDestroyWidget (STextPtr->ScrollBox);
  if (STextPtr) g_free(STextPtr); STextPtr=NULL;/* done with this */
  
  return 1;
} /* end ScrollTextKill */

/**
 * Copy text from TextFile to ScrollText 
 * \param STextPtr  Scrolling text dialog to write
 * \param TFilePtr  Text to insert
 * \return 1 if OK
 */
int ScrollTextFill (ScrollTextPtr STextPtr, TextFilePtr TFilePtr)
{
  int loop, rcode, ccode, HitEof, length, number, maxchar=MAXCHAR_LINE;
  char line[MAXCHAR_LINE+1];
  
  if (!STextPtr) return 0; /* anybody home? */
  
  for (loop=0; loop<=MAXCHAR_LINE; loop++) line[loop] = 0; /* zero fill */
  
  rcode = TextFileOpen (TFilePtr, 1); /* open */
  if (rcode!=1) return 0;
  HitEof = 0; number = 0;
  while (!HitEof && (number<MAX_LINE))
    {rcode = TextFileRead (TFilePtr, line, maxchar); /* next line */
    if (rcode==0) break;
    /* swallow line */
    length = strlen (line);
    STextPtr->lines[number] = (char*)g_malloc (length+1);
    strcpy (STextPtr->lines[number], line);
    HitEof = rcode == -1; /* end of file */
    number++;             /* count entries */
    } /* end of loop reading text file */
  STextPtr->num_lines = number;
  
  ccode = TextFileClose (TFilePtr); /* close */
  if ((ccode!=1) || (rcode==0)) 
    {MessageShow ("Error closing Text/FITS file ");
    return 0;} /* error */
  /* final setup */
  ScrollTextInit (STextPtr);
  return 1;
} /* end ScrollTextFill */

/**
 * Initialization of ScrollText after text strings loaded 
 * \param STextPtr  Scrolling text dialog 
 */
void ScrollTextInit (ScrollTextPtr STextPtr)
{
  Dimension cwid, chei;
  int number;
  int value, increment, slider_size,page_increment;
  
  if (!STextPtr) return; /* anybody home? */
  
  /* find size */
  XtVaGetValues (STextPtr->TextDraw, /* get new size */
		 XmNwidth,  &cwid,
		 XmNheight, &chei,
		 NULL);
  STextPtr->TextDraw_hei = chei;
  STextPtr->TextDraw_wid = cwid;
  
  /* number of lines shown*/
  STextPtr->max_lines = (int)chei / STextPtr->baseskip;
  STextPtr->number = STextPtr->max_lines;
  
  /* set up for the expose callback to draw */
  STextPtr->first = 1;
  STextPtr->number = STextPtr->max_lines;
  if (STextPtr->number > STextPtr->num_lines) 
    STextPtr->number = STextPtr->num_lines;
  
  /* need scroll bars? */
  if (STextPtr->num_lines > STextPtr->max_lines)
    /* set scroll bar */
    {number = STextPtr->num_lines + 5;
    if (number<2) number = 2;
    STextPtr->max_scroll = number; /* maximum scroll value */
    if (STextPtr->first>STextPtr->max_scroll) 
      STextPtr->first = STextPtr->max_scroll;
    slider_size = STextPtr->number; 
    value = STextPtr->first;
    if (value>number-slider_size) value = number-slider_size;
    STextPtr->first = value;
    increment = 1;
    page_increment = STextPtr->max_lines-1; 
    if (page_increment<1) page_increment = 1;
    XtVaSetValues(STextPtr->TextScrollBar, 
		  XmNsliderSize, slider_size,
		  XmNminimum,           1,
		  XmNmaximum,       number,
		  NULL);
    XmScrollBarSetValues (STextPtr->TextScrollBar, value, slider_size, 
			  increment, page_increment, False);
    XtMapWidget (STextPtr->TextScrollBar);}
  else
    XtUnmapWidget (STextPtr->TextScrollBar);
} /* end ScrollTextInit */

/**
 * Move scrolling box to the bottom.
 * \param STextPtr  Scrolling text dialog 
 */
void ScrollTextBottom (ScrollTextPtr STextPtr)
{
  Dimension slider_size=0, slider_max=0;
  if (!STextPtr) return; /* anybody home? */
  
  /* find slider size */
  XtVaGetValues (STextPtr->TextScrollBar,
		 XmNsliderSize, &slider_size,
		 XmNmaximum,    &slider_max,
		 NULL);

 /* need scroll bars? */
  if (STextPtr->num_lines > STextPtr->max_lines) {
    XmScrollBarSetValues (STextPtr->TextScrollBar,
			  slider_max-slider_size, 0, 0, 0, True);
    
    STextPtr->first = STextPtr->num_lines-STextPtr->max_lines+1;
  }
} /* end ScrollTextBottom */

/* internal functions */
/**
 * Callback for expose event
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextExposeCB (Widget w, XtPointer clientData, XtPointer callData)
{
  int loop, start, end, x, y, inc;
  ScrollTextPtr STextPtr = (ScrollTextPtr)clientData;
  Display *display;
  Drawable draw;
  GC       gc;
  
  if (!STextPtr) return; /* anybody home? */
  display = XtDisplay(STextPtr->ScrollBox); /* local copies of variables */
  draw = (Drawable)XtWindow(STextPtr->TextDraw);
  gc = STextPtr->gc;
  
  /* blank it out first, draw in white on black background*/
  XSetForeground (display, gc, 
		  BlackPixelOfScreen(XtScreen(STextPtr->TextDraw)));
  XFillRectangle (display, draw, 
		  gc, 0, 0, STextPtr->TextDraw_wid, 
		  STextPtr->TextDraw_hei); 
  XSetForeground (display, gc, 
		  WhitePixelOfScreen(XtScreen(STextPtr->TextDraw)));
  
  start = STextPtr->first-1;
  end = start + STextPtr->number - 1;
  if (end>=STextPtr->num_lines) end = STextPtr->num_lines - 1;
  inc = STextPtr->baseskip;
  x = 2; y = inc;
  for (loop=start; loop<=end; loop++)
    {XDrawString (display, draw, gc, x, y, STextPtr->lines[loop], 
		  strlen(STextPtr->lines[loop]));
    y += inc;
    } /* end loop loadling rows */
} /* end STextExposeCB */

/**
 * Callback for scrollbar changed
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextScrollCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ScrollTextPtr STextPtr = (ScrollTextPtr)clientData;
  XmScrollBarCallbackStruct *call_data = 
    (XmScrollBarCallbackStruct *)callData;
  
  /* read value of scrollbar */
  STextPtr->first = call_data->value; /* 0 rel */
  
  /* redraw */
  STextExposeCB (w, clientData, callData);
} /* end STextScrollCB */

/**
 * Callback for Dismiss button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextDismissButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ScrollTextPtr STextPtr = (ScrollTextPtr)clientData;
  if (!STextPtr) return; /* anybody home? */
  
  /* call any DismissProc */
  if (STextPtr->DismissProc) STextPtr->DismissProc(clientData);
  
  ScrollTextKill (STextPtr);
} /* end STextDismissButCB */

/**
 * Callback for ScrollText resized
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextResizeCB (Widget w, XtPointer clientData, XtPointer callData)
{
  Dimension cwid, chei;
  int number;
  int value, increment, slider_size,page_increment;
  ScrollTextPtr STextPtr = (ScrollTextPtr)clientData;
  if (!STextPtr) return; /* anybody home? */
  
  /* find new size */
  XtVaGetValues (STextPtr->TextDraw, /* get new size */
		 XmNwidth, &cwid,
		 XmNheight, &chei,
		 NULL);
  STextPtr->TextDraw_hei = chei;
  STextPtr->TextDraw_wid = cwid;
  
  /* new number of lines shown*/
  STextPtr->max_lines = (int)chei / STextPtr->baseskip;
  STextPtr->number = STextPtr->max_lines;
  
  /* need scroll bars? */
  if (STextPtr->num_lines > STextPtr->max_lines)
    /* reset Scroll Bar size, limits */
    {number = STextPtr->num_lines + 5;
    STextPtr->max_scroll = number; /* maximum scroll value */
    if (STextPtr->first>STextPtr->max_scroll) 
      STextPtr->first = STextPtr->max_scroll;
    slider_size = STextPtr->number; 
    value = STextPtr->first;
    if (value>number-slider_size) value = number-slider_size;
    STextPtr->first = value;
    increment = 1;
    page_increment = STextPtr->max_lines-1;
    XtVaSetValues(STextPtr->TextScrollBar, XmNmaximum, number, NULL);
    XmScrollBarSetValues (STextPtr->TextScrollBar, value, slider_size, 
			  increment, page_increment, False);
    
    XtMapWidget (STextPtr->TextScrollBar);}
  else
    XtUnmapWidget (STextPtr->TextScrollBar);
  
  /* redraw */
  STextExposeCB (w, clientData, callData);
} /* end STextResizeCB */


