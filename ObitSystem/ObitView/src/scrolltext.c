/* $Id$  */
/* ScrollText routines for ObitView */
/* Scrollable boxes displaying text */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,2002-2026
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
#include <Xm/Text.h>
#include <X11/Intrinsic.h>
#include "scrolltext.h"
#include "obitview.h"
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
void STextResizeCB(Widget w, XtPointer clientData, XEvent *event, Boolean *continue_to_dispatch);
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
  //printf ("in ScrollTextCopy: before ScrollTextMake file %s\n", TFile->FileName); // debug
  STextPtr = ScrollTextMake (TFile->w, TFile->FileName);
  
  /* copy text */
  //printf ("in ScrollTextCopy: before ScrollTextFill\n"); // debug
  rcode = 0;
  if (STextPtr) rcode = ScrollTextFill(STextPtr, TFile);
  
  //printf ("in ScrollTextCopy: rcode %d\n",rcode); // debug
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
  Widget form, DismissButton;
  
  //printf ("in ScrollTextMake title %s\n", Title); // debug
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
  STextPtr->TextDraw_wid = (int)(SCROLLBOX_WIDTH*sizeFactor);
  STextPtr->TextDraw_hei = (int)(SCROLLBOX_HEIGHT*sizeFactor);
  STextPtr->Title = (char*)g_malloc(strlen(Title)+1);
  strcpy (STextPtr->Title, Title);

  /* create main widget */
  STextPtr->ScrollTop = 
    XtVaCreatePopupShell (STextPtr->Title,
			  xmDialogShellWidgetClass, 
			  STextPtr->Parent,
			  XmNautoUnmanage, False,
			  XmNwidth,  (Dimension)STextPtr->TextDraw_wid,
			  XmNheight, (Dimension)STextPtr->TextDraw_hei, 
			  XmNdeleteResponse, XmDESTROY,
			  XmNfontList,   textFontList, 
			  NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("ScrollTextForm", xmFormWidgetClass,
				  STextPtr->ScrollTop,
				  XmNautoUnmanage, False,
				  XmNwidth,  (Dimension)STextPtr->TextDraw_wid,
				  XmNheight, (Dimension)STextPtr->TextDraw_hei, 
				  XmNx,           0,
				  XmNy,           0,
				  XmNfontList,   textFontList, 
				  NULL);
  XtAddEventHandler(form, StructureNotifyMask, False, STextResizeCB, (XtPointer)STextPtr);
  //XtAddCallback (form, XmNhelpCallback, STextResizeCB, (XtPointer)STextPtr);
  
  /* dismiss button */
  /* Create the Motif compound string for the label */
  XmString dismiss_label = XmStringCreateLocalized("Dismiss");

  DismissButton = 
    XtVaCreateManagedWidget (" Dismiss ", 
			     xmPushButtonWidgetClass, 
			     form,
			     XmNlabelString,     dismiss_label,   /* Explicitly set the visible text */
			     XmNbottomAttachment, XmATTACH_FORM,
			     XmNrightAttachment, XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     XmNfontList,        textFontList, 
			     NULL);
  XtAddCallback (DismissButton, XmNactivateCallback, STextDismissButCB, (XtPointer)STextPtr);
  /* Always free the compound string immediately after the widget is created */
  XmStringFree(dismiss_label);

  STextPtr->num_lines = 31;
  STextPtr->num_cols  = 77;
  // suggested by google:  
  Arg args[15];
  int n = 0;
  // Set resources for the underlying text widget
  XtSetArg(args[n], XmNeditMode, XmMULTI_LINE_EDIT); n++;
  XtSetArg(args[n], XmNrows,     STextPtr->num_lines); n++;
  XtSetArg(args[n], XmNcolumns,  STextPtr->num_cols); n++;
  XtSetArg(args[n], XmNfontList, textFontList); n++;     // set font

  // Google suggests
 // Attach the TOP of the text widget to the TOP of the form
XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM); n++;
XtSetArg(args[n], XmNtopOffset, 10); n++; // 10-pixel gap from top

// Attach the LEFT of the text widget to the LEFT of the form
XtSetArg(args[n], XmNleftAttachment, XmATTACH_FORM); n++;
XtSetArg(args[n], XmNleftOffset, 10); n++; // 10-pixel gap from left

// Attach the RIGHT of the text widget to the RIGHT of the form
XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM); n++;
XtSetArg(args[n], XmNrightOffset, 10); n++; // 10-pixel gap from right

// Attach the Bottom of the text widget to the top  of the dismiss button
XtSetArg(args[n], XmNbottomAttachment, XmATTACH_WIDGET); n++;
XtSetArg(args[n], XmNbottomWidget, DismissButton); n++;
XtSetArg(args[n], XmNbottomOffset, 10); n++; // 10-pixel gap from right

// 3. Create the text widget as a child of the form
STextPtr->TextDraw = XmCreateScrolledText(form, "textDraw", args, n);
// 4. Manage both widgets so they display
XtManageChild(STextPtr->TextDraw);
XtManageChild(form);

  return STextPtr; /* return structure */
} /* end ScrollTextMake */

/**
 * Destroy ScrollText structures
 * \param STextPtr  Scrolling text dialog to delete
 * \return 1 if OK
 */
int ScrollTextKill( ScrollTextPtr STextPtr)
{
  if (!STextPtr) return 0; /* anybody home? */
  
  /* free up text strings */
  if (STextPtr->Title) {g_free (STextPtr->Title);} STextPtr->Title=NULL;
  
  /* kill da wabbit */
  //XtDestroyWidget (STextPtr->ScrollBox);
  XtDestroyWidget (STextPtr->TextDraw);
  //STextPtr->ScrollBox = NULL;
  STextPtr->TextDraw = NULL;
  if (STextPtr) {g_free(STextPtr);} STextPtr=NULL;/* done with this */
  
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
  int loop, rcode, ccode, HitEof, maxchar=MAXCHAR_LINE;
  char line[MAXCHAR_LINE+1];
  
  if (!STextPtr) return 0; /* anybody home? */
  
  for (loop=0; loop<=MAXCHAR_LINE; loop++) line[loop] = 0; /* zero fill */
  
  rcode = TextFileOpen (TFilePtr, 1); /* open */

 // following google suggestions
 XmTextPosition last_pos = XmTextGetLastPosition(STextPtr->TextDraw);
 if (rcode!=1) return 0;
  HitEof = 0; 
  while (!HitEof)
    {rcode = TextFileRead (TFilePtr, line, maxchar); /* next line */
    if (rcode==0) break;
    /* swallow line */
    // 2. Insert the new text at that position
    XmTextInsert(STextPtr->TextDraw, last_pos, line);
    last_pos = XmTextGetLastPosition(STextPtr->TextDraw);
    XmTextInsert(STextPtr->TextDraw, last_pos, "\n");
    last_pos = XmTextGetLastPosition(STextPtr->TextDraw);
    HitEof = rcode == -1; /* end of file */
    } /* end of loop reading text file */
  // 3. Optional: Automatically scroll down to show the new text
  XmTextShowPosition(STextPtr->TextDraw, last_pos);
  // grumble XtManageChild(STextPtr->ScrollBox); // Show it
  //XtManageChild(XtParent(STextPtr->ScrollBox)); // And its parent too
  XtManageChild(XtParent(STextPtr->TextDraw)); // And its parent too

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
  int it[5];
  
  if (!STextPtr) return; /* anybody home? */
  
  /* find size */
  XtVaGetValues (STextPtr->TextDraw, /* get new size */
		 XmNwidth,  &it[0],
		 XmNheight, &it[2],
		 NULL);
  cwid = (Dimension)it[0];
  chei = (Dimension)it[2];
  STextPtr->TextDraw_hei = chei;
  STextPtr->TextDraw_wid = cwid;
  
} /* end ScrollTextInit */

/**
 * Move scrolling box to the bottom - not used but needed for link.
 * \param STextPtr  Scrolling text dialog 
 */
void ScrollTextBottom (ScrollTextPtr STextPtr)
{
  return;
} /* end ScrollTextBottom */

/* internal functions */
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
  
  // Permanently deletes the widget and frees memory
  XtDestroyWidget(STextPtr->ScrollTop);
  STextPtr->ScrollTop = NULL;
  
  ScrollTextKill (STextPtr);
} /* end STextDismissButCB */

/**
 * Callback for expose event (unused but needed to link)
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextExposeCB (Widget w, XtPointer clientData, XtPointer callData)
{
  return;  // not needed
} /* end STextExposeCB */

/**
 * Event handler for ScrollText resized
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void STextResizeCB(Widget w, XtPointer clientData, XEvent *event, Boolean *continue_to_dispatch)
{
  Dimension cwid, chei;
  int it[5];
  ScrollTextPtr STextPtr = (ScrollTextPtr)clientData;

  if (!STextPtr) return; /* anybody home? */

  // Any change?
  /* find new size */
  XtVaGetValues (STextPtr->TextDraw, /* get new size */
		 XmNwidth,  &it[0],
		 XmNheight, &it[2],
		 NULL);
  cwid = (Dimension)it[0];
  chei = (Dimension)it[2];
  //printf ("in STextResizeCB, size %d %d\n", (int)cwid, (int)chei); // debug
  // not needed?if (((int)cwid==STextPtr->TextDraw_wid) && ((int)chei==STextPtr->TextDraw_hei)) return;

  // From google
  short new_columns, new_rows;
  int font_height,  font_width;
  XFontStruct *fontStruct = NULL;
  XmFontContext context;
  XmFontListEntry entry;
  XmStringCharSet charset;
  /* Ensure this is a geometry change event */
  if (event->type == ConfigureNotify) {
    XConfigureEvent *cevent = &event->xconfigure;
    
    /* 1. Extract your ScrolledText widget pointer passed via client_data */
    //Widget scrolledText = STextPtr->ScrollBox;
    Widget scrolledText = STextPtr->TextDraw;
    
    /* 2. Calculate new rows and columns based on pixel sizes */
    /* Initialize the font list context to read the first entry */
    XmFontListInitFontContext(&context, textFontList);
    entry = XmFontListNextEntry(context);
    
    if (entry != NULL) {
        /* Extract the raw X11 font structures from the entry */
        XtPointer font_ptr = XmFontListEntryGetFont(entry, (XmFontType *)&charset);
        fontStruct = (XFontStruct *)font_ptr;
        
        if (fontStruct != NULL) {
            /* Compute exact character height from font baselines */
            font_height = fontStruct->ascent + fontStruct->descent;
            
            /* Get character width (using average or '0' character width) */
            font_width = fontStruct->max_bounds.width; 
            
            /* Use font_width and font_height for your calculations here */
        }
    }
    XmFontListFreeFontContext(context);
    
    
    /* Subtract padding/scrollbar allowances if needed */
    new_columns = (short)(cevent->width / font_width);
    new_rows = (short)(cevent->height / font_height);
    new_rows -= 3;  // don't eat Dismiss button
    
    /* Enforce a safe minimum size to prevent crashes */
    if (new_columns < 5)  new_columns = 5;
    if (new_rows < 2)     new_rows = 2;
    
    /* 3. Apply the new dimensions using XtVaSetValues */
    XtVaSetValues(scrolledText,
		  XmNrows, new_rows,
		  XmNcolumns, new_columns,
		  NULL);
  } else return; // end if reconfigure

  // change - new size
  STextPtr->TextDraw_hei = chei;
  STextPtr->TextDraw_wid = cwid;
  
  /* new number of lines shown*/
  STextPtr->num_lines = new_rows;
  STextPtr->num_cols  = new_columns;
} /* end STextResizeCB */

