/* $Id: menu.c,v 1.6 2005/07/23 20:51:37 bcotton Exp $ */
/* menu routines for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1997,2002-2008
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
/* Uncomment the following for a version which keeps the file open dialog */
/*#define KEEP_FILE_DIALOG */
#include <Xm/Xm.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <glib.h> 
#include <stdlib.h>
#include <stdio.h>
#include <Xm/Separator.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/FileSB.h>
#include <Xm/CascadeB.h>
#include <Xm/RowColumn.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include "obitview.h"
#include "imagedisp.h"
#include "infobox.h"
#include "moviebox.h"
#include "blinkbox.h"
#include "optionbox.h"
#include "markpos.h"
#include "lookpos.h"
#include "color.h"
#include "graph.h"
#include "textfile.h"
#include "scrolltext.h"
#include "helpbox.h"
#include "aboutbox.h"
#include "logger.h"
#include "messagebox.h"
#include "Image2Pix.h"
#include "menu.h"

/**
 *  \file menu.c
 * menu routines.
 */

/*--------------- file global data ----------------*/
/* global widget arrays for menu items */
Widget Zoom_w[7] = {0, 0, 0, 0, 0, 0, 0};
Widget File_w[5] = {0, 0, 0, 0,0};

/**
 * Routine to create main menu bar
 * \param mainWindow  main window
 * \param image       doesn't seem to be used
 * \param IDdata      Image display
 * \return menu widget
 */
Widget MakeMainMenu (Widget mainWindow, XtPointer image, XtPointer IDdata)
     /* routine to create main menu bar */
{
  Widget menu, pulldown, w;
  /*ImageDisplay *IDDATA = (ImageDisplay*) IDdata;*/
  
  /* initialize help dialog */
  MakeHelpBox (mainWindow, False); /* no immediate display */
  
  /* add menu items */
  menu = XmCreateMenuBar (mainWindow, "menuBar", NULL, 0);
  
  /* "File" */
  pulldown = XmCreatePulldownMenu (menu, "File", IDdata, 0);
  w = XtVaCreateManagedWidget ("File", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'F',
			       NULL);
  /* "Open" */
  File_w[0] = XtVaCreateManagedWidget ("Open", xmPushButtonWidgetClass,
				       pulldown,
				       XmNmnemonic, 'O',
				       NULL);
  /*    XtSetArg (File_w[0], XmNmnemonic, 'O');*/
  XtAddCallback (File_w[0], XmNactivateCallback, OpenCB, IDdata);
  
  /* "Preview" */
  File_w[1] = XtVaCreateManagedWidget ("Preview", xmPushButtonWidgetClass,
				       pulldown,
				       XmNmnemonic, 'P',
				       NULL);
  XtAddCallback (File_w[1], XmNactivateCallback, PreviewCB, IDdata);
  
  /* "Save as" */
  File_w[2] = XtVaCreateManagedWidget ("Save As", xmPushButtonWidgetClass,
				       pulldown,
				       XmNmnemonic, 'S',
				       NULL);
  XtAddCallback (File_w[2], XmNactivateCallback, SaveAsCB, IDdata);
  
  /* "Source info" */
  File_w[3] = XtVaCreateManagedWidget ("Source Info", xmPushButtonWidgetClass,
				       pulldown,
				       XmNmnemonic, 'I',
				       NULL);
  XtAddCallback (File_w[3], XmNactivateCallback, InfoBoxCB, IDdata);
  /* "Position logging" */
  File_w[4] = XtVaCreateManagedWidget ("Start Position Logging", 
				       xmPushButtonWidgetClass,
				       pulldown,
				       XmNmnemonic, 'L',
				       NULL);
  XtAddCallback (File_w[4], XmNactivateCallback, LoggerCB, IDdata);
  /* help */
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown,
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"File menu");
  /* separator */
  w = XtCreateManagedWidget ("sep", xmSeparatorWidgetClass,
			     pulldown, NULL, 0);
  /* "quit" */
  w = XtVaCreateManagedWidget ("Quit", xmPushButtonWidgetClass,
			       pulldown,
			       XmNmnemonic, 'Q',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, QuitCB, IDdata);
  
  /* "Options" */
  pulldown = XmCreatePulldownMenu (menu, "options", IDdata, 0);
  w = XtVaCreateManagedWidget ("Options", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'O',
			       NULL);
  w = XtVaCreateManagedWidget ("Set Options", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'S',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, OptionBoxCB, IDdata);
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Options menu");
  
  /* "Zoom" */
  pulldown = XmCreatePulldownMenu (menu, "zoom", IDdata, 0); 
  w = XtVaCreateManagedWidget ("Zoom", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'Z',
			       NULL);
  /* "Zoom 25%" */
  Zoom_w[0] = XtVaCreateManagedWidget ("25%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '%',
				       NULL);
  XtAddCallback (Zoom_w[0], XmNactivateCallback, Zoom25CB, IDdata);
  /* "Zoom 50%" */
  Zoom_w[1] = XtVaCreateManagedWidget ("50%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '5',
				       NULL);
  XtAddCallback (Zoom_w[1], XmNactivateCallback, Zoom50CB, IDdata);
  /* "Zoom 100%" */
  Zoom_w[2] = XtVaCreateManagedWidget ("100% <=", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '1',
				       NULL);
  XtAddCallback (Zoom_w[2], XmNactivateCallback, Zoom100CB, IDdata);
  /* "Zoom 200%" */
  Zoom_w[3] = XtVaCreateManagedWidget ("200%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '2',
				       NULL);
  XtAddCallback (Zoom_w[3], XmNactivateCallback, Zoom200CB, IDdata);
  /* "Zoom 400%" */
  Zoom_w[4] = XtVaCreateManagedWidget ("400%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '4',
				       NULL);
  XtAddCallback (Zoom_w[4], XmNactivateCallback, Zoom400CB, IDdata);
  /* "Zoom 800%" */
  Zoom_w[5] = XtVaCreateManagedWidget ("800%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '8',
				       NULL);
  XtAddCallback (Zoom_w[5], XmNactivateCallback, Zoom800CB, IDdata);
  /* "Zoom 1600%" */
  Zoom_w[6] = XtVaCreateManagedWidget ("1600%", xmPushButtonWidgetClass,
				       pulldown, 
				       XmNmnemonic, '6',
				       NULL);
  XtAddCallback (Zoom_w[6], XmNactivateCallback, Zoom1600CB, IDdata);
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Zoom menu");
  
  /* "Position" */
  pulldown = XmCreatePulldownMenu (menu, "position", IDdata, 0);
  w = XtVaCreateManagedWidget ("Position", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'P',
			       NULL);
  /* "Set Equinox" */
  w = XtVaCreateManagedWidget ("Set Equinox", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'S',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, SetEquCB, IDdata);
  /* "Mark Position" */
  w = XtVaCreateManagedWidget ("Mark Position", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'M',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, MarkPosCB, IDdata);
  /* "Lookup Position" */
  w = XtVaCreateManagedWidget ("Lookup Position", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'L',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, LookPosCB, IDdata);
  /* "Lookup Position" */
  w = XtVaCreateManagedWidget ("Fit Position", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'F',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, FitPosCB, IDdata);
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Position menu");
  
  /* "Blink" */
  pulldown = XmCreatePulldownMenu (menu, "blink", IDdata, 0);
  w = XtVaCreateManagedWidget ("Blink", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'B',
			       NULL);
  w = XtVaCreateManagedWidget ("Swap Blink and Current", 
			       xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'S',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, BlinkSwapCB, IDdata);
  w = XtVaCreateManagedWidget ("Blink Images", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'B',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, BlinkBlinkCB, IDdata);
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Blink menu");
  
  
  /* "Movie" */
  pulldown = XmCreatePulldownMenu (menu, "movie", IDdata, 0);
  w = XtVaCreateManagedWidget ("Movie", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'M',
			       NULL);
  w = XtVaCreateManagedWidget ("Movie control", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'M',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, MovieBoxCB, IDdata);
  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Movie menu");
  
  /* "Colorize" */
  pulldown = XmCreatePulldownMenu (menu, "colorize", IDdata, 0);
  w = XtVaCreateManagedWidget ("Colorize", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'C',
			       NULL);
  /* "Color Contour" */
  w = XtVaCreateManagedWidget ("Color Contour", xmPushButtonWidgetClass,
			       pulldown,
			       XmNmnemonic, 'L',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, ColContourCB, IDdata);
  
  /* "Pseudo Flame" */
  w = XtVaCreateManagedWidget ("Pseudo Flame", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'P',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, PhlameCB, IDdata);
  
  /* "Grayscale" */
  w = XtVaCreateManagedWidget ("Grayscale", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'G',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, GrayscaleCB, IDdata);
  
  /* "Reverse colors" */
  w = XtVaCreateManagedWidget ("Reverse colors", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'v',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, ReverseColCB, IDdata) ;
  
  /* "Reset" */
  w = XtVaCreateManagedWidget ("Reset Colors", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'R',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, ResetColCB, IDdata);
  
  
  /* "Clear graphics" */
  w = XtVaCreateManagedWidget ("Clear Graphics", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'Z',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, ResetGphCB, IDdata);
  
  /* "Set Graphics Color" */
  w = XtVaCreateManagedWidget ("Set Graphics Color", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'S',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, GphCSetCB, IDdata);

  w = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"Colorize menu");
  
  /* "Help menu" */
  pulldown = XmCreatePulldownMenu (menu, "help", IDdata, 0);
  w = XtVaCreateManagedWidget ("Help", xmCascadeButtonWidgetClass,
			       menu, XmNsubMenuId, pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  w = XtVaCreateManagedWidget ("About ObitView", xmPushButtonWidgetClass,
			       pulldown, 

			       XmNmnemonic, 'A',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpAboutCB, IDdata);
  w = XtVaCreateManagedWidget ("Help me", xmPushButtonWidgetClass,
			       pulldown, 
			       XmNmnemonic, 'H',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, HelpBoxTopicCB, 
		 (XtPointer)"How to use help");
  
  XtManageChild (menu);
  return (menu);
} /* end of MakeMainMenu */

/* callback functions */
/**
 * Callback for Quit button - bail out
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void QuitCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* turn off logging if necessary */
  if (doLog) LoggerCB (w, clientData, callData);
  
  /* terminate program */
  exit (0);
} /* end QuitCB */

/**
 * Loads selected FITS file into the ImageData attached to the display 
 * resets display when image is read
 * \param filebox     widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void FileOKCB (Widget filebox, XtPointer clientData, XtPointer callData)
{
  gchar *filename, *directory, szErrMess[121];
  XmFileSelectionBoxCallbackStruct *cbs;
  ImageDisplay  *IDdata;
  
  cbs = (XmFileSelectionBoxCallbackStruct *) callData;
  IDdata = (ImageDisplay *)clientData;
  
  /* get file name */
  if (!XmStringGetLtoR (cbs->value, XmSTRING_DEFAULT_CHARSET, &filename))
    return; /* error */
  /* read FITS file to pixmap */
  FStrngFill (image[CurImag].FileName, filename);
  image[CurImag].DataType = OBIT_IO_FITS;    /* FITS */
  image[CurImag].reLoad   = TRUE;            /* Force load */
  /* read FITS file to PixMap */
  if (Image2Pix (&image[CurImag], IDdata, TRUE))  {
    /* error */
    g_snprintf (szErrMess, 120, "Error reading FITS file = %s", filename);
    MessageShow (szErrMess);
  }
  
  /* get directory name */
  if (!XmStringGetLtoR (cbs->dir, XmSTRING_DEFAULT_CHARSET, &directory))
    return; /* error */
  if (!FITS_dir) FITS_dir = MakeFStrng(" ");
  FStrngFill (FITS_dir, directory);
  if (filename) XtFree(filename); filename = NULL; 
  if (directory) XtFree(directory); directory = NULL;
  
  /* Shazam disappear and die */
#ifndef KEEP_FILE_DIALOG  /* option to let it hang around */
  XtUnmanageChild (filebox);
  XtPopdown (XtParent (filebox));
  XtDestroyWidget(filebox); 
#endif
  
  /* reset display */
  ResetDisplay(IDdata);
  
} /* end FileOKCB */

void FileCancelCB (Widget filebox, XtPointer clientData, XtPointer callData)
     /* cancel file selection dialog box */
{
  /* Shazam disappear and die */
  XtUnmanageChild (filebox); 
  XtPopdown (XtParent (filebox)); 
  XtDestroyWidget(filebox); 
} /* end FileCancelCB */

/**
 * Callback for Open: Load new FITS file
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void OpenCB (Widget w, XtPointer clientData, XtPointer callData)
{
  Widget       filebox;
  ImageDisplay *IDdata;
  XmString     wierdstring = NULL;
  Arg          wargs[5]; 
  
  IDdata = (ImageDisplay *)clientData;
  /* sometimes c is too stupid to live */
  filebox = (Widget) XmCreateFileSelectionDialog (w, "file_open", NULL, 0);
  XtAddCallback (filebox, XmNokCallback, FileOKCB, clientData);
  XtAddCallback (filebox, XmNcancelCallback, FileCancelCB, clientData);
  XtAddCallback (filebox, XmNhelpCallback, HelpBoxTopicCB, 
		 (XtPointer)"Browser");
  
  /* set directory if it is defined */
  if (FITS_dir) {
    wierdstring = XmStringCreateSimple (FITS_dir->sp);
    XtSetArg (wargs[0], XmNdirectory, wierdstring);
    XtSetValues (filebox, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }
  /* Shazam appear: */
  XtManageChild (filebox);
  XtPopup (XtParent (filebox), XtGrabNone);
  /* all the action is in the callback routine */
} /* end OpenCB */

/**
 * Callback for Preview: show text or FITS header in scroll box
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void PreviewCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata;
  TextFilePtr TFilePtr;
  
  IDdata = (ImageDisplay *)clientData;
  
  /* make dummy text file */
  TFilePtr = TextFileMake (IDdata->shell, NULL, NULL);
  
  /* set directory if known*/
  if (FITS_dir)
    {if(TFilePtr->directory) g_free(TFilePtr->directory); TFilePtr->directory=NULL;
    TFilePtr->directory = (char*)g_malloc(FITS_dir->length+1);
    strcpy (TFilePtr->directory, FITS_dir->sp);}
  
  /* ask which file and copy to ScrollText */
  TextFileFind (1, TFilePtr, ScrollTextCopy, (TextFileProc)TextFileKill);
} /* end PreviewCB */

/**
 * Set Zoom menu labels to indicate which zoom is on
 * number is the menu item number 0 = 25% to 6=1600% 
 * \param number  which zoom factor
 */
void MenuMarkZoom (int number)
{
  XmString  wierdstring;
  olong       loop;
  gchar      *clear[] = 
    {"25%", "50%", "100%", "200%", "400%", "800%", "1600%"};
  gchar      *set[] = 
    {"25% <=", "50% <=", "100% <=", "200% <=", "400% <=", "800% <=", 
     "1600% <="};
  
  /* sanity checks */
  if ((number<0) || (number>6)) return;
  if (!Zoom_w[0]) return;
  
  /* clear all first */
  for (loop=0; loop<7; loop++)
    {
      wierdstring = XmStringCreateSimple (clear[loop]);
      XtVaSetValues(Zoom_w[loop], 
		    XmNlabelString,  wierdstring,
		    NULL);
      if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
    }
  /* set the specified one */
  wierdstring = XmStringCreateSimple (set[number]);
  XtVaSetValues(Zoom_w[number], 
		XmNlabelString,  wierdstring,
		NULL);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
} /* end MenuMarkZoom */

/**
 * Set position logging menu label to indicate if hitting it will
 * start or stop position logging.
 * \param onoff 1 = turn on, 2 = turn off 
 */
void MenuMarkLogger (int onoff)
{
  XmString  wierdstring;
  
  if (onoff==1) /* turn on */
    wierdstring = XmStringCreateSimple ("Start Position Logging");
  else if (onoff==2) /* turn off */
    wierdstring = XmStringCreateSimple ("Stop Position Logging");
  else /* can't deal with this */
    return;
  
  XtVaSetValues(File_w[4], XmNlabelString,  wierdstring, NULL);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
} /* end MenuMarkLogger */

/**
 * Callback for Copy current displayed FITS file to specified name 
 * \param filebox     widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void CopyOKCB (Widget filebox, XtPointer clientData, XtPointer callData)
{
  gchar *filename;
  XmFileSelectionBoxCallbackStruct *cbs;
  ImageDisplay  *IDdata;
  FStrng *FileName;
  
  cbs = (XmFileSelectionBoxCallbackStruct *) callData;
  IDdata = (ImageDisplay *)clientData;
  
  /* get file name */
  if (!XmStringGetLtoR (cbs->value, XmFONTLIST_DEFAULT_TAG, &filename))
    return; /* error */
  FileName = MakeFStrng (filename);
  /* copy current file to output */
  /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&   Probably not right &&&&&&&&&&&&&&&&&&&&&&&&*/
  FStrngFill (image[CurImag].FileName, filename);
/*    { */
/*       *//* error */
/*      sprintf (szErrMess, "Error copying FITS file %s",  */
/*	       image[CurImag].FileName->sp); */
/*      MessageShow (szErrMess); */
/*      sprintf (szErrMess, "to %s", filename); */
/*      MessageShow (szErrMess); */
 /*   } */
  if (FileName) KillFStrng (FileName); FileName = NULL;
  if (filename) XtFree(filename); filename = NULL;
  
  /* Shazam disappear and die */
  XtUnmanageChild (filebox); 
  XtPopdown (XtParent (filebox)); 
  XtDestroyWidget(filebox); 
  
  /* reset display */
  ResetDisplay(IDdata);
  
} /* end CopyOKCB */

/**
 * Callback for save copy of FITS file
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void SaveAsCB (Widget w, XtPointer clientData, XtPointer callData)
{
  Widget       filebox;
  ImageDisplay *IDdata;
  XmString     wierdstring;
  Arg          wargs[5]; 
  
  IDdata = (ImageDisplay *)clientData;
  
  filebox = (Widget) XmCreateFileSelectionDialog (w, "save_as", NULL, 0);
  XtAddCallback (filebox, XmNokCallback, CopyOKCB, clientData);
  XtAddCallback (filebox, XmNcancelCallback, FileCancelCB, clientData);
  XtAddCallback (filebox, XmNhelpCallback, HelpBoxTopicCB, 
		 (XtPointer)"Browser");
  
  /* set directory if it is defined */
  if (FITS_dir) {
    wierdstring = XmStringLtoRCreate (FITS_dir->sp, 
				      XmSTRING_DEFAULT_CHARSET);
    XtSetArg (wargs[0], XmNdirectory, wierdstring);
    XtSetValues (filebox, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }
  /* Shazam appear: */
  XtManageChild (filebox);
  XtPopup (XtParent (filebox), XtGrabNone);
  /* all the action is in the callback routine */
} /* end SaveAsCB */
