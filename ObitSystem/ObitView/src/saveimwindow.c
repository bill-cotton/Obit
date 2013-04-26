/* $Id$  */
/*  Save image window to text file for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 2013
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
#include "saveimwindow.h"
#include "menu.h"
#include "imagedisp.h"
#include "messagebox.h"

/*--------------- file global data ----------------*/
  TextFilePtr   ImWinText = NULL; /* logging TextFile */

/*---------------Private function prototypes----------------*/
/* initialize log file, argument not used */
void SaveImWindowInit (XPointer arg);
/* file selection canceled, arg not used */
void SaveImWindowCancel (XPointer arg);

/*----------------------Public functions---------------------------*/

/**
 * Callback for toggle position logging on/off 
 * \param w           widget activated
 * \param clientData  ImageDisplay
 * \param callData    call data
 */
void SaveImWindowCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata= (ImageDisplay *)clientData;

  /* Startup */
  if (ImWinText) TextFileKill (ImWinText); /* shouldn't need this */
  if (!log_dir) log_dir = MakeFStrng("");
  /* ask for file name and initialize file */
  ImWinText = TextFileMake (IDdata->shell, NULL, log_dir->sp);
  TextFileFind (2, ImWinText, SaveImWindowInit, SaveImWindowCancel);
} /* end SaveImWindowCB */

/**
 * Write current image window to text file
 * \param arg   ignored
 */
void SaveImWindowInit (XPointer arg)
{
  char line[100];
  olong field=1, Id, nId, *win;
  ImageData  *Image = &image[CurImag];
  ObitDConCleanWindowType type;
  ObitDConCleanWindow *edtWindow=Image->edtWindow;

  /* Check that window exists */
  if (edtWindow==NULL) {
    MessageShow ("No Image window given");
    TextFileKill(ImWinText);
    ImWinText = NULL;
    return;
   }
  
  /* Open/create Text File */
  if (TextFileOpen (ImWinText, 2) != 1) {
    MessageShow ("Error opening logging file");
    TextFileKill(ImWinText);
    ImWinText = NULL;
    return;
  } /* end of error trap */
  
  /* Add first entry */
  sprintf (line, "CLEANBox=[ \\");
  if (!TextFileWrite(ImWinText, line)) {
    MessageShow ("Error writing logging file");
    TextFileClose(ImWinText);
    TextFileKill(ImWinText);
    ImWinText = NULL;
    return;
  } /* end of error trap */

  /* Loop */
  nId = edtWindow->maxId[field-1];
  for (Id=1; Id<=nId; Id++) {
    if (ObitDConCleanWindowInfo(edtWindow, field, Id, &type, &win, err)) {
      if (type==OBIT_DConCleanWindow_rectangle)    { /* rectangle */
	sprintf (line, "    %d, %d, %d, %d, \\", win[0], win[1], win[2], win[3]);
      } else if (type==OBIT_DConCleanWindow_round) { /* circle */
 	sprintf (line, "    -1, %d, %d, %d, \\", win[0], win[1], win[2]);
     } else if (type==OBIT_DConCleanWindow_unrectangle) { /* unrectangle */
	sprintf (line, "    %d, %d, %d, %d, # unrectangle \\", win[0], win[1], win[2], win[3]);
      } else if (type==OBIT_DConCleanWindow_unround)     {  /* uncircle */
 	sprintf (line, "    -1, %d, %d, %d, # uncircle \\", win[0], win[1], win[2]);
      } else { /* unknown */
	sprintf (line, "Unknown window type %d", type);
      }
       if (!TextFileWrite(ImWinText, line)) {
	MessageShow ("Error writing logging file");
	TextFileClose(ImWinText);
	TextFileKill(ImWinText);
	ImWinText = NULL; if (win) g_free(win);
	return;
      } /* end of error trap */
    } /* end if window exists */
  } /* end loop over boxes */

  /* Add last entry */
  sprintf (line, "    ]");
  if (!TextFileWrite(ImWinText, line)) {
    MessageShow ("Error writing logging file");
    TextFileClose(ImWinText);
    TextFileKill(ImWinText);
    ImWinText = NULL;
    return;
  } /* end of error trap */
  TextFileClose(ImWinText); /* close file */
  
} /* end SaveImWindowInit */

/**
 * File selection canceled
 * \param arg   ignored
 */
void SaveImWindowCancel (XPointer arg)
{
  /* shut down, bail out */
  /*if (ImWinText&&ImWinText->State) TextFileClose(ImWinText);  close file if open*/
  if (ImWinText) TextFileKill(ImWinText); /* delete structure */
  ImWinText = NULL;
} /* end SaveImWindowCancel */
