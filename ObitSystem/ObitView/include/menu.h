/* $Id$ */
/* headers for XFITSview menu routines  */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2013
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

#ifndef MENU_H
#define MENU_H
/* create main menu */
Widget MakeMainMenu (Widget mainWindow, XtPointer image, XtPointer IDdata);

/* set Zoom menu labels to indicate which zoom is on */
/* number is the menu item number 0 = 25% to 6=1600% */
void MenuMarkZoom (int number);

/* set position logging menu label to indicate if hitting it will */
/* start or stop position logging, 1 = turn on, 2 = turn off */
void MenuMarkLogger (int onoff);

/* callback prototypes */
void QuitCB (Widget w, XtPointer clientData, XtPointer callData);
void OpenCB (Widget w, XtPointer clientData, XtPointer callData);
void OpenACB (Widget w, XtPointer clientData, XtPointer callData);
void PreviewCB (Widget w, XtPointer clientData, XtPointer callData);
void FileCancelCB (Widget filebox, XtPointer clientData, XtPointer callData);
void FileOKCB (Widget filebox, XtPointer clientData, XtPointer callData);
void SaveAsCB (Widget w, XtPointer clientData, XtPointer callData);

#endif
