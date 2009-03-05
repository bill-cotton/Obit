/* $Id: $ */
/* menu routines for ObitMess */
/*-----------------------------------------------------------------------
*  Copyright (C) 2009
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
*  Correspondence concerning this software should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
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
#include "obitmess.h"
#include "helpbox.h"
#include "aboutbox.h"
#include "messagebox.h"
#include "TMessMenu.h"

/**
 *  \file TMessMenu.c
 * menu routines.
 */

/*--------------- file global data ----------------*/
/**
 * Routine to create main menu bar
 * \param mainWindow  main window
 * \param image       doesn't seem to be used
 * \param IDdata      Image display
 * \return menu widget
 */
Widget MakeMainTMessMenu (Widget mainWindow)
     /* routine to create main menu bar */
{
  Widget menu, pulldown, w;
  
  /* initialize help dialog */
  MakeHelpBox (mainWindow, False); /* no immediate display */
  
  /* add menu items */
  menu = XmCreateMenuBar (mainWindow, "menuBar", NULL, 0);
  
  /* "File" */
  pulldown = XmCreatePulldownMenu (menu, "File", NULL, 0);
  w = XtVaCreateManagedWidget ("File", xmCascadeButtonWidgetClass,
                               menu, XmNsubMenuId, pulldown,
                               XmNmnemonic, 'F',
                               NULL);
  /* "quit" */
  w = XtVaCreateManagedWidget ("Quit", xmPushButtonWidgetClass,
			       pulldown,
			       XmNmnemonic, 'Q',
			       NULL);
  XtAddCallback (w, XmNactivateCallback, TMessQuitCB, NULL);
  
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
void TMessQuitCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* terminate program */
  exit (0);
} /* end TMessQuitCB */

