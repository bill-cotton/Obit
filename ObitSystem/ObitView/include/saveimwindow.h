/* $Id$  */
/* Header file for saveimwindow.c, image window to text file routines */
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
#include "obitview.h"
#include "textfile.h"
#include "ObitDConCleanWindow.h"
#ifndef SAVEIMWINDOW_H
#define SAVEIMWINDOW_H


/* public function prototypes */
/* toggle start/stop logging callback clientData = ImageDisplay */
void SaveImWindowCB (Widget w, XtPointer clientData, XtPointer callData);

#endif /* SAVEIMWINDOW_H */ 
