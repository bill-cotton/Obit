/* $Id$  */
/* Header file for logger.c, position logging routines */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2008
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
#ifndef LOGGER_H
#define LOGGER_H


/* public function prototypes */
/* toggle start/stop logging callback clientData = ImageDisplay */
void LoggerCB (Widget w, XtPointer clientData, XtPointer callData);

/* log position */
void DoLogger (void);

#endif /* LOGGER_H */ 
