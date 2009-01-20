/* $Id$ */
/* function prototypes for drawbox.c */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005
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
#include "imagedisp.h"
#include "ObitDConCleanWindow.h"   /* DEBUG */
#ifndef DRAWBOX_H
#define DRAWBOX_H

void DrawBoxInit (ImageDisplay *IDdata);
olong DrawBox (ImageDisplay *IDdata, ObitDConCleanWindow* window, 
	      olong field);
olong EditBox (void);
olong ClearBox (void);

#endif  /* GRAPH_H */
