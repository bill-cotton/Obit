/* $Id$ */
/* function prototypes for graph.c */
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
#ifndef GRAPH_H
#define GRAPH_H

void GraphDrawLine (ImageData *ImgData, olong blc[2], olong trc[2], olong gcolor);
void GraphDrawCircle (ImageData *ImgData, olong center[2], olong radius, olong gcolor);
void GraphClear (ImageData *ImgData);
void ResetGphCB (Widget w, XtPointer clientData, XtPointer callData);
void GphCSetCB (Widget parent, XtPointer clientData, XtPointer callData);

#endif  /* GRAPH_H */
