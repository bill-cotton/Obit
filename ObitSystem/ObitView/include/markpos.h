/* $Id$ */
/* function prototypes for markpos.c */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2016
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
#ifndef MARKPOS_H
#define MARKPOS_H

void MarkPosInit (ImageDisplay *IDdata);
void MarkPosCB (Widget parent, XtPointer clientData, XtPointer callData);
void MarkPosFromFile (char* filename, ImageDisplay* IDdata);
void MarkPix (ImageDisplay *IDdata, olong iX, olong iY, olong iInner, olong iOuter);
int MarkPosXML (char *pos);
#endif
