/* $Id: color.h,v 1.3 2005/07/20 14:09:18 bcotton Exp $ */
/* function prototypes for color.c */
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
#ifndef COLOR_H
#define COLOR_H

void SetupColorMap (Widget shell, ImageDisplay *IDdata);
void ColContourCB (Widget w, XtPointer clientData, XtPointer callData);
void PhlameCB (Widget w, XtPointer clientData, XtPointer callData);
void GrayscaleCB (Widget w, XtPointer clientData, XtPointer callData);
void ReverseColCB (Widget w, XtPointer clientData, XtPointer callData);
void ResetColCB (Widget w, XtPointer clientData, XtPointer callData);
void SetColorTable (ImageDisplay* IDdata);
void InitColorTable(ImageDisplay* IDdata);
void InitGraphColor(ImageDisplay* IDdata);

#endif
