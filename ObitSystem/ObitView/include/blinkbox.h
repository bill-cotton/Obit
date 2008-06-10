/* $Id: blinkbox.h,v 1.1.1.1 2005/06/09 12:45:29 bcotton Exp $ */
/* function prototypes for blinkbox.c */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996
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
#ifndef BLINKBOX_H
#define BLINKBOX_H

void BlinkSwapCB (Widget parent, XtPointer clientData, XtPointer callData);
void BlinkBlinkCB (Widget parent, XtPointer clientData, XtPointer callData);
#endif
