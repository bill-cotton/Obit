/* $Id: menu.h 2 2008-06-10 15:32:27Z bill.cotton $ */
/* headers for ObitMess menu routines  */
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

#ifndef TMESSMENU_H
#define TMESSMENU_H
/* create main menu */
Widget MakeMainTMessMenu (Widget mainWindow);

/* callback prototypes */
void TMessQuitCB (Widget w, XtPointer clientData, XtPointer callData);

#endif /* TMESSMENU_H */
