/* $Id$ */
/* headers for XFITSview message box routines  */
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
#include <glib.h>
#ifndef MESSAGEBOX_H
#define MESSAGEBOX_H
/* display message in message box, creating if necessary */
void MessageShow (char *message);

/* redraw message box if it exists */
void MessageRefresh (void);
/* Obit g_log message handler */
void MessageObitHandler (const gchar *log_domain,
			 GLogLevelFlags log_level,
			 const gchar *message,
			 gpointer user_data);

#endif
