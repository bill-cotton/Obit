/* $Id: $ */
/* headers for Obit task message window routines  */
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
*  Correspondence concerning ObitView should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#include <glib.h>
#ifndef TASKMESSAGE_H
#define TASKMESSAGE_H
/** Init task message windowing */
void TaskMessageInit(void);
/** Create task message window */
glong TaskMessageCreate(Widget parent, gchar *taskname, glong *status);
/** Display task message */
void TaskMessageDisplay(glong taskID, gchar *message, glong *status);
/** Request user input  */
void TaskMessageUserResponseSetup(glong taskID, glong *status);
/** Request user input  */
void TaskMessageUserResponseReply(glong taskID, gchar *response);
/** Set task status  */
void TaskMessageStatus(glong taskID, gchar *taskstatus, glong *status);
/** Check if abort requested  */
gboolean TaskMessagedoAbort(glong taskID, glong *status);
/** Delete task message window */
void TaskMessageDelete(glong taskID, glong *status);
#endif
