/* $Id:  $ */
/* Task Message box */
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
#ifndef TASKMESSAGEBOX_H
#define TASKMESSAGEBOX_H
/*---------------Class Structure---------------------------*/
/** TaskMessageBox Class structure. */
typedef struct {
  /** text box to display messages */
  Widget  task_dialog;
  Widget  task_text;
  Widget  user_input;
  Widget  status;
  glong    taskID;
  gboolean doAbort;
  gboolean doUserInput;
  gchar    userString[200];
} TaskMessageBox;

/*---------------Class Functions---------------------------*/
/* TaskMessageBox constructor */
TaskMessageBox* TaskMessageBoxMake (Widget parent, gchar *taskname, 
				    glong taskID); 

/* Add a message to display */
void TaskMessageBoxDisplay (TaskMessageBox* tmb, gchar *message);

/* Setup user input   */
void TaskMessageBoxUserRequest (TaskMessageBox* tmb);

/* Save contents of text box to a file */
void TaskMessageBoxSaveAs(TaskMessageBox* tmb, gchar *filename);

/* Set status */
void TaskMessageBoxStatus(TaskMessageBox* tmb, gchar *status);

/* Destructor*/
void TaskMessageBoxDelete(TaskMessageBox* tmb);

#endif /* end TASKMESSAGEBOX_H */


