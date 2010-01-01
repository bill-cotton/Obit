/* $Id:  $ */
/* Obit task message window routines  */
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

#include <Xm/Xm.h> 
#include "taskmessage.h"
#include "taskmessagebox.h"
#include "obitmess.h"
#include <stdio.h> 

/**
 *  \file taskmessage.c
 * displays task messages.
 */

/*---------------File Structures---------------------------*/
/**  TaskWindow list element structure. */
typedef struct {
  /** Task ID */
  glong taskID;
  /* Task name */
  gchar *taskname;
  /** Task Message Window. */
  TaskMessageBox *taskWindow;
} TaskWindowListElem;

/**  TaskWindow List structure. */
typedef struct {
  /** Number of elements in list */
  glong nTaskWin;
  /** Highest current task ID */
  glong highID;
  /** glib singly linked list */
  GSList* list;
} TaskWindowList;

/*--------------- file global data ----------------*/
TaskWindowList *windowList;

/*---------------Private function prototypes----------------*/
/** TaskWindowListElem constructor */
TaskWindowListElem*  newTaskWindowListElem (glong taskID, 
					    gchar *taskname, 
					    Widget parent);

/** TaskWindowListElem destructor */
TaskWindowListElem*  freeTaskWindowListElem (TaskWindowListElem *in);

/* Add item to task list */
glong TaskWindowListAdd (TaskWindowList *in, gchar *taskname, 
			 Widget parent);

/* Remove  item from task list */
void TaskWindowListRemove  (TaskWindowList *in, glong taskID);

/* Lookup task window */
TaskWindowListElem* TaskWindowListFind  (TaskWindowList *in, 
					 glong taskID);

/* shutdown task message window, arg not used */
void TaskMessageDismiss(XtPointer arg);

/* save content of task message window, arg not used */
void TaskMessageSaveAs(XtPointer arg);

/*-----------------public functions--------------*/

/**
 * Initialize Task Message Windowing
 * \param taskname  name of task, used as label
 * \param status    [out] status, 0=> OK, else failed
 * \return task ID to be used in future calls
 */
void TaskMessageInit(void)
{
  windowList = g_malloc0(sizeof(TaskWindowList));
  windowList->nTaskWin = 0;
  windowList->highID   = 0;
  windowList->list     = NULL;
} /* end TaskMessageInit */

/**
 * Create Task Message Window
 * \param taskname  name of task, used as label
 * \param status    [out] status, 0=> OK, else failed
 * \return task ID to be used in future calls
 */
glong TaskMessageCreate(Widget parent, gchar *taskname, glong *status)
{
  glong taskID=-1;
  TaskWindowListElem* elem;

  *status = 0;

  /* Create widget and add to list */
  taskID = TaskWindowListAdd (windowList, taskname, parent);

  /* Get widget */
  elem = TaskWindowListFind (windowList, taskID);
  if (elem==NULL) {  /* Check */
    *status = 1;
    return -1;
  }

  return taskID;
} /* end  TaskMessageCreate */

/**
 * Display message in task message window
 * \param taskID    ID of task, returned from Create
 * \param message   message to display
 * \param status    [out] status, 0=> OK, else failed
 */
void TaskMessageDisplay(glong taskID, gchar *message, glong *status)
{
  TaskWindowListElem* tmb;

  *status = 0;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL) {  /* Check */
    *status = 1;
    return;
  }

  /* Update display */
  TaskMessageBoxDisplay (tmb->taskWindow, message);

} /* end TaskMessageDisplay */

/**
 * Request user input for task taskID
 * \param taskID    ID of task, returned from Create
 * \param status    [out] status, 0=> OK, else failed
 */
void TaskMessageUserResponseSetup(glong taskID, glong *status)
{
  TaskWindowListElem* tmb;

  *status = 0;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL) {  /* Check */
    *status = 1;
    return;
  }

  /*Setup for response */
  TaskMessageBoxUserRequest (tmb->taskWindow);

} /* end TaskMessageUserResponseSetup */

/**
 * Return user input for task taskID
 * \param taskID    ID of task, returned from Create
 * \param response  user response string
 * \param status    [out] status, 0=> OK, else failed
 */
void TaskMessageUserResponseReply(glong taskID, gchar *response)
{
  TaskWindowListElem* tmb;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL)   /* Check */
    return;

  /* Send reply */
  XMLRPCTWSetReturn((void*)response);
    
} /* end TaskMessageUserResponseReply */

/**
 * Set task status for task taskID
 * \param taskID     ID of task, returned from Create
 * \param taskstatus New status label strioing
 * \param status    [out] status, 0=> OK, else failed
 */
void TaskMessageStatus(glong taskID, gchar *taskstatus, glong *status)
{
  TaskWindowListElem* tmb;

  *status = 0;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL) {  /* Check */
    *status = 1;
    return;
  }

  /* Update widget */
  TaskMessageBoxStatus (tmb->taskWindow, taskstatus);

} /* end  TaskMessageStatus */

/**
 * Check if abort requested
 * \param taskID    ID of task, returned from Create
 * \param status    [out] status, 0=> OK, else failed
 * \return TRUE iff abort requested
 */
gboolean TaskMessagedoAbort(glong taskID, glong *status)
{
  gboolean out = FALSE;
  TaskWindowListElem* tmb;

  *status = 0;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL) {  /* Check */
    *status = 1;
    return out;
  }

  /* Requested? */
  out = tmb->taskWindow->doAbort;

  return out;
} /* end TaskMessagedoAbort */

/**
 * Delete task message window for task taskID
 * \param taskID    ID of task, returned from Create
 * \param status    [out] status, 0=> OK, else failed
 */
void TaskMessageDelete(glong taskID, glong *status)
{
  TaskWindowListElem* tmb;

  *status = 0;

  /* Get widget */
  tmb = TaskWindowListFind (windowList, taskID);
  if (tmb==NULL) {  /* Check */
    *status = 1;
    return;
  }

  /* Delete object & dialog */
  tmb = freeTaskWindowListElem(tmb); 
} /* end TaskMessageDelete */

/*---------------Private functions-------------------------*/
/** 
 * TaskWindowListElem Constructor.
 * \param taskID   Task ID
 * \param taskname Name of task
 * \param parent   Parent widget
 * \return the new object.
 */
TaskWindowListElem*  newTaskWindowListElem (glong taskID, 
					    gchar *taskname, 
					    Widget parent)
{
  TaskWindowListElem* me;

  /* allocate structure */
  me = g_malloc0(sizeof(TaskWindowListElem));

  /* initialize values */
  me->taskID = taskID;
  me->taskname = g_strdup(taskname);
  me->taskWindow = TaskMessageBoxMake (parent, taskname, taskID);

  return me;
} /* end newTaskWindowListElem */

/**
 * TaskWindowListElem destructor.
 * \param in Object to delete
 * \return NULL.
 */
TaskWindowListElem* freeTaskWindowListElem (TaskWindowListElem *in)
{
  g_assert (in != NULL);

  /* Delete scroll text */
  TaskMessageBoxDelete(in->taskWindow);

  /* Task name */
  if (in->taskname) g_free(in->taskname);

  /* deallocate structure */
  g_free (in);
  
  return NULL;
} /* end freeTaskWindowListElem */

/**
 * Add an item to the task window list stack
 * \param in       Pointer to object to add task window to.
 * \param taskname Name of task
 * \param parent   parent widget
 * \return newly assigned task ID;
 */
glong TaskWindowListAdd (TaskWindowList *in, gchar *taskname, 
			 Widget parent)
{
  glong taskID;
  TaskWindowListElem* elem;

  /* assign task ID */
  taskID = in->highID+1;
  in->highID = taskID;

  /* create TaskWindowListElem with message */
  elem = newTaskWindowListElem(taskID, taskname, parent);

  /* add to list */
  in->list = g_slist_append (in->list, elem); 
  in->nTaskWin++;   /* add one to the count */

  return taskID;
} /* end TaskWindowListAdd */

/**
 * Remove task window for taskID from the list
 * \param in       Input TaskWindowList.
 * \param taskID   Task ID to remove from list
 */
void TaskWindowListRemove  (TaskWindowList *in, glong taskID)
{
  TaskWindowListElem* elem;

  /* anything there */
  if (in->nTaskWin <1) return;

   /* retrieve from stack */
  elem =  TaskWindowListFind (in, taskID);
  if (elem==NULL) return; /* nothing in the stack */

  /* remove from list */
  in->list = g_slist_remove(in->list, elem);
  in->nTaskWin--;

  /* deallocate TaskWindowListElem */
  freeTaskWindowListElem(elem);

} /* end TaskWindowListRemove */

/**
 * Lookup element in Task list
 * \param in       Input TaskWindowList.
 * \param taskID   Task ID to find
 */
TaskWindowListElem* TaskWindowListFind (TaskWindowList *in, 
					glong taskID)
{
  GSList *tmp;
  TaskWindowListElem* elem;

  /* anything there? */
  if (in->nTaskWin <1) return NULL;

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (TaskWindowListElem*)tmp->data;
    /* check if this is a match */
    if (elem->taskID==taskID) return elem;
    tmp = g_slist_next(tmp);
  }

  return NULL; /* didn't find */
} /* end TaskWindowListFind */

