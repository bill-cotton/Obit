/* $Id$         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2010                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <sys/types.h>
#include <string.h>
#include "ObitErr.h"
#include "ObitMem.h"
#include "ObitSystem.h"

/** name of the class defined in this file */
static gchar *myClassName = "ObitErr";

/** Has the default message log handler been set? */
static gboolean defaultHandler = FALSE;

/**
 * \file ObitErr.c
 * ObitErr Error stack class function definitions.
 * 
 * This is an error stack class for obtaining tracebacks for error
 * conditions.
 * When an error is detected, it should be entered onto the ObitErr
 * and the function returns.  
 * If an error level of OBIT_Traceback or higher has been entered
 * the error member is set TRUE.
 * Each function with an ObitErr argument should check at the beginning 
 * to see if an error condition already exists (error=TRUE) and if so
 * return.
 * Any function calling a function which encounters an error should add
 * its message to the stack and return.
 *
 *    This class is a member of the Obit class and therefore cannot be 
 * derived from it.  Also is not structured to be derived from.
 *
 * \section ObitErrUsage Usage of member pointers.
 * The Ref and Unref member functions should always be used to make a 
 * copy of an object pointer or to release it.
 * The ref function increments its reference count and returns a pointer.
 * The unref function decrements the reference count, deleted the object
 * if the value is below 1 and returns NULL.
 * Unreferenced pointers should always be NULLed or set to another 
 * valid value.
 */

/*---------------Private function prototypes----------------*/
/** Private: Create ObitErrElem. */
static ObitErrElem* newObitErrElem (ObitErrCode errLevel, gchar *errMsg);

/** Private: Destroy ObitErrElem. */
static ObitErrElem* freeObitErrElem (ObitErrElem *in);

/** Private: Destructor. */
static ObitErr* freeObitErr (ObitErr *in);

/** Private: default log handler */
static void DefaultLogHandler (const gchar *log_domain,
			       GLogLevelFlags log_level,
			       const gchar *message,
			       gpointer user_data);

/*---------------Public functions---------------------------*/
/* constructor */
/**
 * ObitErr Constructor.
 * \return the new object.
 */
ObitErr* newObitErr (void)
{
  ObitErr* me;

  /* allocate structure */
  me = ObitMemAlloc0Name(sizeof(ObitErr),"ObitErr");

  /* initialize values */
  me->className = g_strdup(myClassName);
  me->number    = 0;
  me->error     = FALSE;
  me->stack     = g_queue_new();
  me->ReferenceCount = 1;
  me->prtLv     = 0;
  me->pgmName   = NULL;

  /* Set default handler if not done yet */
  if (!defaultHandler) {
    defaultHandler = TRUE;
    g_log_set_handler (NULL,  G_LOG_LEVEL_WARNING | G_LOG_FLAG_FATAL |
		       G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_MESSAGE |
		       G_LOG_LEVEL_INFO | G_LOG_LEVEL_DEBUG |
		       G_LOG_FLAG_RECURSION, 
		       (GLogFunc)DefaultLogHandler, (gpointer)me);
  }


  return me;
} /* end newObitErr */

/**
 * Unconditional ObitErr destructor.
 * \param in Object to delete
 * \return NULL pointer.
 */
ObitErr* freeObitErr (ObitErr *in)
{
  /* error checks */
  g_assert (ObitErrIsA(in));

  /* clear stack */
  ObitErrClear(in);
  g_queue_free (in->stack); /* destroy queue */

  /* deallocate */
  g_free (in->className);
  ObitMemFree (in);

  return NULL;
} /* end freeObitErr */


/**
 * To reference pointer, incrementing ReferenceCount and returning 
 * the pointer.
 * This function should always be used to copy pointers as this 
 * will ensure a proper reference count.
 * \param in Pointer to object to link.
 * \return the pointer to in.
 */
ObitErr* ObitErrRef (ObitErr* in)
{
  /* error checks */
  g_assert (ObitErrIsA(in));

  /* increment reference count */
  in->ReferenceCount++;

  return in;
} /* end ObitErrRef */

/**
 * Always use this function to dismiss an object as it will
 * ensure that the object is only deleted when there are no more 
 * pointers to it.
 * \param  in Pointer to object to unlink.
 * \return NULL pointer.
 */
ObitErr* ObitErrUnref (ObitErr* in)
{
  if (in==NULL) return NULL;
  /* error checks */
  g_assert (ObitErrIsA(in));

  /* decrement reference count, delete if non positive */
  in->ReferenceCount--;
  if (in->ReferenceCount<=0) freeObitErr(in);

  return NULL;
} /* end ObitErrUnref */

/**
 * Removes all entries in stack and deallocate them
 * \param  in Pointer to object to clear.
 * \return NULL pointer.
 */
void ObitErrClear (ObitErr* in)
{
  gpointer next;

  /* error checks */
  g_assert (in != NULL);

  /* loop through list */
  next = g_queue_pop_head (in->stack);
  while (next!=NULL) {
    freeObitErrElem((ObitErrElem*)next); /* delete */
    next = g_queue_pop_head (in->stack);
  }

  in->number = 0;    /* reset counter */
  in->error = FALSE; /* clear error status */
} /* end ObitErrClear */

/**
 * Removes error entries in stack and clears error status
 * \param  in Pointer to object to clear.
 * \return NULL pointer.
 */
void ObitErrClearErr (ObitErr* in)
{
  gpointer next;
  GQueue *stack = NULL;

  /* error checks */
  g_assert (in != NULL);

  /* Create new, replacement Queue */
  stack = g_queue_new();

  /* loop through list */
  next = g_queue_pop_tail (in->stack);
  while (next!=NULL) {
    /* Is this an error message? */
    if (((ObitErrElem*)next)->errLevel >= OBIT_Traceback) {
      /* error - drop */
      freeObitErrElem((ObitErrElem*)next); /* delete */
      in->number--;        /* decrement counter */
    } else { /* message copy to replacement queue */
      g_queue_push_head (stack, next);
    }
    next = g_queue_pop_tail (in->stack);
  } /* end loop over old queue */

  g_queue_free (in->stack); /* destroy old queue */
  in->stack = stack;        /* replace queue */

  in->error = FALSE; /* clear error status */
} /* end ObitErrClearErr */

/**
 * Add an item to the stack
 * \param in       Pointer to object to add message to.
 * \param errLevel Error level code.
 * \param errMsg   Error message.
 */
void ObitErrPush (ObitErr *in, ObitErrCode errLevel, gchar *errMsg)
{
  ObitErrElem* elem;

  /* error checks */
  g_assert (ObitErrIsA(in));
  g_assert (errMsg != NULL);

  /* create ObitErrElem with message */
  elem = newObitErrElem(errLevel, errMsg);

  /* add to stack */
  g_queue_push_head (in->stack, (gpointer)elem);
  in->number++;   /* add one to the count */
  /* Is this an error or info?*/
  in->error = in->error || (errLevel >= OBIT_Traceback);
} /* end ObitErrPush */

/**
 * Fetch the top item on the stack.
 * This item is then removed.
 * \param in         Input ObitErr.
 * \param errLevel   (output) Error level code.
 * \param errMsg     (output) Error message (deallocate with g_free when done).
 *                   NULL if there are no more messages.
 * \param errTimeTag (output) Error message (deallocate with g_free when done).
 */
void ObitErrPop  (ObitErr *in, ObitErrCode *errLevel, gchar **errMsg, 
		  time_t *errTimeTag)
{
  ObitErrElem* elem;

  /* error checks */
  g_assert (ObitErrIsA(in));

  /* default output */
   *errLevel   = OBIT_None;
   *errMsg     = NULL;
   *errTimeTag = 0;

  /* anything there */
  if (in->number <1) return;

  /* one fewer  when done */
  in->number--;

  /* retrieve from stack */
  elem = (ObitErrElem*)g_queue_pop_tail (in->stack);
  if (elem==NULL) return; /* nothing in the stack */

  /* set output */
  *errLevel   = elem->errLevel;
  *errMsg     = g_strdup(elem->errMsg);
  *errTimeTag = elem->timeTag;

  /* deallocate ObitErrElem */
  freeObitErrElem(elem);

} /* end ObitErrPop */

/**
 * Writes all entries to the g_log.
 * Stack will be cleared when done
 * \param in       Input ObitErr.
 *                 NULL if there are no more messages.
 */
void ObitErrLog  (ObitErr *in)
{
  ObitErrCode errLevel;
  gchar *errMsg, *errLevelStr;
  time_t timeTag;
  struct tm *lp;
/*
 * Human readable versions of the ObitErr codes.
 * Should be coordinated with enum definition.
 */
gchar *ObitErrorLevelString[] = {
  "no msg ",    /* OBIT_None        */
  "info   ",    /* OBIT_InfoErr     */
  "warn   ",    /* OBIT_InfoWarn    */
  "trace  ",    /* OBIT_Traceback   */
  "MildErr",    /* OBIT_MildError   */
  "Error  ",    /* OBIT_Error       */
  "Serious",    /* OBIT_StrongError */
  "Fatal  ",    /* OBIT_Fatal       */
};

  /* error checks */
  g_assert (ObitErrIsA(in));

  /* Make sure program name set */
  if ((in->pgmName==NULL) && (ObitSystemIsInit()))
    in->pgmName = ObitSystemGetPgmName();

  /* loop logging messages */
  ObitErrPop (in, &errLevel, &errMsg, &timeTag);
  while (errMsg!=NULL) {
    /* convert error level to something human readable */
    errLevelStr = ObitErrorLevelString[errLevel];
    /* Time */
    lp = localtime (&timeTag);
    if (lp->tm_year<1000) lp->tm_year += 1900;
    g_log (NULL, G_LOG_LEVEL_MESSAGE,
	   "%s %4d%2.2d%2.2dT%2.2d%2.2d%2.2d %s", 
	   errLevelStr, 
	   lp->tm_year, lp->tm_mon+1, lp->tm_mday,
	   lp->tm_hour, lp->tm_min, lp->tm_sec,
	   errMsg);
    if (errMsg) g_free(errMsg); errMsg=NULL;
    ObitErrPop (in, &errLevel, &errMsg, &timeTag);
  }
  if (errMsg) g_free(errMsg); errMsg=NULL;

  /* Clear any error condition */
  in->error = FALSE;
} /* end ObitErrLog */

/**
 * Adds a message with the current data and time
 * \param in       Input ObitErr.
 */
void ObitErrTimeStamp  (ObitErr *in)
{
  struct tm *lp;
  time_t clock;
  olong timea[3], datea[3];

  /* error checks */
  g_assert (ObitErrIsA(in));

  /* date and time info */
 /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to local arrays */
  datea[0] = lp->tm_year;
  if (datea[0]<1000) datea[0] += 1900; /* full year */
  datea[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  datea[2] = lp->tm_mday;
  timea[0] = lp->tm_hour;
  timea[1] = lp->tm_min;
  timea[2] = lp->tm_sec;

  /* Compose  */
  Obit_log_error(in, OBIT_InfoErr, 
		 "Date/Time: %4d-%2.2d-%2.2d  %2.2d:%2.2d:%2.2d",
		 datea[0],datea[1],datea[2],timea[0], timea[1],timea[2]);

} /* end ObitErrTimeStamp */

/**
 * Adds a log message with the current data and time
 * \param in       Input ObitErr.
 * \param errMsg   Message string
 */
void ObitErrTimeLog  (ObitErr *in, gchar *errMsg)
{
  struct tm *lp;
  time_t clock;
  olong timea[3], datea[3];

  /* error checks */
  g_assert (ObitErrIsA(in));

  /* date and time info */
 /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to local arrays */
  datea[0] = lp->tm_year;
  if (datea[0]<1000) datea[0] += 1900; /* full year */
  datea[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  datea[2] = lp->tm_mday;
  timea[0] = lp->tm_hour;
  timea[1] = lp->tm_min;
  timea[2] = lp->tm_sec;

  /* Compose  */
  Obit_log_error(in, OBIT_InfoErr, 
		 "%4d-%2.2d-%2.2d  %2.2d:%2.2d:%2.2d %s",
		 datea[0],datea[1],datea[2],timea[0],timea[1],timea[2],errMsg);

} /* end ObitErrTimeLog */

/**
 * Determines if the input object is a member of this class
 * \param in Pointer to object to test.
 * \return TRUE if member else FALSE.
 */
gboolean ObitErrIsA (ObitErr* in)
{
  gboolean out;

  /* error checks */
  if (in == NULL) return FALSE;
  if (in->className == NULL) return FALSE;

  /* compare class name member */
  out = !strcmp(in->className, myClassName);

  return out;
} /* end ObitErrIsA */

/*---------------Private functions-------------------------*/
/**
 * ObitErrElem Constructor.
 * \param errLevel Error level code.
 * \param errMsg   Error message.
 * \return the new object.
 */
ObitErrElem*  newObitErrElem (ObitErrCode errLevel, gchar *errMsg)
{
  ObitErrElem* me;

  /* error checks */
  g_assert (errMsg != NULL);

  /* allocate structure */
  me = ObitMemAlloc0Name(sizeof(ObitErrElem),"ObitErrElem");

  /* initialize values */
  me->errLevel = errLevel;
  me->errMsg   = g_strdup(errMsg);
  time(&me->timeTag); /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */

  return me;
} /* end newObitErrElem */

/**
 * ObitErrElem destructor.
 * \param in Object to delete
 * \return NULL.
 */
ObitErrElem* freeObitErrElem (ObitErrElem *in)
{
  g_assert (in != NULL);
  /* deallocate messsage string */
  if (in->errMsg) g_free (in->errMsg);

  /* deallocate structure */
  ObitMemFree (in);
  
  return NULL;
} /* end freeObitErrElem */

/**
 * Default log handler
 * \param message message to log
 */
/** Private: default log handler */
static void DefaultLogHandler (const gchar *log_domain,
			       GLogLevelFlags log_level,
			       const gchar *message,
			       gpointer user_data)
{
  ObitErr* me;
  /* Pass program name */
  if (ObitErrIsA((ObitErr*)user_data)) {
    me  = (ObitErr*)user_data;
    if (me->pgmName)
      fprintf (stdout, "%s: %s\n", me->pgmName, message);
    else
    fprintf (stdout, "Obit: %s\n", message);
  } else {
    fprintf (stdout, "Obit: %s\n", message);
  }
} /* end DefaultLogHandler */
  
