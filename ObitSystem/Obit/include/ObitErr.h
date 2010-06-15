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
#ifndef OBITERR_H 
#define OBITERR_H 
#include <glib.h>
#include <time.h>
#include "ObitTypes.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitErr.h
 * ObitErr Error stack class definition.
 * 
 * This is an error stack class for obtaining tracebacks for error
 * conditions.
 * When an error is detected, it should be entered onto the ObitErr
 * and the function returns.
 * Each function with an ObitErr argument should check at the beginning 
 * to see if an error condition already exists (number > 0) and if so
 * return.
 * Any function calling a function which encounters an error should add
 * its message to the stack and return.
 * Numerous macroes assist with use of the ObitErr.
 *
 *  A number of macroes are defined to simplify managing the error stack.
 * Implementation uses a glib GQueue.
 * 
 * \section ObitErrUsage Usage
 * Instances can be obtained using the #newObitErr constructor or a 
 * pointer duplicated using the #ObitErrRef function.
 * When an instance is no longer needed, use the #ObitErrUnref function
 * to release it.
 */

/*----------------- enums ---------------------------*/
/**
 * enum for error codes.
 * Should be coordinated with ObitErrorLevelString in ObitErr.c.
 */
enum ObitErrCode {
  /** no message */
  OBIT_None = 0,
  /** information */
  OBIT_InfoErr,
  /**  warning */
  OBIT_InfoWarn,
  /** no new error - use for traceback */
  OBIT_Traceback,
  /** mild error */
  OBIT_MildError,
  /** error */
  OBIT_Error,
  /** serious error */
  OBIT_StrongError,
  /**Fatal condition  */
  OBIT_Fatal
};/* end enum ObitErrCode */
/**
 * typedef for enum for error codes
 */
typedef enum ObitErrCode ObitErrCode;

/** Size of message buffer in characters*/
#define OBITERRBUFSIZE 120

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to log error message
 * adds errMsg/traceback to err for errors
 * \li err     = ObitErr to receive messages.
 * \li errCode = ObitErr level for message
 * \li format  = formatting string and arguments.
 */
#define Obit_log_error(err,errCode,format...)        G_STMT_START{ \
        g_snprintf (err->buffer, OBITERRBUFSIZE, format);          \
	ObitErrPush (err, errCode, err->buffer);                   \
        if (errCode>OBIT_Traceback) {                              \
           g_snprintf (err->buffer, OBITERRBUFSIZE,                \
                       " This occured at file %s: line %d (%s)",   \
                        __FILE__, __LINE__, __PRETTY_FUNCTION__);  \
	   ObitErrPush (err, OBIT_Error, err->buffer);             \
	  }                                                        \
        }G_STMT_END  

/** 
 * Macro to log error message without tracepacb
 * Evaluates expression and if false, adds errMsg/traceback to err
 * \li err     = ObitErr to receive messages.
 * \li errCode = ObitErr level for message
 * \li format  = formatting string and arguments.
 */
#define Obit_log_error_no_trace(err,errCode,format...) G_STMT_START{ \
        g_snprintf (err->buffer, OBITERRBUFSIZE, format);          \
	ObitErrPush (err, errCode, err->buffer);                   \
        }G_STMT_END  

/** 
 * Macro to evaluate expression, log in err and return on failure.
 * Evaluates expression and if false, adds errMsg to err
 * log the location and then return (no return value).
 *\ li expr = expression to evaluate
 * \li err = ObitErr to receive messages.
 * \li format  = formatting string and arguments.
 */
#define Obit_return_if_fail(expr,err,format...)   G_STMT_START{  \
     if (expr) { } else  {                                       \
        g_snprintf (err->buffer, OBITERRBUFSIZE, format);        \
	ObitErrPush (err, OBIT_Error, err->buffer);              \
        g_snprintf (err->buffer, OBITERRBUFSIZE,                 \
                    " Occured at file %s: line %d (%s)",         \
                     __FILE__, __LINE__, __PRETTY_FUNCTION__);   \
	ObitErrPush (err, OBIT_Error, err->buffer);              \
        return;                                                  \
     } }G_STMT_END  

/** 
 * Macro to evaluate expression, log in err and return a value 
 * on failure.
 * Evaluates expression and if false, adds errMsg to err
 * log the location and then return, passing out if nonNULL.
 * \li expr = expression to evaluate
 * \li err = ObitErr to receive messages.
 * \li out = return value.
 * \li format  = formatting string and arguments.
 */
#define Obit_retval_if_fail(expr,err,out,format...) G_STMT_START{ \
     if (expr) { } else  {                                        \
        g_snprintf (err->buffer, OBITERRBUFSIZE, format);         \
	ObitErrPush (err, OBIT_Error, err->buffer);               \
        g_snprintf (err->buffer, OBITERRBUFSIZE,                  \
                    " Occured at file %s: line %d (%s)",          \
                     __FILE__, __LINE__, __PRETTY_FUNCTION__);    \
	ObitErrPush (err, OBIT_Error, err->buffer);               \
        return out;                                               \
     } }G_STMT_END  

/** 
 * Macro for traceback when an error in a called routine is encountered.
 * Writes traceback info and returns (no return value).
 * \li err = ObitErr structure to log to.
 * \li me = name of calling routine;
 * \li name = object name.
 */
#define Obit_traceback_msg(err,me,name) G_STMT_START{           \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                "Traceback: routine %s:  object %s", me, name); \
     ObitErrPush (err, OBIT_Traceback, err->buffer);            \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                 " Occured at file %s: line %d (%s)",           \
                  __FILE__, __LINE__, __PRETTY_FUNCTION__);     \
      ObitErrPush (err, OBIT_Traceback, err->buffer);           \
      return;                                                   \
     }G_STMT_END

/** 
 * Macro for traceback logging with a return value
 * Writes traceback info and returns, optionally passing a value.
 * \li err = ObitErr structure to log to.
 * \li me = name of calling routine;
 * \li name = object name.
 * \li out = return code
 */
#define Obit_traceback_val(err,me,name,out) G_STMT_START{       \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                "Traceback: routine %s:  object %s", me, name); \
     ObitErrPush (err, OBIT_Traceback, err->buffer);            \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                 " at file %s: line %d (%s)",                   \
                  __FILE__, __LINE__, __PRETTY_FUNCTION__);     \
      ObitErrPush (err, OBIT_Traceback, err->buffer);           \
      return out;                                               \
     }G_STMT_END

/** 
 * Macro to dump cfitsio error stack to an ObitErr.
 * Clears cfitsio error stack when done
 * \li err = ObitErr to recieve messages.
 */
#define Obit_cfitsio_error(err)   G_STMT_START{                       \
      while (fits_read_errmsg(err->buffer))                           \
	ObitErrPush (err, OBIT_Error, err->buffer);                   \
      fits_clear_errmsg();     }G_STMT_END  

/*---------------Class Structures---------------------------*/
/**  ObitErr Class structure. */
typedef struct {
  /** class name for verification */
  gchar *className;
  /** Number of error messages. */
  gint16 number;
  /** Has an error been entered */
  gboolean error;
  /** Stack structure. */
  GQueue *stack;
  /** Reference count of pointers to this object. */
  gint32 ReferenceCount;
  /** message buffer */
  gchar buffer[OBITERRBUFSIZE+1];
  /** Print level flag */
  olong prtLv;
  /** Program name */
  gchar *pgmName;
  /** logfile name */
  gchar *logFile;
} ObitErr;

/**  ObitErr Stack Element structure. */
typedef struct {
  /** Error level as enum ObitErrCode */
  ObitErrCode errLevel;
  /** Error message string. */
  gchar *errMsg;
  /** Time tag Unix seconds */
  time_t timeTag;
} ObitErrElem;

/* Private functions defined only in ObitErr.c */
/*---------------Public functions---------------------------*/

/** Public: Constructor. */
ObitErr* newObitErr (void);

/** Public: Reference to object, update reference count. */
ObitErr* ObitErrRef (ObitErr* in);

/** Public: Unreference object, destroy if no more references. */
ObitErr* ObitErrUnref (ObitErr* in);

/** Public: Initialize message/error stack. */
void ObitErrInit (ObitErr* in,  gpointer info);

/** Public: Clear error stack. */
void ObitErrClear (ObitErr* in);

/** Public: Clear only error messages and status in stack. */
void ObitErrClearErr (ObitErr* in);

/** Public: Add entry in error stack. */
void ObitErrPush (ObitErr *in, ObitErrCode errLevel, gchar *errMsg);

/** Public: Pop last entry from top of stack. */
void ObitErrPop  (ObitErr *in, ObitErrCode *errLevel, gchar **errMsg, 
		  time_t *errTimeTag);

/** Public: Write all entries in log file. */
void ObitErrLog  (ObitErr *in);

/** Public: Add timestamp message. */
void ObitErrTimeStamp  (ObitErr *in);

/** Public: Add message with timestamp. */
void ObitErrTimeLog  (ObitErr *in, gchar *errMsg);

/** Public: Returns true if input is a  ObitErr* */
gboolean ObitErrIsA (ObitErr* in);

#endif /* OBITERR_H */ 

