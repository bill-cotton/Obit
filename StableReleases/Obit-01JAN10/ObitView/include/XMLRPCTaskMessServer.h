/* $Id:  $  */
/* XMLRPC server for ObitMess */
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
#include <pthread.h>
#include "glib.h"
#include "ObitTypes.h"

#ifndef XMLRPCTASKMESSSERVER_H
#define XMLRPCTASKMESSSERVER_H
/*------------------ Structures -----------------------------*/


/*-------------- enumerations -------------------------------------*/
/**
 * \enum xmlrpcFunc
 * enum for functions.
 * This specifies the status
 */
enum exmlrpcTWFunc {
  /** inactive */
  XMLRPC_Inactive, 
  /** Ping argument ignored */
  XMLRPC_Ping, 
  /** Create a task message window */
  XMLRPC_CreateWindow,
  /** Display task message */
  XMLRPC_DisplayMessage,
  /** Request User response */
  XMLRPC_UserResponse,
  /** Set task status */
  XMLRPC_SetStatus,
}; /* end enum exmlrpcTWFunc */
/** typedef for enum for ObitIO object status. */
typedef enum exmlrpcTWFunc XMLRPCFunc;

/*----------------- Globals ---------------------------*/
#ifdef OBITMESSMAIN   /* Definitions for main program */
/** mutex to lock request data */
pthread_mutex_t request_lock;

/** flag for request recieved,  1=ready, 0=wait */
glong receive_flag;

/** flag for request return ready, 1=ready, 0=wait */
glong return_flag;

/** User request timeout in sec, <0 -> forever */
ofloat  ERtime_out;  

/** Request function code */
XMLRPCFunc xmlrpcTWFunc;

/** Communication between threads pointer */
void *xmlrpcData;

/** TaskID for request */
glong taskID;

/** Status from handling request, 0=OK else fail */
glong xmlStatus;

/** Abort request */
gboolean doAbort;

#else                 /* Definitions for routines */
extern pthread_mutex_t request_lock;
extern glong    receive_flag;
extern glong    return_flag;
extern ofloat   ERtime_out;
extern XMLRPCFunc xmlrpcTWFunc;
extern void     *xmlrpcData;
extern glong    taskID;
extern glong xmlStatus;
extern gboolean doAbort;
#endif  /* OBITMESSMAIN  */

/*---------------RPC functions--------------------------- */
/**
 * ObitMess RPC interface:
 * \li ping  
 * Ping to see if server is alive and taskID window active
 * argument: XML int - taskID
 * \li CreateWindow  
 * Create new task message window
 * argument: XML string = task name
 * returns taskID in reply
 * \li DisplayMessage 
 * Display a message in display server for task taskID
 * Argument:  XML structure {"taskID":taskID_int,"message":message_string}
 * \li SetStatus  
 * Sets task status string
 * Argument: structure {"taskID":taskID_int,"status":status_string}
 * \li UserResponse  
 * Solicit user response string
 * argument: XML int - taskID from Create Window
 * string returned as "Result"
 *
 * Returns XML structure:
 * {
 * "Status":{'reason':reason_for failure_string("OK"=OK),"code":return_code(0=OK)},
 * "Abort":user_request_abort_boo (0=false),
 * "Result":any_reply_string,
 * "taskID":task_ID_for_task_window
 * }
 */
/*---------------Public functions---------------------------*/
/** Start server function - to run in separate thread */
void* start_myTWServer (void* port);
/** Watch for requests from XMLRPC thread */
void* XMLRPCTWWatcher   (void* arg);
void  XMLRPCTWRelease   (void* data);
void  XMLRPCTWSetReturn (void* valueP);

#endif /* end XMLRPCTASKMESSSERVER_H */


