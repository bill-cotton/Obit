/* How to pass OK/failed status */
/* $Id:  $ */
/* XMLRPC server for ObitMess */
/* Much of this material directly adapted from xmlrpc-c-1.2/examples */
/*-----------------------------------------------------------------------
*  Copyright (C) 2009-2012
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
*  Correspondence concerning ObitView should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#include <unistd.h> 
#include <errno.h>    /* global C error include */
#include <stdio.h>    /* i/o include */
#include <stdlib.h>   /* i/o include */
#include <signal.h>
#include <xmlrpc_client.h>
#include <xmlrpc_server.h>
/*#include <xmlrpc_server_abyss.h>*/
/*#include <xmlrpc-c/base.h>
  #include <xmlrpc-c/server.h>
  #include <xmlrpc-c/server_abyss.h>*/
#include "XMLRPCTaskMessServer.h"
#include "ObitRPC.h"
#include "obitmess.h"
#include "messagebox.h"
#include "taskmessage.h"
#include "xml.h"
#include "ObitMem.h"
#include "ObitInfoList.h"
#include <Xm/Xm.h> 

/*---------------Private function prototypes----------------*/
static void server_loop(gint port);
static xmlrpc_value*
ping(xmlrpc_env *  const envP, 
     xmlrpc_value * const paramArrayP,
     void *const userData);

static xmlrpc_value *
createWindow(xmlrpc_env *   const envP, 
	     xmlrpc_value * const paramArrayP,
	     void *         const userData);

static xmlrpc_value *
displayMessage(xmlrpc_env *   const envP, 
	       xmlrpc_value * const paramArrayP,
	       void *         const userData);

static xmlrpc_value *
userResponse(xmlrpc_env *   const envP, 
	     xmlrpc_value * const paramArrayP,
	     void *         const userData);

static xmlrpc_value *
setStatus(xmlrpc_env *   const envP, 
	     xmlrpc_value * const paramArrayP,
	     void *         const userData);

/*----------------- Globals ---------------------------*/
/* Are we busy doing something? */
gboolean WeAreBusy;

/*---------------Public functions ----------------*/
/**
 * start server thread
 * \param port  TCP port to watch
 */
void* start_myTWServer (void* port)
{

  /* initialization */
  WeAreBusy    = FALSE;
  doAbort      = FALSE;
  xmlStatus    = 0;
  taskID       = -1;
  receive_flag = 0;
  return_flag  = 0;
  xmlrpcTWFunc = XMLRPC_Inactive;
  if (pthread_mutex_init (&request_lock, NULL)) {
    fprintf (stderr, "Failed to initialize Request lock \n");
    return NULL;
  }

  /* Start looping */
  server_loop (*(gint*)port);

  return NULL;
} /* end start_server */

/**
 * start server and begin looping
 * \param port  TCP port to watch
 */
static void server_loop(gint port) 
{
  ObitRPC *server=NULL;
  ObitErr *err=NULL;

  /* Create Server */
  err = newObitErr();
  server = ObitRPCCreateServer("ObitMess", err);

  ObitRPCAddMethod (server, "ping",           &ping,           NULL, err);
  ObitRPCAddMethod (server, "CreateWindow",   &createWindow,   NULL, err);
  ObitRPCAddMethod (server, "DisplayMessage", &displayMessage, NULL, err);
  ObitRPCAddMethod (server, "UserResponse",   &userResponse,   NULL, err);
  ObitRPCAddMethod (server, "SetStatus",      &setStatus,      NULL, err);

  /* Loop forever 'neath the streets of Boston */
  /*ObitRPCServerLoop(server, port, "/tmp/xmlrpc_log");*/
  ObitRPCServerLoop(server, port, "/dev/null");

} /* end server_loop */

/**
 * Function to run in X thread to act on requests.
 * Communication with the XMLRPC thread is through globals.
 * (time out loop/timer in ObitMess )
 * \li request_lock - mutex to protect global values
 * \li xmlrpcTWFunc - code of function called
 * \li xmlrpcData   - data pointer
 * \li receive_flag - boolean flag indicating a request
 * \li return_flag  - boolean flag indicating a return value is ready
 */
void* XMLRPCTWWatcher (XtPointer clientData) 
{
  static gchar *IAm = "ObitMess";
  Widget parent = (Widget)clientData;

  /*while(1) {*/

  if (receive_flag) {
    pthread_mutex_lock(&request_lock); /* lock mutex */
    receive_flag=0;  /* Clear receive flag */
    
    xmlStatus = 0;  /* unless proved otherwise */
    /* Branch by function */
    switch (xmlrpcTWFunc) {
    case XMLRPC_Inactive:        /* Called by mistake */
      break;
    case XMLRPC_Ping:            /* Ping arg taskID */
      xmlrpcData = IAm;
      return_flag = 1;   /* Set return flag */
      break;
    case XMLRPC_CreateWindow:    /* Create Window, arg taskname, return taskID */
      taskID = TaskMessageCreate(parent, xmlrpcData, &xmlStatus);
      return_flag = 1;   /* Set return flag */
      break;
    case XMLRPC_DisplayMessage:  /* Display message args taskID, message */
      TaskMessageDisplay(taskID, xmlrpcData, &xmlStatus);
      return_flag = 1;   /* Set return flag */
      break;
    case XMLRPC_UserResponse:    /* request user input arg taskID */
      TaskMessageUserResponseSetup(taskID, &xmlStatus);
      /* Reply sent from call back */
      break;
    case XMLRPC_SetStatus:  /* Set task status args taskID, status */
      TaskMessageStatus(taskID, xmlrpcData, &xmlStatus);
      return_flag = 1;   /* Set return flag */
      break;

    } /* end switch */
    pthread_mutex_unlock(&request_lock); /* unlock mutex */
  } 
  /* Let's not burn too many cycles in this event loop */
  usleep(250000);  /* 250 msec */
  /* } end of loop forever */
  return NULL;
} /* end XMLRPCTWWatcher */


/**
 * Set return ready flag
 * \param data  if nonNULL use as xmlrpcData
 */
void XMLRPCTWRelease (void *data)
{  
  /* Something going on? */
  if (xmlrpcTWFunc == XMLRPC_Inactive) return;
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  return_flag = 1;                     /* Send result if requested via XMLRPC */
  if (data) xmlrpcData = data;         /* Set return value */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */
  xmlrpcTWFunc = XMLRPC_Inactive;

} /* end XMLRPCTWRelease */

/**
 * Set return ready flag and value pointer
 * \param valueP  pointer to the item to be returned
 */
void XMLRPCTWSetReturn (void* valueP)
{  
  /* Something going on? */
  if (xmlrpcTWFunc == XMLRPC_Inactive) return;
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  if (valueP) xmlrpcData = valueP;     /* Set return value */
  return_flag = 1;                     /* Set return flag */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */
  xmlrpcTWFunc = XMLRPC_Inactive;

} /* end XMLRPCTWSetReturn */

/*------------------- Private functions -------------------*/

/**
 * Handle ping function call 
 * Argument is taskID
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
ping(xmlrpc_env *   const envP, 
     xmlrpc_value * const paramArrayP,
     void *         const userData) 
{  
  xmlrpc_int32 xmlTaskID;
  gboolean ldoAbort;
  glong status=0;

  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy");
  WeAreBusy = TRUE; /* busy now */
  
  /* Parse our argument array. */
  xmlrpc_decompose_value(envP, paramArrayP, "(i)", &xmlTaskID);
  if (envP->fault_occurred)
    return NULL;

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock); /* lock mutex */
  receive_flag = 1;          /* now have a request */
  xmlrpcTWFunc = XMLRPC_Ping;  /* Function called */
  taskID = (glong)xmlTaskID;              /* TaskID */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result. */
      ldoAbort =  TaskMessagedoAbort (taskID, &status);
      WeAreBusy = FALSE; /* no longer busy */
      xmlrpcTWFunc = XMLRPC_Inactive;
      /* It work? */
      if (status==0) { /* OK */
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)0,"reason","OK",
				  "Result", xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)ldoAbort);
      } else { /* FAILED */ 
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)status,
				     "reason","Not Alive",
				  "Result", xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end ping */

/**
 * Handle CreateWindow function call 
 * Argument is task name
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
createWindow(xmlrpc_env *   const envP, 
	     xmlrpc_value * const paramArrayP,
	     void *         const userData)
{  
  gchar *taskname;
  gchar savestr[200]; /* Storage for string arg*/
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
 				"Result", "Busy",
				"taskID", (xmlrpc_int32)taskID,
				"Abort",  (xmlrpc_int32)doAbort);
   WeAreBusy = TRUE; /* busy now */
 

  /* Parse our argument array. */
  xmlrpc_decompose_value(envP, paramArrayP, "(s)", &taskname);
  if (envP->fault_occurred)
    return NULL;
  /* Save argument */
  strncpy(savestr, taskname, 199);

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  receive_flag = 1;                    /* now have a request */
  xmlrpcTWFunc = XMLRPC_CreateWindow;  /* Function called */
  xmlrpcData   = savestr;              /* pass task name */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
      g_free (taskname);   /* Cleanup */
     
      /* Return our result. */
      WeAreBusy    = FALSE; /* no longer busy */
      xmlrpcTWFunc = XMLRPC_Inactive;
      /* It work? */
      if (xmlStatus!=0)  /* FAILED */
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)xmlStatus,
				     "reason","Failed",
				  "Result", (char*)xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      else { /* OK */
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)0,"reason","OK",
				  "Result", (char*)xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end createWindow */

/**
 * Handle DisplayMessage function call 
 * Argument is {"taskID":taskID, "message":message} message should include \n
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
displayMessage(xmlrpc_env *   const envP, 
	       xmlrpc_value * const paramArrayP,
	       void *         const userData)
{  
  xmlrpc_int32 xmlTaskID;
  xmlrpc_value *strt;
  xmlrpc_type xmlType;
  char *xmlMessage;
  gboolean ldoAbort;
  glong status=0;
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy",
				"taskID", (xmlrpc_int32)taskID,
				"Abort",  (xmlrpc_int32)doAbort);
  WeAreBusy = TRUE; /* busy now */
 
  /* Parse our argument structure. */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &strt);

  /* Be sure paramArrayP an array */
  xmlType = xmlrpc_value_type (strt);
  xmlMessage = NULL;
  if (xmlType==XMLRPC_TYPE_STRUCT) {
    xmlrpc_decompose_value(envP, strt, "{s:i,s:s,*}",
			   "taskID", &xmlTaskID,
			   "message", &xmlMessage);
  } 

  if (envP->fault_occurred) {
    ObitErrLog (err);
    WeAreBusy = FALSE; /* no longer busy */
    xmlrpcTWFunc = XMLRPC_Inactive;
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", "Bad",
			      "taskID", (xmlrpc_int32)0,
			       "Abort", (xmlrpc_int32)0);

  }

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);      /* lock mutex */
  receive_flag = 1;                       /* now have a request */
  taskID = (glong)xmlTaskID;              /* TaskID */
  xmlrpcData   = xmlMessage;              /* pass message */
  xmlrpcTWFunc = XMLRPC_DisplayMessage;   /* Function called */
  pthread_mutex_unlock(&request_lock);    /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result.- string set in XMLRPCWatcher */
      /* It work? */
      if (xmlStatus==0) { /* OK */
	ldoAbort =  TaskMessagedoAbort (taskID, &status);
	WeAreBusy = FALSE; /* no longer busy */
	xmlrpcTWFunc = XMLRPC_Inactive;
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", "Displayed",
				    "taskID", (xmlrpc_int32)taskID,
				    "Abort",  (xmlrpc_int32)ldoAbort);
      } else { /* FAILED */ 
	WeAreBusy = FALSE; /* no longer busy */
	xmlrpcTWFunc = XMLRPC_Inactive;
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)xmlStatus,
				     "reason","Failed",
				  "Result", "Failed",
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end displayMessage */


/**
 * Handle UserResponse function call 
 * Argument is taskID
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
userResponse(xmlrpc_env *   const envP, 
	       xmlrpc_value * const paramArrayP,
	       void *         const userData)
{  
  xmlrpc_int32 xmlTaskID;
  gboolean ldoAbort;
  glong status=0;
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy",
				"taskID", (xmlrpc_int32)taskID,
				"Abort",  (xmlrpc_int32)doAbort);
  WeAreBusy = TRUE; /* busy now */
 
  /*fprintf (stderr, "DEBUG in loadImage\n");*/

  /* parse information */
  xmlrpc_decompose_value(envP, paramArrayP, "(i)", &xmlTaskID);
  if (envP->fault_occurred) {
    ObitErrLog (err);
    WeAreBusy = FALSE; /* no longer busy */
    xmlrpcTWFunc = XMLRPC_Inactive;
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", "Bad",
			      "taskID", (xmlrpc_int32)0,
			      "Abort",  (xmlrpc_int32)0);

  }

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);  /* lock mutex */
  receive_flag = 1;                   /* now have a request */
  taskID = (glong)xmlTaskID;              /* TaskID */
  xmlrpcTWFunc = XMLRPC_UserResponse;  /* Function called */
  pthread_mutex_unlock(&request_lock);/* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result.- string set in XMLRPCTWWatcher */
      WeAreBusy = FALSE; /* no longer busy */
      /* It work? */
      if (xmlStatus==0) { /* OK */
	ldoAbort =  TaskMessagedoAbort (taskID, &status);
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)0,"reason","OK",
				  "Result", (char*)xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)ldoAbort);
      } else {/* FAILED */ 
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)xmlStatus,
				     "reason","Failed",
				  "Result", "Failed",
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end userResponse */

/**
 * Handle setStatus function call 
 * Argument is {"taskID":taskID, "status":status_string}
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
setStatus(xmlrpc_env *   const envP, 
	  xmlrpc_value * const paramArrayP,
	  void *         const userData)
{  
  xmlrpc_value *strt;
  xmlrpc_type xmlType;
  xmlrpc_int32 xmlTaskID;
  char *xmlMessage;
  gboolean ldoAbort;
  glong status=0;
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy",
				"taskID", (xmlrpc_int32)taskID,
				"Abort",  (xmlrpc_int32)doAbort);
  WeAreBusy = TRUE; /* busy now */
 
  /* Parse our argument structure. */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &strt);
  /* Be sure paramArrayP an array */
  xmlType = xmlrpc_value_type (strt);
  xmlMessage = NULL;
  if (xmlType==XMLRPC_TYPE_STRUCT) {
    xmlrpc_decompose_value(envP, strt, "{s:i,s:s,*}",
			   "taskID", &xmlTaskID,
			   "status", &xmlMessage);
  } 

  if (envP->fault_occurred) {
    ObitErrLog (err);
    WeAreBusy = FALSE; /* no longer busy */
    xmlrpcTWFunc = XMLRPC_Inactive;
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", "Invalid argument",
			      "taskID", (xmlrpc_int32)taskID,
			       "Abort", (xmlrpc_int32)doAbort);

  }

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);      /* lock mutex */
  receive_flag = 1;                       /* now have a request */
  taskID = (glong)xmlTaskID;              /* TaskID */
  xmlrpcData   = xmlMessage;              /* pass message */
  xmlrpcTWFunc = XMLRPC_SetStatus;        /* Function called */
  pthread_mutex_unlock(&request_lock);    /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result.- string set in XMLRPCWatcher */
      WeAreBusy = FALSE; /* no longer busy */
      /* It work? */
      if (xmlStatus==0) { /* OK */
	ldoAbort =  TaskMessagedoAbort (taskID, &status);
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", (char*)xmlrpcData,
				    "taskID", (xmlrpc_int32)taskID,
				    "Abort",  (xmlrpc_int32)ldoAbort);
      } else { /* FAILED */ 
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s,s:i,s:i}", 
				  "Status","code", (xmlrpc_int32)xmlStatus,
				     "reason","Failed",
				  "Result", (char*)xmlrpcData,
				  "taskID", (xmlrpc_int32)taskID,
				  "Abort",  (xmlrpc_int32)doAbort);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end setStatus */



