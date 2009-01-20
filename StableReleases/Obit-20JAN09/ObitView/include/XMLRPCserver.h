/* $Id$  */
/* XMLRPC server for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005-2008
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

#ifndef XMLRPCSERVER_H
#define XMLRPCSERVER_H
/*------------------ Structures -----------------------------*/


/*-------------- enumerations -------------------------------------*/
/**
 * \enum xmlrpcFunc
 * enum for functions.
 * This specifies the status of the connection a disk resident data.
 */
enum xmlrpcFunc {
  /** inactive */
  XMLRPC_Inactive, 
  /** Ping argument ignored */
  XMLRPC_Ping, 
  /** Load Fits file - argument = filename */
  XMLRPC_LoadFits,
  /** Load Image file - arguments = file info */
  XMLRPC_LoadImage,
  /** Edit Clean window - argument = window */
  XMLRPC_EditWindow,
  /** Copy File - argument = binary blob and description */
  XMLRPC_FileCopy
}; /* end enum xmlrpcFunc */
/** typedef for enum for ObitIO object status. */
typedef enum xmlrpcFunc XMLRPCFunc;

/*----------------- Globals ---------------------------*/
#ifdef OBITVIEWMAIN   /* Definitions for main program */
/** mutex to lock request data */
pthread_mutex_t request_lock;

/** flag for request recieved,  1=ready, 0=wait */
olong receive_flag;

/** flag for request return ready, 1=ready, 0=wait */
olong return_flag;

/** Edit request timeout in sec, <0 -> forever */
ofloat  ERtime_out;  

/** Request function code */
XMLRPCFunc xmlrpcFunc;

/** Communication between threads pointer */
void *xmlrpcData;

#else                 /* Definitions for routines */
extern pthread_mutex_t request_lock;
extern olong receive_flag;
extern olong return_flag;
extern ofloat  ERtime_out;
extern XMLRPCFunc xmlrpcFunc;
extern void *xmlrpcData;
#endif  /* OBITVIEWMAIN  */

/*---------------Public functions---------------------------*/
/** Start server function - to run in separate thread */
void* start_myServer (void* port);
/** Watch for requests from XMLRPC thread */
void* XMLRPCWatcher (void* arg);
void XMLRPCRelease (void *data);
void XMLRPCSetReturn (void* valueP);

#endif /* end XMLRPCSERVER_H */


