/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITRPC_H 
#define OBITRPC_H 

#include <xmlrpc_client.h>
#include <xmlrpc_server.h>
#include <xmlrpc_server_abyss.h>
#include <xmlrpc-c/base.h>
#include <xmlrpc-c/client.h>
#include <xmlrpc-c/server.h>
#include <xmlrpc-c/server_abyss.h>
#include "Obit.h"
#include "ObitErr.h"
#include "ObitXML.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRPC.h
 * ObitRPC Remote Procedure Call class
 *
 * This class is derived from the #Obit class.
 * Related functions are in the 
 * \link ObitRPCUtil.h ObitRPCUtil 
 * \endlink module.
 *
 * This class handles Remote Procedure Calls.
 * The implementation is based on xmlrpc which only allows a single 
 * client or server handler at a time.
 * 
 * \section ObitRPCaccess Creators and Destructors
 * An ObitRPC will usually be created using ObitRPCCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitRPC should always be made using the
 * #ObitRPCRef function which updates the reference count in the object.
 * Then whenever freeing an ObitRPC or changing a pointer, the function
 * #ObitRPCUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitRPCType
 * enum for object type.
 * This specifies which known type of RPC object
 */
enum obitRPCType {
  /** Client */
  OBIT_RPC_Client = 0, 
  /** Server */
  OBIT_RPC_Server 
}; /* end enum obitRPCType */
/** typedef for enum for ObitRPCType object status. */
typedef enum obitRPCType ObitRPCType;

/**
 * \enum obitRPCRequestType
 * enum for request sent from server back to client
 * This list needs to be compatible with #ObitDisplayRequest
 * defined in ObitDisplay.h
 */
enum obitRPCRequestType {
  /** Continue the program */
  OBIT_RPC_Request_Continue = 0, 
  /** Abort the program */
  OBIT_RPC_Request_Abort, 
  /** Quit - graceful shutdown */
  OBIT_RPC_Request_Quit,
  /** Stop displaying intermediate results */
  OBIT_RPC_Request_NoTV,
  /** Send another field in the current mosaic */
  OBIT_RPC_Request_Field,
  /** Edit Clean window in current field */
  OBIT_RPC_Request_EditWin
}; /* end enum obitRPCRequestType */
/** typedef for enum for ObitRPCRequestType  */
typedef enum obitRPCRequestType ObitRPCRequestType;

typedef xmlrpc_response_handler ObitRPC_response_handler;

/*--------------Class definitions-------------------------------------*/
/** ObitRPC Class structure. */
typedef struct {
#include "ObitRPCDef.h"   /* this class definition */
} ObitRPC;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitRPC
 * returns a ObitRPC*.
 * in = object to unreference
 */
#define ObitRPCUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitRPC.
 * returns a ObitRPC*.
 * in = object to reference
 */
#define ObitRPCRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitRPCIsA(in) ObitIsA (in, ObitRPCGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitRPCClassInit (void);

/** Public: Default Constructor. */
ObitRPC* newObitRPC (gchar* name);

/** Public: Create/initialize client ObitRPC structures */
ObitRPC* ObitRPCCreateClient (gchar* name, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitRPC* (*ObitRPCCreateClientFP) (gchar* name, ObitErr *err);

/** Public: Create/initialize server ObitRPC structures */
ObitRPC* ObitRPCCreateServer (gchar* name, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitRPC* (*ObitRPCCreateServerFP) (gchar* name, ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitRPCGetClass (void);

/** Public: Send synchronous RPC request */
ObitXML* ObitRPCCall (ObitRPC* client, gchar *serverURL, ObitXML* arg, 
		      ObitInfoList **status, ObitInfoList **request,
		      ObitErr *err);

/** Public: Send asynchronous RPC request */
void ObitRPCCallSnd (ObitRPC* client, gchar *serverURL, ObitXML* arg, 
		     ObitRPC_response_handler callback, gpointer user_data,
		     ObitErr *err);

/** Public: Process asynchronous RPC response */
ObitXML* ObitRPCCallRcv (ObitRPC* client, ObitXMLValue *result, 
			 ObitInfoList **status, ObitInfoList **request,
			 ObitErr *err);

/** Add method callback to server */
void  ObitRPCAddMethod (ObitRPC* server, gchar *method_name, 
			xmlrpc_method method, gpointer user_data,
			ObitErr *err);

/** Start Server loop */
void  ObitRPCServerLoop (ObitRPC* server, olong port, gchar *log_file);

/** Run client async event loop */
void  ObitRPCClientAsyncLoop (olong timeout);

/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitRPCClassDef.h"
} ObitRPCClassInfo; 

#endif /* OBITFRPC_H */ 
