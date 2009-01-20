/* $Id$    */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitRPC.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRPC.c
 * ObitRPC class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitRPC";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitRPCClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitRPCClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitRPCInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitRPCClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitRPCClassInfoDefFn (gpointer inClass);

static void 
die_if_fault_occurred (xmlrpc_env * const envP) {
    if (envP->fault_occurred) {
        fprintf(stderr, "Something failed. %s (XML-RPC fault code %d)\n",
                envP->fault_string, envP->fault_code);
        exit(1);
    }
}
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitRPC* newObitRPC (gchar* name)
{
  ObitRPC* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitRPCClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitRPC));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitRPCInit((gpointer)out);

 return out;
} /* end newObitRPC */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitRPCGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitRPCClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitRPCGetClass */

/**
 * Creates an ObitRPC Client
 * \param name  An optional name for the object.
 * \param err   Obit Error message
 * \return the new object.
 */
ObitRPC* ObitRPCCreateClient (gchar* name, ObitErr *err)
{
  ObitRPC* out;
  struct xmlrpc_clientparms clientparmsP;
  gchar *routine = "ObitRPCCreateClient";

  /* Only initialize actual client once - there is only one
   probably need locking for multi threaded operation */
  if (myClassInfo.numberClient<1) {
    /* Start up our XML-RPC client library.
    xmlrpc_client_init(XMLRPC_CLIENT_NO_FLAGS, "Obit XML client", 
		       VERSION);  done in ClassInit */
    
  }
  
  /* Create basic structure */
  out = newObitRPC (name);
  
  /* Type dependent initialization */
  out->type    = OBIT_RPC_Client;

  /* Create private client */
  clientparmsP.transport = "curl";  /* Force using cUrl */
  xmlrpc_client_create(&out->envP, XMLRPC_CLIENT_NO_FLAGS, "Obit", "1.0", 
		       &clientparmsP, sizeof(clientparmsP.transport),   /* Force transport */ 
		       &out->clientP);
  /* Make sure it worked */
  Obit_retval_if_fail((!out->envP.fault_occurred ),
		      err, out, "%s: XML-RPC Fault: %s (%d)",
		      routine, out->envP.fault_string, 
		      out->envP.fault_code);

  
  myClassInfo.numberClient++;  /* Keep track */
  
  return out;
} /* end ObitRPCCreateClient */

/**
 * Creates an ObitRPC Server
 * \param name  An optional name for the object.
 * \param err   Obit Error message
 * \return the new object.
 */
ObitRPC* ObitRPCCreateServer (gchar* name, ObitErr *err)
{
  ObitRPC* out;
  /*gchar *routine = "ObitRPCCreateServer";*/


  /* Create basic structure */
  out = newObitRPC (name);

  /* Type dependent initialization */
  out->type = OBIT_RPC_Server;
  
  /* define registry */
  out->registryP = xmlrpc_registry_new(&out->envP);
    
  /* Set server parameters */
  out->serverparm.config_file_name = NULL;
  out->serverparm.registryP        = out->registryP;
  out->serverparm.port_number      = 0;
  out->serverparm.log_file_name    = NULL;
 
  myClassInfo.numberServer++;

  return out;
} /* end ObitRPCCreateServer */

/**
 * Make synchronous remote procedure call
 * Uses private client
 * Return value from RPC Call expects an xml struct with up to 3 parts:
 * \li "Status", 
 *       "code" an integer code (0=OK) 
 *       "reason" a status string
 * \li "Result" any result, depends on call
 * \li "Request", a request to the client
 *       "code" obitRPCRequestType
 *        parameters as needed for request
 * \param client     Client ObitRPC
 * \param serverURL  URL of service, e.g. "http://localhost:8765/RPC2"
 * \param arg        Argument of call (includes method name )
 * \param status     [out] if non NULL, Status in form of Info list, 
 *                   entries "code", "reason"
 *                   Should be Unrefed when done
 * \param request    [out] if non NULL, Any in form of Info list, 
 *                   entries "code", and case dependent parameters.
 *                   Returns NULL if no request
 *                   Should be Unrefed when done
 * \param err        Obit Error message
 * \return XML object returned, NULL on communications failure,
 *            even if this is defined the function may have failed.
 */
ObitXML* ObitRPCCall (ObitRPC* client, gchar *serverURL, ObitXML* arg, 
		      ObitInfoList **status, ObitInfoList **request,
		      ObitErr *err)
{
  ObitXML* out=NULL, *tempXML=NULL;
  xmlrpc_value *returnP=NULL, *tempP;
  xmlrpc_type xmlType;
  gchar *routine = "ObitRPCCall";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  Obit_retval_if_fail((client->type==OBIT_RPC_Client), 
		      err, out, "%s: RPC NOT a client", routine);
  XMLRPC_FAIL_IF_FAULT(&client->envP);

  /* Make the remote procedure call */
  xmlrpc_client_call2f(&client->envP, client->clientP, serverURL, 
		       arg->func, &returnP, "(V)", arg->parmP);

  /* Make sure it worked */
  Obit_retval_if_fail((!client->envP.fault_occurred && (returnP!=NULL)),
		      err, out, "%s: XML-RPC Fault: %s (%d)",
		      routine, client->envP.fault_string, 
		      client->envP.fault_code);

  /* Get returned result - better be a struc with key "Result" */
  xmlType = xmlrpc_value_type (returnP);
  Obit_retval_if_fail((xmlType==XMLRPC_TYPE_STRUCT),
		      err, out, "%s: return NOT a struct", routine);
  Obit_retval_if_fail((xmlrpc_struct_has_key(&client->envP, returnP, "Result")),
		      err, out, "%s: return has no Result member", routine);

  xmlrpc_struct_find_value (&client->envP, returnP, "Result", &tempP);
  XMLRPC_FAIL_IF_FAULT(&client->envP);
 
  /* Construct return object */
  out = ObitXMLReturn (arg->func, tempP, err);
  xmlrpc_DECREF(tempP);

  /* Get status if requested */
  if (status) {
    if (xmlrpc_struct_has_key(&client->envP, returnP, "Status")) {
      xmlrpc_struct_find_value (&client->envP, returnP, "Status", &tempP);
      XMLRPC_FAIL_IF_FAULT(&client->envP);
      tempXML = ObitXMLReturn (arg->func, tempP, err);
      xmlrpc_DECREF(tempP);
      tempXML->type = OBIT_XML_InfoList;
      *status = ObitXMLXML2InfoList (tempXML, err);
      tempXML = ObitXMLUnref(tempXML);
    }
    if (err->error) Obit_traceback_val (err, routine, "Get status", out);  
  } /* end get status */

  /* Get request if requested - should be an InfoList on the other end */
  if (request) {
    if (xmlrpc_struct_has_key(&client->envP, returnP, "Request")) {
      xmlrpc_struct_find_value (&client->envP, returnP, "Request", &tempP);
      XMLRPC_FAIL_IF_FAULT(&client->envP);
      tempXML = ObitXMLReturn (arg->func, tempP, err);
      xmlrpc_DECREF(tempP);
      tempXML->type = OBIT_XML_InfoList;
      *request = ObitXMLXMLInfo2List (tempXML, err);
      tempXML = ObitXMLUnref(tempXML);
    }
  } /* end get request */


  /* Make sure everything is cool */
  cleanup:
  if (client->envP.fault_occurred) {
    Obit_log_error(err, OBIT_StrongError, "XML-RPC Fault: %s (%d)",
		   client->envP.fault_string, client->envP.fault_code);
  } else {
    xmlrpc_DECREF(returnP);
  }

  return out;
} /* end ObitRPCCall */

/**
 * Make asynchronous remote procedure call request
 * Uses global client.
 * Request is aynchronously sent; the reply will be to callback
 * \param client     Client ObitRPC
 * \param serverURL  URL of service, e.g. "http://localhost:8765/RPC2"
 * \param arg        Argument of call (includes method name )
 * \param callback   Callback function
 * typedef void (*xmlrpc_response_handler) (const char *server_url,
 *                                          const char *method_name,
 *                                          xmlrpc_value *param_array,
 *                                          void *user_data,
 *                                          xmlrpc_env *fault,
 *                                          xmlrpc_value *result);
 * \param user_data  Pointer to be passed to Callback with result
 * \param err        Obit Error message
 * \return XML object returned, NULL on communications failure,
 *            even if this is defined the function may have failed.
 */
void ObitRPCCallSnd (ObitRPC* client, gchar *serverURL, ObitXML* arg, 
		     ObitRPC_response_handler callback, gpointer user_data,
		     ObitErr *err)
{
  gchar *routine = "ObitRPCCallSnd";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  Obit_return_if_fail((client->type==OBIT_RPC_Client), 
		      err, "%s: RPC NOT a client", routine);
  XMLRPC_FAIL_IF_FAULT(&client->envP);

  /* Make asynchronous call */
  xmlrpc_client_call_asynch (serverURL, arg->func, callback, user_data,
    "(V)", arg->parmP);

cleanup:
  if (client->envP.fault_occurred) {
    Obit_log_error(err, OBIT_StrongError, "XML-RPC Fault: %s (%d)",
		   client->envP.fault_string, client->envP.fault_code);
  } 
} /* end ObitRPCCallSnd */

/**
 * Process response from asynchronous call
 * Also kills xmlrpc_client_event_loop
 * Return value from RPC Call expects an xml struct with up to 3 parts:
 * \li "Status", 
 *       "code" an integer code (0=OK) 
 *       "reason" a status string
 * \li "Result" any result, depends on call
 * \li "Request", a request to the client
 *       "code" obitRPCRequestType
 *        parameters as needed for request
 * \param client     Client ObitRPC
 * \param result     Response from callback as gpointer
 * \param status     [out] if non NULL, Status in form of Info list, 
 *                   entries "code", "reason"
 *                   Should be Unrefed when done
 * \param request    [out] if non NULL, Any in form of Info list, 
 *                   entries "code", and case dependent parameters.
 *                   Returns NULL if no request
 *                   Should be Unrefed when done
 * \param err        Obit Error message
 * \return XML object returned, NULL on communications failure,
 *            even if this is defined the function may have failed.
 */
ObitXML* ObitRPCCallRcv (ObitRPC* client, ObitXMLValue *result, 
			 ObitInfoList **status, ObitInfoList **request,
			 ObitErr *err)
{
  ObitXML* out=NULL, *tempXML=NULL;
  xmlrpc_value *returnP=NULL, *tempP;
  xmlrpc_type xmlType;
  gchar *routine = "ObitRPCCallRcv";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;

  returnP = (xmlrpc_value*)result;
  
  /* Make sure it worked */
  Obit_retval_if_fail((!client->envP.fault_occurred && (returnP!=NULL)),
		      err, out, "%s: XML-RPC Fault: %s (%d)",
		      routine, client->envP.fault_string, 
		      client->envP.fault_code);

  /* Get returned result - better be a struc with key "Result" */
  xmlType = xmlrpc_value_type (returnP);
  Obit_retval_if_fail((xmlType==XMLRPC_TYPE_STRUCT),
		      err, out, "%s: return NOT a struct", routine);
  Obit_retval_if_fail((xmlrpc_struct_has_key(&client->envP, returnP, "Result")),
		      err, out, "%s: return has no Result member", routine);

  xmlrpc_struct_find_value (&client->envP, returnP, "Result", &tempP);
  XMLRPC_FAIL_IF_FAULT(&client->envP);
 
  /* Construct return object */
  out = ObitXMLReturn ("XMLReturn", tempP, err);
  xmlrpc_DECREF(tempP);

  /* Get status if requested */
  if (status) {
    if (xmlrpc_struct_has_key(&client->envP, returnP, "Status")) {
      xmlrpc_struct_find_value (&client->envP, returnP, "Status", &tempP);
      XMLRPC_FAIL_IF_FAULT(&client->envP);
      tempXML = ObitXMLReturn ("XMLReturn", tempP, err);
      xmlrpc_DECREF(tempP);
      tempXML->type = OBIT_XML_InfoList;
      *status = ObitXMLXML2InfoList (tempXML, err);
      tempXML = ObitXMLUnref(tempXML);
    }
    if (err->error) Obit_traceback_val (err, routine, "Get status", out);  
  } /* end get status */

  /* Get request if requested - should be an InfoList on the other end */
  if (request) {
    if (xmlrpc_struct_has_key(&client->envP, returnP, "Request")) {
      xmlrpc_struct_find_value (&client->envP, returnP, "Request", &tempP);
      XMLRPC_FAIL_IF_FAULT(&client->envP);
      tempXML = ObitXMLReturn ("XMLReturn", tempP, err);
      xmlrpc_DECREF(tempP);
      tempXML->type = OBIT_XML_InfoList;
      *request = ObitXMLXMLInfo2List (tempXML, err);
      tempXML = ObitXMLUnref(tempXML);
    }
  } /* end get request */


  /* Make sure everything is cool */
  cleanup:
  if (client->envP.fault_occurred) {
    Obit_log_error(err, OBIT_StrongError, "XML-RPC Fault: %s (%d)",
		   client->envP.fault_string, client->envP.fault_code);
  } else {
    xmlrpc_DECREF(returnP);
  }

  return out;
} /* end ObitRPCCallRcv */

/**
 * Adds method callback to server
 * \param server       Server ObitRPC
 * \param method_name  name of the method
 * \param method       function pointer to xmlrpc_method callback
 * \param user_data    Additional data to be passed to method, may be NULL
 * \param err          Obit Error message
 */
void  ObitRPCAddMethod (ObitRPC* server, gchar *method_name, 
			xmlrpc_method method, gpointer user_data,
			ObitErr *err)
{
  gchar *routine = "ObitRPCAddMethod";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  Obit_return_if_fail((server->type==OBIT_RPC_Server), 
		      err, "%s: RPC NOT a server",  routine);

  /* Add to registry */
  xmlrpc_registry_add_method(&server->envP, server->registryP, NULL, 
			     method_name, method, user_data);

  /* Make sure it worked */
  Obit_return_if_fail((!server->envP.fault_occurred),
		      err, "%s: XML-RPC Fault: %s (%d)",
		      routine, server->envP.fault_string, 
		      server->envP.fault_code);

} /* end ObitRPCAddMethod */

/**
 * Starts the server loop (never returns)
 * \param server  Server ObitRPC
 * \param port    Port number to listen to
 * \param         logging file, NULL=>/tmp/xmlrpc_log
 */
void  ObitRPCServerLoop (ObitRPC* server, olong port, gchar *log_file)
{
	
  server->serverparm.port_number   = port;
  if (log_file) server->serverparm.log_file_name = log_file;
  else server->serverparm.log_file_name = "/tmp/xmlrpc_log";
  
  printf("Running XML-RPC server...port %d\n",port);

  /* Loop forever 'neath the streets of Boston */
  xmlrpc_server_abyss(&server->envP, &server->serverparm, 
		      XMLRPC_APSIZE(log_file_name));
  
} /* end ObitRPCServerLoop */

/** 
 * Run client async event loop
 * \param timeout  timeout in  microseconds, <= -> forever
 */
void  ObitRPCClientAsyncLoop (olong timeout)
{
  unsigned long RPCtimeout = (unsigned long)timeout;
  
  if (timeout>0) xmlrpc_client_event_loop_finish_asynch_timeout(RPCtimeout); 
  else xmlrpc_client_event_loop_finish_asynch();
} /* end ObitRPCClientAsyncLoop */



/**
 * Initialize global ClassInfo Structure.
 */
void ObitRPCClassInit (void)
{
  xmlrpc_env env;

  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitRPCClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */

  /* Initialize our error environment variable */
  xmlrpc_env_init(&env);

  /* Init global xmlrpc client stuff -
     Required before any use of Xmlrpc-c client library: */
  xmlrpc_client_setup_global_const(&env);
  die_if_fault_occurred(&env);
  
} /* end ObitRPCClassInit */

/**
 * Shutdown global ClassInfo Structure.
 */
void ObitRPCClassShutdown (void)
{

  if (!myClassInfo.initialized) return;  /* only once */
  
  myClassInfo.initialized = FALSE; /* Deinitialized */

  /* shutdown global xmlrpc client stuff */
  xmlrpc_client_teardown_global_const();
   
} /* end ObitRPCClassShutdown */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitRPCClassInfoDefFn (gpointer inClass)
{
  ObitRPCClassInfo *theClass = (ObitRPCClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */
  theClass->initialized = TRUE;  /* done here */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitRPCClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitRPCClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitRPCGetClass;
  theClass->newObit       = (newObitFP)newObitRPC;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitRPCClear;
  theClass->ObitInit      = (ObitInitFP)ObitRPCInit;
  theClass->ObitRPCCreateClient = 
    (ObitRPCCreateClientFP)ObitRPCCreateClient;
  theClass->ObitRPCCreateServer = 
    (ObitRPCCreateServerFP)ObitRPCCreateServer;
  theClass->numberClient  = 0;
  theClass->numberServer  = 0;

} /* end ObitRPCClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitRPCInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitRPC *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->registryP = NULL;
  in->clientP   = NULL;

  /* Initialize our error-handling environment. */
  xmlrpc_env_init(&in->envP);

} /* end ObitRPCInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitRPC* cast to an Obit*.
 */
void ObitRPCClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitRPC *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread  = ObitThreadUnref(in->thread);
  xmlrpc_env_clean(&in->envP);
  if (in->clientP) xmlrpc_client_destroy(in->clientP);
  in->clientP = NULL;

  /* Update type dependent stuff */
  if (in->type==OBIT_RPC_Client) {
    /*  NO myClassInfo.numberClient--;*/
    /* Shutdown our XML-RPC client library. */
    /*  NO if (myClassInfo.numberClient==0) xmlrpc_client_cleanup();*/
  }
  if (in->type==OBIT_RPC_Server) {
    /*  NO xmlrpc_registry_free(in->registryP);*/
    /* There seems to be no way to delete the server structures */
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitRPCClear */

