/* $Id:  $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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

#include <unistd.h>
#include "ObitMultiProc.h"
#include "ObitXML.h"
#include "ObitRPC.h"
#include "Obit.h"
#include "ObitReturn.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitMultiProc.c
 * ObitMultiProc class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitMultiProc";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitMultiProcClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitMultiProcClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/** MultiProc Thread routine argument for starting asynchronous processes */
typedef struct {
  /** Comamnd line  */
  const gchar *cmd_line;
  /** child status */
   gint status;
} ObitMultiProcAPFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitMultiProcInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitMultiProcClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitMultiProcClassInfoDefFn (gpointer inClass);

/** Private: Threaded multiprocessing send */
static gpointer MultiProcSnd (gpointer arg);

/** Private: Threaded multiprocessing receive */
static void MultiProcRcv (const char *server_url,
				const char *method_name,
				ObitXMLValue *param_array,
				void *user_data,
				ObitXMLEnv  *fault,
				ObitXMLValue *result);

/** Private: Split RPC URL into host and port */
static void splitURL (const gchar *RPCURL, gchar *host, olong *port);

/** Private: spawn asynchronous process */
static gboolean spawn (ObitThread *thread, const gchar *cmd_line);

/** Private: Threaded spawn asynchronous process */
static gpointer ThreadSpawn (gpointer arg);

/** Public: Runs multiple copies of a function in asynchronous processes */
static gboolean MultiProcExec (ObitThread *in, olong nstreams, 
			       ofloat timeout, ObitThreadFunc func, 
			       olong njobs, gboolean *doJob,
			       ObitMultiProcFuncArg *args[]);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitMultiProc* newObitMultiProc (gchar* name)
{
  ObitMultiProc* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitMultiProcClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitMultiProc));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitMultiProcInit((gpointer)out);

 return out;
} /* end newObitMultiProc */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitMultiProcGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitMultiProcClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitMultiProcGetClass */

/**
 * Make a deep copy of an ObitMultiProc.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitMultiProc* ObitMultiProcCopy  (ObitMultiProc *in, ObitMultiProc *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, oldNjobs;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitMultiProc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  oldNjobs       = out->njobs;
  out->njobs     = in->njobs;
  if (out->MPfunc) g_free(out->MPfunc);
  out->MPfunc    = g_strdup(in->MPfunc);
  out->localFunc = in->localFunc;
  if (out->doJob) g_free(out->doJob);
  out->doJob = g_malloc0(out->njobs*sizeof(gboolean));
  for (i=0; i<out->njobs; i++) out->doJob[i] = in->doJob[i];
  if (out->args) {
    for (i=0; i<oldNjobs; i++) {
      if (out->args[i]->client) ObitRPCUnref(out->args[i]->client);
      if (out->args[i]->URL)    g_free(out->args[i]->URL);
      if (out->args[i]->MPfunc) g_free(out->args[i]->MPfunc);
      if (out->args[i]->args)   ObitInfoListUnref(out->args[i]->args);
      if (out->args[i]->retVal) ObitInfoListUnref(out->args[i]->retVal);
      if (out->args[i]->thread) ObitThreadUnref(out->args[i]->thread);
      if (out->args[i]->err)    ObitErrUnref(out->args[i]->err);
      g_free(out->args[i]);
    }
    g_free(out->args);
  } /* End free old args */

  /* Create/copy args */
  out->args  = g_malloc0(out->njobs*sizeof(ObitMultiProcFuncArg*));
  for (i=0; i<out->njobs; i++) {
    out->args[i]  = g_malloc0(sizeof(ObitMultiProcFuncArg));
    out->args[i]->client    = ObitRPCRef(in->args[i]->client);
    if (in->args[i]->URL) 
      out->args[i]->URL     = g_strdup(in->args[i]->URL);
    out->args[i]->MPfunc    = g_strdup(in->args[i]->MPfunc);
    out->args[i]->localfunc = in->args[i]->localfunc;
    out->args[i]->args      = ObitInfoListRef(in->args[i]->args);
    out->args[i]->retVal    = newObitInfoList();
    out->args[i]->thread    = newObitThread();
    out->args[i]->ithread   = i;
    out->args[i]->err       = ObitErrRef(in->args[i]->err);
    out->args[i]->user_data = in->args[i]->user_data;
  } /* end init arrays */

  return out;
} /* end ObitMultiProcCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an MultiProc similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitMultiProcClone  (ObitMultiProc *in, ObitMultiProc *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, oldNjobs;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  oldNjobs       = out->njobs;
  out->njobs     = in->njobs;
  if (out->MPfunc) g_free(out->MPfunc);
  out->MPfunc    = g_strdup(in->MPfunc);
  out->localFunc = in->localFunc;
  if (out->doJob) g_free(out->doJob);
  out->doJob = g_malloc0(out->njobs*sizeof(gboolean));
  for (i=0; i<out->njobs; i++) out->doJob[i] = in->doJob[i];
  if (out->args) {
    for (i=0; i<oldNjobs; i++) {
      if (out->args[i]->client) ObitRPCUnref(out->args[i]->client);
      if (out->args[i]->URL)    g_free(out->args[i]->URL);
      if (out->args[i]->MPfunc) g_free(out->args[i]->MPfunc);
      if (out->args[i]->args)   ObitInfoListUnref(out->args[i]->args);
      if (out->args[i]->retVal) ObitInfoListUnref(out->args[i]->retVal);
      if (out->args[i]->thread) ObitThreadUnref(out->args[i]->thread);
      if (out->args[i]->err)    ObitErrUnref(out->args[i]->err);
      g_free(out->args[i]);
    }
    g_free(out->args);
  } /* End free old args */

  /* Create/copy args */
  out->args  = g_malloc0(out->njobs*sizeof(ObitMultiProcFuncArg*));
  for (i=0; i<out->njobs; i++) {
    out->args[i]  = g_malloc0(sizeof(ObitMultiProcFuncArg));
    out->args[i]->client    = ObitRPCRef(in->args[i]->client);
     if (in->args[i]->URL) 
      out->args[i]->URL       = g_strdup(in->args[i]->URL);
    out->args[i]->MPfunc    = g_strdup(in->args[i]->MPfunc);
    out->args[i]->localfunc = in->args[i]->localfunc;
    out->args[i]->args      = ObitInfoListRef(in->args[i]->args);
    out->args[i]->retVal    = newObitInfoList();
    out->args[i]->thread    = newObitThread();
    out->args[i]->ithread   = i;
    out->args[i]->err       = ObitErrRef(in->args[i]->err);
    out->args[i]->user_data = in->args[i]->user_data;
  } /* end init arrays */


} /* end ObitMultiProcClone */

/**
 * Creates an ObitMultiProc 
 * \param name      An optional name for the object.
 * \param njobs     Number of "jobs" to be processed
 * \param MPfunc    Name of remote FuncContainer function
 * \param localFunc Local pointer to job function which should be the
 *                  same as that executed in FuncContainer.
 * \return the new object.
 */
ObitMultiProc* ObitMultiProcCreate (gchar* name, olong njobs,
				    gchar *MPfunc, ObitMultiProcFunc localFunc,
				    ObitErr *err)
{
  ObitMultiProc* out;
  ObitRPC *client=NULL;
  olong i;
  gchar *routine = "ObitMultiProcCreate";

  /* Create basic structure */
  out = newObitMultiProc (name);

  /* Initialize */
  ObitMultiProcInit (out);

  /* Set members */
  out->njobs = njobs;
  out->MPfunc = g_strdup(MPfunc);
  out->localFunc = localFunc;

  /* Arrays */
  out->doJob = g_malloc0(out->njobs*sizeof(gboolean));
  out->args  = g_malloc0(out->njobs*sizeof(ObitMultiProcFuncArg*));

  /* ObitRPC client for all */
  client =  ObitRPCCreateClient ("RPCclient", err);
  if (err->error) Obit_traceback_val (err, routine, "Create", out); 

  /* Initialize arrays */
  for (i=0; i<out->njobs; i++) {
    out->doJob[i] = TRUE;   /* Until indicated otherwise */
    out->args[i]  = g_malloc0(sizeof(ObitMultiProcFuncArg));
    out->args[i]->client    = ObitRPCRef(client);
    out->args[i]->URL       = NULL;
    out->args[i]->MPfunc    = g_strdup(MPfunc);
    out->args[i]->localfunc = localFunc;
    out->args[i]->args      = newObitInfoList();
    out->args[i]->retVal    = NULL;
    out->args[i]->thread    = ObitThreadRef(out->thread);
    out->args[i]->ithread   = i;
    out->args[i]->user_data = NULL;
  } /* end init arrays */

   /*ObitRPCUnref(client);  Cleanup */
  return out;
} /* end ObitMultiProcCreate */

/**
 * Enables multiprocessing based on values in myInput and the availability
 * of threading. 
 * Starts auxillary executions of FuncContainer and saves the necessary information
 * in class structures where they are available to any ObitMultiProc instance.
 * \param myInput Task parameter input list.  Significant entries:
 * \li "RPCURL"  If present and has more than one nonblank entry and RPCURL has
 *                no non-blank entries then a set processes running FuncContainer
 *                will be started the on each URL watching each port in the URL..
 * \li "RPCports" If present and has more than one positive entry and RPCURL has
 *                no non-blank entries then a set processes running FuncContainer
 *                will be started on localhost watching each port listed.
 * \li "AIPSuser" AIPS user number
 * \li "nAIPS"    Number of AIPS data directories
 * \li "AIPSdirs" List of AIPS data directories
 * \li "nFITS"    Number of FITS data directories
 * \li "FITSdirs" List of FITS data directories
 * \li "nThreads" Number of threads functions are allowed to use.
 * \param err Obit error stack object.
 */
void ObitMultiProcStart (ObitInfoList *myInput, ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  gchar *strTemp;
  oint *itemp, i, j, nURL;
  olong pgmNumber, port;
  ObitInfoList *taskList=NULL;
  gchar tempURL[300], cmd_line[300], host[300], taskFile[300];
  gchar        *taskParms[] = {  /* task parameters */
    "AIPSuser", "nAIPS", "AIPSdirs", "nFITS", "FITSdirs", 
    "nThreads", "noScrat", "prtLv",
    NULL
  };
  gchar *routine = "ObitMultiProcStart";
  
  if (err->error) return;
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitMultiProcClassInit();
  
  /* Check for RPCURL */
  nURL = 0;
  ObitInfoListGetP(myInput, "RPCURL",  &type, dim, (gpointer)&strTemp);
  if (strTemp!=NULL) {
    for (i=0; i<dim[1]; i++) { /* loop over array */
      /* Something here? */     
      if (strTemp[0]!=' ') {
	if (!myClassInfo.URL) myClassInfo.URL = g_malloc0(dim[1]*sizeof(gchar*));
	for (j=0; j<dim[0]; j++) {
	  if (strTemp[j]==' ') break;  /* drop blanks */
	  tempURL[j] = strTemp[j];
	}
	tempURL[j] = 0;
	myClassInfo.URL[nURL] = g_strdup(tempURL);
	nURL++;
	strTemp += dim[0];
      } else break;
    }
  } /* end if RPCURL present */
  
  if (nURL<=0) {
    ObitInfoListGetP(myInput, "RPCport",  &type, dim, (gpointer)&itemp);
    if (itemp!=NULL) {
      for (i=0; i<dim[0]; i++) { /* loop over array */
	if (itemp[i]<=0) break;
	if (!myClassInfo.URL) myClassInfo.URL = g_malloc0(dim[0]*sizeof(gchar*));
	sprintf (tempURL,"http://localhost:%d/RPC2",itemp[i]);
	myClassInfo.URL[nURL] = g_strdup(tempURL);
 	nURL++;
     }
    }
  } /* End get URL from local port number */
  
  if (nURL<=0) return;  /* Nothing to do */
  
  myClassInfo.nProcessor = nURL;  /* save Number found */
  
  /* Create basic InfoList to pass to processes */
  taskList = newObitInfoList();
  ObitInfoListCopyList (myInput, taskList, taskParms);   

  /* Create threads for running  processes */
  myClassInfo.ASThreads = g_malloc0(nURL*sizeof(ObitThread*));
  for (i=0; i<nURL; i++) myClassInfo.ASThreads[i] = newObitThread();
  
  /* Place to save command lines */
  myClassInfo.cmd_lines = g_malloc0(nURL*sizeof(gchar*));

  for (i=0; i<nURL; i++) { /* loop over processes */
    pgmNumber = i+1; /* Program number */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (taskList, "pgmNumber", OBIT_int, dim, &pgmNumber);
    
    /* Get port and host URL */
    splitURL (myClassInfo.URL[i], host, &port);
    ObitInfoListAlwaysPut (taskList, "port", OBIT_int, dim, &port);
    /* write parameter file */
    sprintf (taskFile, "/tmp/FuncContainerInput%d.inp", pgmNumber);
    ObitReturnDumpRetCode (0, taskFile, taskList, err);  
    if (err->error) Obit_traceback_msg (err, routine, "Starting MultiProcs");

    /* set command line */
    sprintf (cmd_line, "FuncContainer -input %s -port %d", taskFile, port);
    myClassInfo.cmd_lines[i] = g_strdup(cmd_line);
    if (!spawn (myClassInfo.ASThreads[i], myClassInfo.cmd_lines[i])) {
      Obit_log_error(err, OBIT_InfoWarn, "Error starting process port %d on %s - already running?", 
		     port, host);
    }
  } /* end loop over processes */
  
  /* If successful */
  myClassInfo.haveAsynchProc = TRUE;  
  ObitErrLog(err); /* show any error messages on err */
  
} /* end ObitMultiProcStart */

/**
 * Disables multiprocessing in auxillary FuncContainer processes.
 * Sequential processing using the ObitMultiProc will still be supported.
 * \param err Obit error stack object.
 * \return NULL
 */
void ObitMultiProcShutdown (ObitErr *err)
{ 
  olong i;
  ObitXML *arg=NULL, *reply=NULL;
  ObitInfoList *argList=NULL;
  ObitRPC *client=NULL;
  gchar *routine = "ObitMultiProcShutdown";
  
  /* Class initialization if needed - then clearly nothing to shutdown */
  if (!myClassInfo.initialized) {
    ObitMultiProcClassInit();
    return;
  }
  
  /* ObitRPC client for all */
  client =  ObitRPCCreateClient ("ShutdownClient", err);
  if (err->error) goto done;

  /* Send processes the "terminate" command */
  for (i=0; i<myClassInfo.nProcessor; i++) {
    
    /* Create argument */
    argList = newObitInfoList();
    arg = ObitXMLSetCallArg ("terminate", argList, err);
    if (err->error) goto done;

    /* termination call */
    reply = ObitRPCCall (client, myClassInfo.URL[i], arg, NULL, NULL, err);
    if (err->error) goto done;

    /* cleanup */
    argList = ObitInfoListUnref(argList);
    reply   = ObitXMLUnref(reply);
    arg     = ObitXMLUnref(arg);
    myClassInfo.ASThreads[i] = ObitThreadUnref(myClassInfo.ASThreads[i]);
    if (myClassInfo.URL[i]) g_free(myClassInfo.URL[i]);
  } /* end loop over remote processes */

 done:
   ObitRPCUnref(client);  /* Cleanup */
  
  /* No longer allow use of asynchronous processes */
  myClassInfo.haveAsynchProc = FALSE; 

  /* Something go wrong? */
  if (err->error) Obit_traceback_msg (err, routine, routine); 
  
} /* end ObitMultiProcShutdown */

/**
 * Copies argument values for a given job
 * of threading 
 * \param in     MultiProc object
 * \param jobNo  0-rel job number 
 * \param arg    Argument list to copy
 */
void ObitMultiProcSetFuncArg(ObitMultiProc* in, olong jobNo, ObitInfoList *arg)
{ 
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if ((jobNo<0) || (jobNo>=in->njobs)) return;

  in->args[jobNo]->args = ObitInfoListCopyData (arg, in->args[jobNo]->args);
} /* end ObitMultiProcSetFuncArg */

/**
 * Set execution flag for a given job
 * \param in     MultiProc object
 * \param jobNo  0-rel job number 
 * \param flag   If TRUE (default) this this job will be executed by
 *               ObitMultiProcExecute, otherwise, not.
 */
/* Set execute flag for a given job */
void ObitMultiProcSetExecFlag(ObitMultiProc* in, olong jobNo, gboolean flag)
{ 
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if ((jobNo<0) || (jobNo>=in->njobs)) return;
 
  in->doJob[jobNo] = flag;
} /* end ObitMultiProcSetExecFlag */

/**
 * Execute selected jobs
 * Can run synchronously or asynchronously if enabled.
 * \param in      MultiProc object to execute
 * \param timeout Timeout (min) per function call, <= 0 -> forever 
 * \param err     Obit error/message structure to be used for process
 *                messages and error handling
 */
void ObitMultiProcExecute (ObitMultiProc* in, ofloat timeout, ObitErr *err)
{ 
  olong i, nstreams, ndo;
  gboolean OK;
  gchar *routine = "ObitMultiProcExecute";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* Set err members of args to err */
    for (i=0; i<in->njobs; i++) {
      if (in->args[i]->err) ObitErrUnref(in->args[i]->err);
      in->args[i]->err = ObitErrRef(err);
    }

  /* Is MultiProcessing with asynchronous processes enabled? */
  if (myClassInfo.haveAsynchProc) { /* do asynchronously */

    /* Assign URLs to first nstreams */
    nstreams = myClassInfo.nProcessor;
    ndo = 0;
    for (i=0; i<in->njobs; i++) {
      if (in->doJob[i]) {
	if (in->args[i]->URL) g_free(in->args[i]->URL);
	in->args[i]->URL = g_strdup(myClassInfo.URL[ndo]);
	ndo++;
	if (ndo>= nstreams) break;
      }
    } /* end asigning URLs */
    

    /* Execute */
    OK = MultiProcExec (in->thread, nstreams, timeout, (ObitThreadFunc)MultiProcSnd, 
			in->njobs, in->doJob, in->args);
    if (!OK) {  /* Something went wrong */
      /* Error checks */
      for (i=0; i<in->njobs; i++) {
	if ((err->error) || (in->args[i]->retVal==NULL))
	  Obit_traceback_msg (err, routine, in->name);  
      }
    }
    
  } else {  /* do synchronously */
    for (i=0; i<in->njobs; i++) {
      if (in->doJob[i]) {
	in->args[i]->retVal = in->localFunc (in->args[i]->args, err);
	/* Error check */
 	if ((err->error) || (in->args[i]->retVal==NULL))
	  Obit_traceback_msg (err, routine, in->name);  
     }
    }
  } /* end synchronous/asynchronous divide */
 
} /* end ObitMultiProcExecute */

/**
 * Return reference to job return value ObitInfoList
 * \param  in     MultiProc object
 * \param  jobNo  0-rel job number 
 * \return pointer to return value ObitInfoList, NULL on error
 */
ObitInfoList* ObitMultiProcGetFuncRet(ObitMultiProc* in, olong jobNo)
{ 
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
 
  if ((jobNo<0) || (jobNo>=in->njobs)) return NULL;
  return ObitInfoListRef(in->args[jobNo]->retVal);
} /* end ObitMultiProcGetFuncRet */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitMultiProcClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitMultiProcClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitMultiProcClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitMultiProcClassInfoDefFn (gpointer inClass)
{
  ObitMultiProcClassInfo *theClass = (ObitMultiProcClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitMultiProcClassInit;
  theClass->newObit       = (newObitFP)newObitMultiProc;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitMultiProcClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitMultiProcGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitMultiProcCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitMultiProcClear;
  theClass->ObitInit      = (ObitInitFP)ObitMultiProcInit;
  theClass->ObitMultiProcCreate = (ObitMultiProcCreateFP)ObitMultiProcCreate;

  /* Other class data */
  theClass->haveAsynchProc = FALSE;
  theClass->nProcessor     = 0;
  theClass->URL            = NULL;

} /* end ObitMultiProcClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitMultiProcInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMultiProc *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread = newObitThread();
  in->njobs     = 0;
  in->MPfunc    = NULL;
  in->localFunc = NULL;
  in->args      = NULL;
  in->doJob     = NULL;

} /* end ObitMultiProcInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitMultiProc* cast to an Obit*.
 */
void ObitMultiProcClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitMultiProc *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->thread) ObitThreadUnref(in->thread);
  if (in->MPfunc) g_free(in->MPfunc);
  if (in->doJob)  g_free(in->doJob);
  if (in->args) {
    for (i=0; i<in->njobs; i++) {
      if (in->args[i]->client) ObitRPCUnref(in->args[i]->client);
      if (in->args[i]->URL)    g_free(in->args[i]->URL);
      if (in->args[i]->MPfunc) g_free(in->args[i]->MPfunc);
      if (in->args[i]->args)   ObitInfoListUnref(in->args[i]->args);
      if (in->args[i]->retVal) ObitInfoListUnref(in->args[i]->retVal);
      if (in->args[i]->thread) ObitThreadUnref(in->args[i]->thread);
      if (in->args[i]->err)    ObitErrUnref(in->args[i]->err);
      g_free(in->args[i]);
    }
    g_free(in->args);
  } /* End free args */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitMultiProcClear */

/** 
 * Multiprocessing call passed through thread
 * This will forward the processing request to the FuncContainer
 * asynchronously with arguments specified in args.
 * Arguments are given in the structure passed as arg
 * \param arg  Pointer to ObitMultiProcFuncArg argument with elements:
 * \li client     RPC client
 * \li URL        URL for FuncContainer
 * \li MPfunc     Function name to call in FuncContainer for remote RPC call
 * \li localfunc  Function pointer to call in local address space (unused here)
 * \li args       Arguments as InfoList 
 * \li retVal     Return values as InfoList
 * \li thread     thread object
 * \li ithread    thread number
 * \li err        Obit error/logging object
 *  \return NULL
 */
static gpointer MultiProcSnd (gpointer arg)
{
  /* Get arguments from structure */
  ObitMultiProcFuncArg *largs = (ObitMultiProcFuncArg*)arg;
  ObitRPC *client      = largs->client;
  gchar *URL           = largs->URL;
  gchar *MPfunc        = largs->MPfunc;
  ObitInfoList *args   = largs->args;
  ObitThread *thread   = largs->thread;
  olong ithread        = largs->ithread;
  ObitErr *err         = largs->err;

  ObitXML *XMLarg=NULL;
  gboolean doCall;
  ObitErr *terr = newObitErr();
  gchar *routine = "ThreadMultProcSnd";
 
  doCall = FALSE;
   
  /* Create XML argument */
  XMLarg = ObitXMLSetCallArg(MPfunc, args, terr);
  if (terr->error) goto done;
  
  /* Do call  */
  ObitRPCCallSnd (client, URL, XMLarg, 
		  (ObitRPC_response_handler)MultiProcRcv, 
		  (gpointer)arg, terr);
   if (terr->error) goto done;
   else   doCall = TRUE;

  ObitThreadUnlock(thread);

  /* finished */
 done:
  
  if (terr->error) {
    ObitErrLog(err); /* show any error messages on err */
    Obit_log_error(terr, OBIT_Error,   "%s %d: Error making RPC call", routine, ithread);
}

  /* Cleanup */
  if (XMLarg) ObitXMLUnref(XMLarg);

  /* If something went wrong indicate completion */
  if (!doCall)  ObitThreadPoolDone (thread, (gpointer)&URL);
  return NULL;

} /* end MultiProcSnd */

/** 
 * Multiprocessing call passed through thread
 * This callback function processed the result of the resuest.
 * Logging messages are extracted from the return InfoList
 * and added to the arg-err member.
 * If the return value is other than "OK" and error is indicated in arg->err.
 * Arguments are given in the structure passed as arg
 * At completion, sends URL of the FuncContainer used to ObitThreadPoolDone.
 * \param server_url   Server called
 * \param method_name  Function called 
 * \param param_array  Request sent
 * \param user_data  Pointer to ObitMultiProcFuncArg argument with elements:
 * \li client     RPC client
 * \li URL        URL for FuncContainer
 * \li MPfunc     Function name to call in FuncContainer for remote RPC call
 * \li localfunc  Function pointer to call in local address space (unused here)
 * \li args       Arguments as InfoList 
 * \li retVal     Return values as InfoList
 * \li thread     thread object
 * \li ithread    thread number
 * \li err        Obit error/logging object
 * \li user_data  Unused
 * \param fault   RPC environment
 * \param result  Result of RPC call
 */
static void MultiProcRcv (const char *server_url,
				const char *method_name,
				ObitXMLValue *param_array,
				void *user_data,
				ObitXMLEnv  *fault,
				ObitXMLValue *result)
{
  /* Get arguments from structure */
  ObitMultiProcFuncArg *largs = (ObitMultiProcFuncArg*)user_data;
  ObitRPC *client      = largs->client;
  gchar *URL           = largs->URL;
  ObitThread *thread   = largs->thread;
  olong ithread        = largs->ithread;
  ObitErr *err         = largs->err;

  ObitXML *XMLreply=NULL;
  ObitInfoType type;
  olong nlog, ilog;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  gchar *strTemp, mlabel[40];
  ObitErrCode ECode;
  ObitErr *terr = newObitErr();
  gchar *routine = "ThreadMultProcRcv";

  /* Convert reply to XML object */
  XMLreply = ObitRPCCallRcv (client, result, NULL, NULL, terr);

  /* Convert reply to InfoList */
  largs->retVal = ObitXMLGetServerResult(XMLreply, terr);
  if (terr->error) goto done;

  /* Copy log to err */
  if (ObitInfoListGetTest (largs->retVal, "NumMessage",  &type, dim, &nlog)) {
    ObitThreadLock(thread); 
    for (ilog=1; ilog<=nlog; ilog++) {
      sprintf (mlabel,"M%8.8d",ilog);
      ObitInfoListGetP(largs->retVal, mlabel,  &type, dim, (gpointer)&strTemp);
      if (strTemp==NULL) break;
      /* Convert beginning into ObitErr code */
      if (*strTemp=='i') ECode = OBIT_InfoErr;
      else if (*strTemp=='w') ECode = OBIT_InfoWarn;
      else if (*strTemp=='t') ECode = OBIT_Traceback;
      else if (*strTemp=='E') ECode = OBIT_Error;
      else if (*strTemp=='M') ECode = OBIT_MildError;
      else if (*strTemp=='S') ECode = OBIT_StrongError;
      else if (*strTemp=='F') ECode = OBIT_Fatal;
      else  ECode = OBIT_None;
      Obit_log_error(err, ECode, "%s", &strTemp[12]);
    }
    ObitThreadUnlock(thread);
  } /* End if mesages */
 
  /* Check Status for error */
   ObitInfoListGetP(largs->retVal, "Status",  &type, dim, (gpointer)&strTemp);
   if ((strTemp==NULL) || (strTemp[0]!='O') || (strTemp[1]!='K') || (strTemp[2]!=0)) {
     /* Oh Bother! it failed */
     ObitThreadLock(thread);
     Obit_log_error(err, OBIT_Error,   "%s %d: Error in RPC call, status %s", routine, ithread, strTemp);
     ObitThreadUnlock(thread);
  }
 
  /* finished */
 done:
  
  /* Cleanup */
  if (XMLreply)  ObitXMLUnref(XMLreply);
  if (client)    ObitRPCUnref(client);

  /* Indicate completion */
  ObitThreadPoolDone (thread, (gpointer)URL);
  
} /* end MultiProcRcv */

/**
 * Extract host URL and port number from RPCURL
 * in class structures where they are available to any ObitMultiProc instance.
 * \param RPCURL   input RPC URL, e.g. "http://localhost:8766/RPC2"
 * \param hostURL  host   e.g. localhost
 * \param port     Port number e.g. 8766
 */
static void splitURL (const gchar *RPCURL, gchar *host, olong *port)
{
  olong i, j, n, doubleslash=0, colon=0, singleslash=0;
  gchar temp[20];

  /* find double slash */
  n = strlen (RPCURL);
  for (i=1; i<n; i++) {
    if ((RPCURL[i-1]=='/')&&(RPCURL[i]=='/')) {doubleslash = i; break;}
  }

  /* find colon */
  for (i=doubleslash; i<n; i++) {
    if (RPCURL[i]==':') {colon = i; break;}
  }

  /* find single slash after colon */
  for (i=colon; i<n; i++) {
    if (RPCURL[i]=='/') {singleslash = i; break;}
  }
    
  /* Extract host */
  j = 0;
  for (i=doubleslash+1; i<colon; i++) host[j++] = RPCURL[i];
  host[j] = 0;
  
  /* Extract port */
  j = 0;
  for (i=colon+1; i<singleslash; i++) temp[j++] = RPCURL[i];
  temp[j] = 0;
  *port = (olong)strtol(temp, NULL, 0);
  
  }  /* end splitURL */

/**
 * spawn asynchronous process 
 * \param thread to use
 * \param cmd_line command line to execute
 */
static gboolean spawn (ObitThread *thread, const gchar *cmd_line)
{
  ObitMultiProcAPFuncArg arg;
  arg.cmd_line = cmd_line;
  arg.status   = -999;
  

  ObitThreadStart1 (thread, (ObitThreadFunc)ThreadSpawn, (gpointer)&arg);
  /* Give this some time to start up */
  usleep(500000); /* 500 msec */
  return (arg.status ==-999);
} /* end spawn */

/**
 * Thread spawn asynchronous process 
 * glib will run this routine in a separate thread
 * \param argument as ObitMultiProcAPFuncArg
 * \return NULL
 */
static gpointer ThreadSpawn (gpointer arg)
{
  ObitMultiProcAPFuncArg *larg = (ObitMultiProcAPFuncArg*)arg;
  gboolean OK;

  OK = g_spawn_command_line_sync (larg->cmd_line, NULL, NULL, &larg->status, NULL);
  fprintf (stderr, "Starting asynchronous process failed\n");

  return NULL;
} /* end spawn */

/**
 * Loops over a set of function calls in a limited number of job streams.
 * The functions are intended to contain an ObitRPCCall to an asynchronous 
 * process.  The URL/port numbers of the first nstreams elements of args
 * are recycled as previous jobs are completed.
 * If nstreams=1 routine called directly.
 * Waits for operations to finish before returning, 1 min timeout.
 * Initializes asynchronous queue (ObitThreadPoolInit) if not already done.
 * When threaded operations are finished, call ObitThreadPoolFree to release
 * Thread pool.
 * \param thread     Pointer to thread object
 * \param nstreams   Number of threads to create/run
 * \param timeout    Timeout (min) per function call, <= 0 -> forever NYI
 * \param func       Function to call to start thread
 *                   func should call ObitThreadPoolDone to indicate completion 
 *                   and return the URL used for the RPC operation.
 * \param njobs      Number of executions of func each with a separate args[i].
 * \param doJob      Array of flags indicating which of njobs are desired
 * \param args       Array of argument function pointers, njob elements
 * \return TRUE if OK else FALSE.  If false, check the err elements in args.
 */
static gboolean MultiProcExec (ObitThread* thread, olong nstreams, 
			       ofloat timeout, ObitThreadFunc func, 
			       olong njobs, gboolean *doJob,
			       ObitMultiProcFuncArg *args[])
{
  gboolean out = TRUE;
  glong add_time;
  olong i, j, k, nsubmit, isub, ndo, ndone=0;
  gpointer retVal;
  gchar *URL;
  olong RPCtimeout;
  gchar *routine = "MultiProcExec";

  /* error checks */
  g_assert (thread != NULL);

  /* How many? */
  ndo = 0;
  for (i=0; i<njobs; i++) if (doJob[i]) ndo++;
  if (ndo<=0) return TRUE;

  /* If only 1 don't use thread run locally */
  if (njobs==1)  {
    args[0]->retVal = args[0]->localfunc (args[0]->args , args[0]->err);
    if (args[0]->retVal==NULL) out = FALSE;
    if (args[0]->err->error)   {
      out = FALSE;
      Obit_traceback_val (args[0]->err, routine, "execution", out);  
    }
    return out;
  }

  /* If no threading don't use thread - call sequentially */
  if ((ndo<=1) || (!ObitThreadHaveThreads(thread))) {
    for (i=0; i<njobs; i++) {
      if (doJob[i]) {
        args[i]->retVal = args[0]->localfunc (args[i]->args, args[i]->err);
        if (args[i]->retVal==NULL) out = FALSE;
        if (args[i]->err->error)   {
          out = FALSE;
          Obit_traceback_val (args[i]->err, routine, "execution", out);
        }
      }
    }
    return out;
  }

  /* Submit the first nstreams jobs - asynchronous RPC calls */
  nsubmit = 0;
  isub    = 0;
  for (i=0; i<nstreams; i++) {
    if (doJob[i]) {
      func (args[i]);
      isub = i;
      nsubmit++;
    }
  }
 
  /* Start asynchronous message queue */
  ObitThreadQueueInit (thread);

  /* Wait for them to finish, expects each to send a message to the asynchronous 
     queue iff they finish giving the URL actually used */
  /* set timeouts */
  RPCtimeout = 50;            /* RPC event queue timeout 50 msec */
  add_time = 0.05 * 1000000;  /* async_queue timeout 50 msec */
  while (ndone<njobs) {  /* Loop until all done */
    
    /* Impatiently wait for one to finish - run RPC event queue 
       then check async_queue */
    if (nsubmit>=ndo) RPCtimeout = 0;   /* If Done - wait for all */
    ObitRPCClientAsyncLoop (RPCtimeout);                     /* RPC */
    retVal = ObitThreadQueueCheck(thread, add_time);         /* Completion queue */
    
    /* if retVal is nonNULL, something finished */
    if (retVal) ndone++;
    else continue;
      
    /* More to submit? */
    if (retVal && (nsubmit<ndo)) {
      /* find next to submit */
      for (j=isub+1; j<njobs; j++) {
	if (doJob[j]) {
	  isub = j;
	  break;
	}
      }
      /* Recycle URL */
      URL = (gchar*)retVal;
      if (args[isub]->URL) g_free(args[isub]->URL);
      args[isub]->URL = g_strdup(URL);

      /* Stick next into the queue */
      func (args[isub]);
     nsubmit++;
    } /* End finding next to submit */
    
    /* error check on jobs already finished */
    for (k=0; k<i; k++) {
      /* Only want to wait on the actual number submitted */
      if (doJob[k]) {
	if (args[k]->retVal==NULL) out = FALSE;
	if (args[k]->err->error)   {
	  out = FALSE;
	  Obit_traceback_val (args[k]->err, routine, "execution", out);  
	}
      } /* end if doJob */
    } /* end loop over jobs error check */
  } /* end loop till all jobs done*/

  /* Start asynchronous message queue */
  ObitThreadQueueFree (thread);

  return out;
} /* end MultiProcExec */


