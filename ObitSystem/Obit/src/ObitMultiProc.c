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
  olong i;
  /*gchar *routine = "ObitMultiProcCreate";*/

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
  
  /* Initialize arrays */
  for (i=0; i<out->njobs; i++) {
    out->doJob[i] = TRUE;   /* Until indicated otherwise */
    out->args[i]  = g_malloc0(sizeof(ObitMultiProcFuncArg));
    out->args[i]->client    = NULL;
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
  gboolean OK, local=FALSE;
  ObitFile *parmFile=NULL;
  gint status;
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
	if (!myClassInfo.clients) myClassInfo.clients = g_malloc0(dim[1]*sizeof(ObitRPC*));
	for (j=0; j<dim[0]; j++) {
	  if (strTemp[j]==' ') break;  /* drop blanks */
	  tempURL[j] = strTemp[j];
	}
	tempURL[j] = 0;
	myClassInfo.URL[nURL]     = g_strdup(tempURL);
	myClassInfo.clients[nURL] = ObitRPCCreateClient ("RPCclient", err);
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
	if (!myClassInfo.clients) myClassInfo.clients = g_malloc0(dim[1]*sizeof(ObitRPC*));
	sprintf (tempURL,"http://localhost:%d/RPC2",itemp[i]);
	myClassInfo.URL[nURL]     = g_strdup(tempURL);
	myClassInfo.clients[nURL] = ObitRPCCreateClient ("RPCclient", err);
 	nURL++;
     }
    }
  } /* End get URL from local port number */
  
  if (err->error) Obit_traceback_msg (err, routine, "Starting MultiProcs");
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

    /* Local or remote? */
    local = !strcmp(host, "localhost");
    /* write parameter file */
    if (local) {  /* Local */
      sprintf (taskFile, "/tmp/FuncContainerInput%d.inp", pgmNumber);
      ObitReturnDumpRetCode (0, taskFile, taskList, err);  
      if (err->error) Obit_traceback_msg (err, routine, "Starting writing FuncContainer input");
    } else { /* remote */
      sprintf (taskFile, "/tmp/FuncContainerInput%d.inp", pgmNumber);
      ObitReturnDumpRetCode (0, taskFile, taskList, err);  
      if (err->error) Obit_traceback_msg (err, routine, "Writing FuncContainer input");
      sprintf (cmd_line, "scp -q %s %s:/tmp/FuncContainerInput%d.tmp", taskFile, host, pgmNumber);
      /* DEBUG */
      fprintf (stderr, "%s\n", cmd_line);
      OK = g_spawn_command_line_sync (cmd_line, NULL, NULL, &status, NULL);
    } /* end write task parameters */

    /* set command line */
    if (local) {  /* Local */
      sprintf (cmd_line, "FuncContainer -input %s -port %d", taskFile, port);
      myClassInfo.cmd_lines[i] = g_strdup(cmd_line);
      if (!spawn (myClassInfo.ASThreads[i], myClassInfo.cmd_lines[i])) {
	Obit_log_error(err, OBIT_InfoWarn, "Error starting process port %d on %s - already running?", 
		       port, host);
      }
    } else { /* remote */
      sprintf (taskFile, "/tmp/FuncContainerInput%d.tmp", pgmNumber);
      sprintf (cmd_line, "ssh %s /tmp/FuncContainer -input %s -port %d", host, taskFile, port);
      /* DEBUG */
      fprintf (stderr, "%s\n", cmd_line);
      myClassInfo.cmd_lines[i] = g_strdup(cmd_line);
      if (!spawn (myClassInfo.ASThreads[i], myClassInfo.cmd_lines[i])) {
	Obit_log_error(err, OBIT_InfoWarn, "Error starting process port %d on %s - already running?", 
		       port, host);
      }
    } /* end start FuncContainer */
  } /* end loop over processes */

  /* If successful */
  myClassInfo.haveAsynchProc = TRUE;  
  ObitErrLog(err); /* show any error messages on err */

  /* Wait a bit for things to actually start */
  usleep(2000000);  /* 2 sec */
  
  /* delete temporary (remote) input files */
  for (i=0; i<nURL; i++) { /* loop over processes */
    pgmNumber = i+1; /* Program number */
    /* Local or remote? */
    local = !strcmp(host, "localhost");
    if (!local) { /* remote */
      sprintf (taskFile, "/tmp/FuncContainerInput%d.inp", pgmNumber);
      /* delete temp taskFile */
      parmFile = newObitFile(taskFile);
      ObitFileOpen (parmFile, taskFile, OBIT_IO_ReadWrite, OBIT_IO_Text, 0, err);
      ObitFileClose(parmFile, err);
      /*ObitFileZap(parmFile, err); DEBUG*/
      if (err->error) Obit_traceback_msg (err, routine, "Zapping FuncContainer input");
    }
  }
  
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
 * If the URLs are specified in in->args[*]->URL, they are used, otherwise the
 * first myClassInfo.nProcessor are assigned the list of URL in the class global structure.
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

    /* Assign URL/RPC clientss to first nstreams if not given */
    nstreams = myClassInfo.nProcessor;
    ndo = 0;
    for (i=0; i<in->njobs; i++) {
      if ((in->doJob[i]) && (in->args[i]->URL==NULL)) {
	in->args[i]->URL    = g_strdup(myClassInfo.URL[ndo]);
	in->args[i]->client = ObitRPCRef(myClassInfo.clients[ndo]);
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
	Obit_return_if_fail(((!err->error) && (in->args[i]->retVal!=NULL)), 
			    err, "%s: MultiProcExec RPC Call %d failed", routine, i);
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
 * with arguments specified in args.
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
 * \li ijob       job number
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

  ObitXML *XMLarg=NULL, *XMLreply=NULL;
  ObitErr *terr = newObitErr();
  ObitInfoType type;
  olong nlog, ilog;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  gchar *strTemp, mlabel[40];
  ObitErrCode ECode;
  gchar *routine = "MultProcSnd";
  
  /* Create XML argument */
  XMLarg = ObitXMLSetCallArg(MPfunc, args, terr);
  if (terr->error) goto done;
  
  /* Do call  */
  XMLreply = ObitRPCCall (client, URL, XMLarg, NULL, NULL, terr);
  if (terr->error) goto done;
  
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
     Obit_log_error(err, OBIT_Error,   "%s %d: Error in RPC call, status %s", \
		    routine, ithread, strTemp);
     ObitThreadUnlock(thread);
  }
 
  /* finished */
 done:
  
   if (terr->error) {
     ObitThreadLock(thread);
     ObitErrLog(terr); /* show any error messages on terr */
     Obit_log_error(err, OBIT_Error,   "%s %d: Error making RPC call", routine, ithread);
     ObitThreadUnlock(thread);
   }

   /* Cleanup */
   if (XMLarg)   ObitXMLUnref(XMLarg);
   if (XMLreply) ObitXMLUnref(XMLreply);
   
  /* Indicate completion */
  ObitThreadPoolDone (thread, (gpointer)&largs->ijob);
   return NULL;
   
} /* end MultiProcSnd */

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
 * should be set.  If subsequent jobs must be performed on a particular URL,
 * their URLs should be set, otherwise the next available URL will be used.
 * If nstreams=1 routine called directly.
 * Waits for operations to finish before returning,
 * Initializes asynchronous queue (ObitThreadPoolInit) if not already done.
 * When threaded operations are finished, call ObitThreadPoolFree to release
 * Thread pool.
 * \param thread     Pointer to thread object
 * \param nstreams   Number of threads to create/run
 * \param timeout    Timeout (min) per function call, <= 0 -> forever NYI
 * \param func       Function to call to start thread
 *                   func should call ObitThreadPoolDone to indicate completion 
 *                   and return the job number used for the RPC operation.
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
  olong i, j, k, nsubmit, isub, jsub, ndo, ijob, ndone=0;
  gpointer retVal;
  gchar *URL;
  olong myTimeout, istream=0;
  gboolean *done;
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

  /* Make sure pool is using the correct function */
  if ((thread->pool) && (((GFunc)thread->pool->func)!=((GFunc)func))) {
     g_thread_pool_free(thread->pool, FALSE, TRUE);  /* Be patient */
     thread->pool = NULL;
     if (thread->queue) g_async_queue_unref (thread->queue);
     thread->queue = NULL;
  }
   
  /* Threading allowed - has the Thread Pool been started? */
  if (thread->pool==NULL) ObitThreadPoolInit (thread, nstreams, func, (gpointer**)args);

  /* Start asynchronous message queue */
  ObitThreadQueueInit (thread);

  /* Array of flags for jobs done */
  done = g_malloc(njobs*sizeof(gboolean));
  for (i=0; i<njobs; i++) done[i] = !doJob[i];

  /* Submit the first nstreams jobs -  RPC calls in threads */
  nsubmit = 0;
  isub    = 0;
  for (i=0; i<nstreams; i++) {
    if (doJob[i]) {
      args[i]->ijob    = i;
      args[i]->ithread = i;
      g_thread_pool_push (thread->pool, args[i], NULL);
      fprintf (stderr, "Start job %d on thread %d\n",i,i);
      done[i] = TRUE;
      isub = i;
      nsubmit++;
    }
  }
 
  /* Wait for them to finish, expects each to send a message to the asynchronous 
     queue iff they finish giving the URL actually used */
  /* set timeouts */
  myTimeout = timeout;   /* event queue timeout */
  if(myTimeout<=0.0) myTimeout = 10.0;
  while (ndone<njobs) {  /* Loop until all done */
    
    /* Impatiently wait for one to finish - run RPC event queue 
       then check async_queue */
    if (nsubmit>=ndo) myTimeout = 1000000000;            /* If Done - wait for all */
    add_time = (glong)(myTimeout * 1000000);
    retVal = ObitThreadQueueCheck(thread, add_time); /* Completion queue */
    
    /* if retVal is nonNULL, something finished */
    if (retVal) ndone++;
    else continue;
    ijob = *(olong*)retVal;
    istream = args[ijob]->ithread;
      
    /* More to submit? */
    if (nsubmit<ndo) {
      /* find next to submit */
      for (j=isub+1; j<njobs; j++) {
	if (doJob[j]) {
	  isub = j;
	  break;
	}
      }
      /* Recycle URL if any more work for it */
      URL = myClassInfo.URL[istream];
      /* Find next job for URL */
      jsub = -1;
      for (k=isub; k<ndo; k++) {
	/* If URL not give or matches use URL and dispatch this job */
	if (((args[k]->URL==NULL) || !strcmp(URL, args[k]->URL)) &&
	(!done[k])) {
	  jsub = k;
	  if (args[jsub]->URL==NULL) args[jsub]->URL = g_strdup(URL);
	  args[jsub]->client =  myClassInfo.clients[istream];
	  break;
	}
      }

      /* Stick next into the queue */
      if (jsub>=0) {
	args[jsub]->ijob    = jsub;
	args[jsub]->ithread = istream;
	g_thread_pool_push (thread->pool, args[jsub], NULL);
	fprintf (stderr, "Start job %d on thread %d\n",jsub,istream);
	done[jsub] = TRUE;
	nsubmit++;
      }
    } /* End finding next to submit */
    
    /* error check on jobs already finished */
    for (k=0; k<isub; k++) {
      /* Only want to check on the actual number submitted */
      if (doJob[k]) {
	if (args[k]->retVal==NULL) out = FALSE;
	if (args[k]->err->error)   {
	  out = FALSE;
	  Obit_traceback_val (args[k]->err, routine, "execution", out);  
	}
      } /* end if doJob */
    } /* end loop over jobs error check */
  } /* end loop till all jobs done*/

  /* Stop asynchronous message queue */
  ObitThreadQueueFree (thread);

  return out;
} /* end MultiProcExec */


