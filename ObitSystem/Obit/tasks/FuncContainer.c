/* $Id:  $  */
/* FuncContainer Obit function container process                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2010                                          */
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
#include "ObitImage.h"
#include "ObitXML.h"
#include "ObitRPC.h"
#include "ObitData.h"
#include "ObitFile.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* FuncContIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void FuncContOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Start Obit System */
ObitSystem* startupObit (ObitErr *err);
/* Shutdown Obit System */
void shutdownObit (ObitSystem *mySystem, ObitErr *err);
/* Shutdown Container */
gpointer shutdownContainer(gpointer arg);
/* Obit Message g_log handler */
void ObitMessageHandler (const gchar *log_domain,
			 GLogLevelFlags log_level,
			 const gchar *message,
			 gpointer logList);

/*---------------Public RPC callable function prototypes----------------*/

/* "Are you there?" function */
ObitXMLValue*
ping(ObitXMLEnv*   const envP, 
     ObitXMLValue* const paramArrayP,
     void*         const contData
);

/* test function */
ObitXMLValue*
test(ObitXMLEnv*    const envP, 
     ObitXMLValue*  const paramArrayP,
     void*          const contData
);
/* shutdown function */
ObitXMLValue*
terminate(ObitXMLEnv*    const envP, 
	  ObitXMLValue*  const paramArrayP,
	  void*          const contData
);

/* Program globals */
olong  port;              /* Port number to accept requests on */

gchar *pgmName = "FuncCont";       /* Program name */
gchar *infile  = "FuncCont.inp";   /* File with program inputs */
gchar *outfile = "FuncCont.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL;  /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myLog    = NULL; /* Message log list */

/*---------------Private structures----------------*/

/* Shutdown data structure */
typedef struct {
  /* die flag */
  gboolean doDie;
  /* Thread I'm running in */
  ObitThread *dieThread;
  /* ObitSystem */
  ObitSystem *mySystem;
  /* Obit Error/message structure  */
  ObitErr *err;
} ShutdownDataStruct;

/* Container data structure */
typedef struct {
  /* Obit Error/message structure */
  ObitErr *err;
  ShutdownDataStruct *myShutdownDataStruct;
} contDataStruct;

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit Function container - watch for RPC requests and execute         */
/*----------------------------------------------------------------------- */
{
  contDataStruct contData;
  ObitSystem   *mySystem= NULL;
  ShutdownDataStruct myShutdownDataStruct; 
  ObitErr *err;
  ObitRPC *server;
  guint mlog;
  oint  ierr = 0;
  
   /* Startup - parse command line */
  err = newObitErr();
  myInput = FuncContIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

 /* Start up Obit System */
  mySystem = startupObit (err);
  ObitErrLog(err); /* show any error messages on err */
  if (err->error) goto exit;

  /* Initialize deathThread data */
  myShutdownDataStruct.doDie     = FALSE;
  myShutdownDataStruct.dieThread = newObitThread();
  myShutdownDataStruct.err       = ObitErrRef(err);
  myShutdownDataStruct.mySystem  = mySystem;

  /* Initialize container data */
  contData.err = err;
  contData.myShutdownDataStruct = &myShutdownDataStruct;

  /* Start RPC server */
  server = ObitRPCCreateServer ("FuncContainer", err);

  /* Register functions */
  ObitRPCAddMethod (server, "ping",       &ping,      (void*)&contData, err);
  ObitRPCAddMethod (server, "test",       &test,      (void*)&contData, err);
  ObitRPCAddMethod (server, "terminate",  &terminate, (void*)&contData, err);
  if (err->error) goto exit;

  /* Start deathThread */
  ObitThreadStart1 (myShutdownDataStruct.dieThread, (ObitThreadFunc)shutdownContainer,
		    (gpointer)&myShutdownDataStruct);

  /* Initialize Threading */
  ObitThreadInit (myInput);

  /* DEBUG */
  Obit_log_error(err, OBIT_InfoErr, "%s: Started on port %d",pgmName,port);
  ObitErrLog(err); /* show any error messages on err */

  /* Set Obit Message g_log handler */
  myLog = newObitInfoList();
  mlog = g_log_set_handler (NULL,  G_LOG_LEVEL_WARNING | G_LOG_FLAG_FATAL |
			    G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_MESSAGE |
			    G_LOG_LEVEL_INFO | G_LOG_LEVEL_DEBUG |
			    G_LOG_FLAG_RECURSION, 
			    /*(GLogFunc)ObitMessageHandler, &myLog);*/
			    (GLogFunc)ObitMessageHandler, NULL);

  /* Loop forever 'neath the streets of Boston */
  ObitRPCServerLoop(server, port, "/tmp/xmlrpc_log");

  /* Shutdown  */
 exit:
  shutdownObit (mySystem, err);
  ObitErrLog(err); /* show any error messages on err */
  return ierr;
} /* end main */

ObitSystem* startupObit (ObitErr *err)
/*----------------------------------------------------------------------- */
/*   Startup Obit System                                                  */
/*----------------------------------------------------------------------- */
{
  ObitSystem   *mySystem= NULL;

   /* Startup */
  err = newObitErr();

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  return mySystem;
} /* end of startup */

void shutdownObit (ObitSystem *mySystem, ObitErr *err)
/*----------------------------------------------------------------------- */
/*   Shutdown Obit System                                                 */
/*----------------------------------------------------------------------- */
{
  ObitErrLog(err); /* show any error messages on err */

  /* Shutdown  */
  mySystem = ObitSystemShutdown (mySystem);
  myLog = ObitInfoListUnref(myLog);
 
  return;
} /* end of Shutdown */

ObitInfoList* FuncContIn (int argc, char **argv, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Parse control info from command line                                  */
/*   Input:                                                               */
/*      argc   Number of arguments from command line                      */
/*      argv   Array of strings from command line                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   return  ObitInfoList with defaults/parsed values                     */
/*----------------------------------------------------------------------- */
{
  olong ax;
  gchar *arg;
  gboolean init=FALSE;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint    itemp, i, j, k;
  ObitInfoList* list;
  ObitFile *parmFile=NULL;
  gchar *strTemp;
  gchar *routine = "FuncContIn";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* command line arguments */
  /* fprintf (stderr,"DEBUG arg %d %s\n",argc,argv[0]); DEBUG */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {

     /*fprintf (stderr,"DEBUG next arg %s %s\n",argv[ax],argv[ax+1]); DEBUG */
    arg = argv[ax];
    if (strcmp(arg, "-input") == 0){ /* input parameters */
      infile = argv[++ax];
      /* parse input file */
      ObitParserParse (infile, list, err);
      init = TRUE;

    } else if (strcmp(arg, "-output") == 0){ /* output results */
      outfile = argv[++ax];

    } else if (strcmp(arg, "-port") == 0) { /* Port number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListAlwaysPut (list, "port", OBIT_oint, dim, &itemp);
      /*fprintf (stdout,"port %d\n",itemp);*/
      
    } else if (strcmp(arg, "-pgmNumber") == 0) { /* Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-AIPSuser") == 0) { /* AIPS user number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "AIPSuser", OBIT_oint, dim, &itemp, err);
            
    } else { /* unknown argument */
      Usage();
    }
    if (err->error) Obit_traceback_val (err, routine, "GetInput", list);
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (infile, list, err);

  /* delete infile */
  parmFile = newObitFile(infile);
  ObitFileOpen (parmFile, infile, OBIT_IO_ReadWrite, OBIT_IO_Text, 0, err);
  ObitFileClose(parmFile, err);
  ObitFileZap(parmFile, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);
  
  /* Extract basic information to program globals */
  ObitInfoListGet(list, "port",      &type, dim, &port,      err);
  ObitInfoListGet(list, "pgmNumber", &type, dim, &pgmNumber, err);
  ObitInfoListGet(list, "AIPSuser",  &type, dim, &AIPSuser,  err);
  ObitInfoListGet(list, "nAIPS",     &type, dim, &nAIPS,     err);
  ObitInfoListGet(list, "nFITS",     &type, dim, &nFITS,     err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  /* Directories more complicated */
  ObitInfoListGetP(list, "AIPSdirs",  &type, dim, (gpointer)&strTemp);
  if (strTemp) {  /* Found? */
    AIPSdirs = g_malloc0(dim[1]*sizeof(gchar*));
    for (i=0; i<dim[1]; i++) {
      AIPSdirs[i] =  g_malloc0(dim[0]*sizeof(gchar));
      k = 0;
      for (j=0; j<dim[0]; j++) { /* Don't copy blanks */
	if (strTemp[j]!=' ') {AIPSdirs[i][k] = strTemp[j]; k++;}
      }
      AIPSdirs[i][k] = 0;
      strTemp += dim[0];
    }
  }

  ObitInfoListGetP(list, "FITSdirs",  &type, dim, (gpointer)&strTemp);
  if (strTemp)   {  /* Found? */
    FITSdirs = g_malloc0(dim[1]*sizeof(gchar*));
    for (i=0; i<dim[1]; i++) {
      FITSdirs[i] =  g_malloc0(dim[0]*sizeof(gchar));
      k = 0;
      for (j=0; j<dim[0]; j++) { /* Don't copy blanks */
	if (strTemp[j]!=' ') {FITSdirs[i][k] = strTemp[j]; k++;}
      }
      FITSdirs[i][k] = 0;
      strTemp += dim[0];
    }
  }

  return list;
} /* end FuncContIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: FuncCont -input file -output ofile [args]\n");
    fprintf(stderr, "Obit function container\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def FuncCont.in\n");
    fprintf(stderr, "  -output output result file, def FuncCont.out\n");
    fprintf(stderr, "  -port port number to watch \n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     pgmNumber Int        Program number (like POPS number) def 1       */
/*     nFITS     Int        Number of FITS directories [def. 1]           */
/*     FITSdirs  Str [?,?]  FITS directories [def {"./"}]                 */
/*     AIPSuser  Int        AIPS user number [def 2}]                     */
/*     nAIPS     Int        Number of AIPS directories [def. 1]           */
/*     AIPSdirs  Str [?,?]  AIPS directories [def {"AIPSdata/"}]          */
/*     DataType  Str [4]    "AIPS" or "FITS" [def {"FITS"}]               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint   itemp;
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListPut (out, "pgmNumber", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default port */
  itemp = 8764;
  ObitInfoListPut (out, "port", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default FITS directories - same directory */
  dim[0] = 1; dim[1] = 1;
  itemp = 0; /* number of FITS directories */
  ObitInfoListPut (out, "nFITS", OBIT_oint, dim, &itemp, err);

  /* AIPS user number */
  dim[0] = 1; dim[1] = 1;
  itemp = 2;
  ObitInfoListPut (out, "AIPSuser", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Default AIPS directories */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; /* number of AIPS directories */
  ObitInfoListPut (out, "nAIPS", OBIT_oint, dim, &itemp, err);

  return out;
} /* end defaultInputs */

/*------------------- Private functions -------------------*/

/**
 * Squirrels away ObitMessages in myLog with keys in form Mnnnnnnnn
 * Run ag g_log handler
 * \param log_domain
 * \param log_level
 * \param message message to be logged
 * \param logList
 */
void ObitMessageHandler (const gchar *log_domain,
			 GLogLevelFlags log_level,
			 const gchar *message,
			 gpointer logList)
{
  olong number;
  static gchar mlabel[40];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  if (myLog==NULL) return;

  number = myLog->number+1;
  /* Unique label for message */
  sprintf (mlabel,"M%8.8d",number);
  
  /* Add entry */
  dim[0] = strlen (message);
  ObitInfoListAlwaysPut (myLog, mlabel, OBIT_string, dim, (gchar*)message);
} /* end ObitMessageHandler */

/**
 * Lies in wait for doDie to turn positive, then shuts down the container
 * Runs in a thread.
 * \param arg  pointer to ShutdownDataStruct with:
 * \li doDie      Flag to indicate when to shutdown
 * \li dieThread  Thread I'm running in (for mutex)
 * \li mySystem   System to shutdown
 * \li err        Obit Error/message structure
 */
gpointer shutdownContainer (gpointer arg) 
{
  ShutdownDataStruct *largs = (ShutdownDataStruct*)arg;
  gboolean done = FALSE;

  /* Idle loop waiting for doDie */
  while (!done) {
      usleep(250000);  /* 250 msec */
      /* test */
      ObitThreadLock(largs->dieThread);
      if (largs->doDie) done = TRUE;
      ObitThreadUnlock(largs->dieThread);
  } /* end idle loop */

  /* Wait a decent interval */
  usleep(2000000);  /* 2 sec */
  fprintf (stderr, "%s on Port %d: I am shutting down\n", pgmName, port);

  /* Yeah - I'm free - kill everything in sight */
  shutdownObit (largs->mySystem, largs->err);
  exit(0);
}/* end shutdownContainer */

/*------------------- RPC callable functions -------------------*/
/**
 * Handle ping function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml, 
 *                      contains argument InfoList
 * \param contData      containter data
 * \return xml return value for function
 * This consists of an ObitInfoList containing:
 * \li Status  Completion status "OK", failed"
 * \li log     Message log as string
 * \li ping    42
 */
ObitXMLValue *
ping(ObitXMLEnv *   const envP, 
     ObitXMLValue * const paramArrayP,
     void *         const contData) 
{  
  contDataStruct *largs = (contDataStruct*)contData;
  ObitErr *err = largs->err;
  ObitInfoList *myInput = ObitXMLGetCallArg (envP, paramArrayP, err);
  ObitInfoList *outList=NULL;
  ObitXML *outXML=NULL;
  ObitXMLValue *outXMLValue=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *OK = "OK", *failed="Failed", *log = "I am here";
  olong ping=42;

   /* Initialize Return results */
  outList = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = strlen(failed);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, failed);
  dim[0] = strlen(log);  
  ObitInfoListAlwaysPut (outList, "log",    OBIT_string, dim, log);
  dim[0] = 1;  
  ObitInfoListAlwaysPut (outList, "ping",    OBIT_int, dim, &ping);
  
  /* error? */
  if (err->error) goto finish;

  /* Do operation */
  /* Nothing for ping */

  /* Must be OK */
  dim[0] = strlen(OK);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, OK);

  /* Done */
 finish:
  ObitErrLog(err); /* show any error messages on err */
  /* Create return value with status, message log and any return parameters */
  if (outXML) outXML = ObitXMLUnref(outXML);
  outXML = ObitXMLServerReturn (outList, err);
  ObitErrLog(err); /* show any error messages on err */
  outXMLValue = ObitXMLGetValue(outXML);

  /* Cleanup ;*/
  myInput = ObitInfoListUnref(myInput);
  outList = ObitInfoListUnref(outList);
  if (outXML) outXML = ObitXMLUnref(outXML);
  
  return outXMLValue;
} /* end ping */

/**
 * Handle terminate function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml, 
 *                      contains argument InfoList
 * \param contData      containter data
 * \return xml return value for function
 * This consists of an ObitInfoList containing:
 * \li Status  Completion status "OK", failed"
 * \li log     Message log as string
 */
ObitXMLValue *
terminate(ObitXMLEnv *   const envP, 
	  ObitXMLValue * const paramArrayP,
	  void *         const contData) 
{  
  contDataStruct *largs = (contDataStruct*)contData;
  ObitErr *err = largs->err;
  ObitInfoList *myInput = ObitXMLGetCallArg (envP, paramArrayP, err);
  ObitInfoList *outList=NULL;
  ObitXML *outXML=NULL;
  ObitXMLValue *outXMLValue=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *OK = "OK", *failed="Failed", *log = "terminating";

   /* Initialize Return results */
  outList = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = strlen(failed);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, failed);
  dim[0] = strlen(log);  
  ObitInfoListAlwaysPut (outList, "log",    OBIT_string, dim,log);

  /* error? */
  if (err->error) goto finish;

  /* Leave note for deathThread */
  ObitThreadLock(largs->myShutdownDataStruct->dieThread);
  largs->myShutdownDataStruct->doDie = TRUE;
  ObitThreadUnlock(largs->myShutdownDataStruct->dieThread);

  /* Must be OK */
  dim[0] = strlen(OK);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, OK);

  /* Done */
 finish:
  /* Create return value with status, message log and any return parameters */
  if (outXML) outXML = ObitXMLUnref(outXML);
  outXML = ObitXMLServerReturn (outList, err);
  outXMLValue = ObitXMLGetValue(outXML);

  /* Cleanup */
  myInput = ObitInfoListUnref(myInput);
  outList = ObitInfoListUnref(outList);
  if (outXML) outXML = ObitXMLUnref(outXML);

  /* DEBUGGING */
  outXML  = ObitXMLReturn ("Nonesuch", (ObitXMLValue*)outXMLValue, err);
  outList = ObitXMLGetServerResult(outXML, err);

  return outXMLValue;
} /* end terminate */

/**
 * Handle test function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml, 
 *                      contains argument InfoList
 * \param contData      containter data
 * \return xml return value for function
 * This consists of an ObitInfoList containing:
 * \li Status  Completion status "OK", "failed"
 * \li NumMessage  Number of log messages
 * \li Mnnnnnnnn   Log message number nnnnnnnn
 * \li any function dependent return values
 */
ObitXMLValue *
test(ObitXMLEnv *   const envP, 
     ObitXMLValue * const paramArrayP,
     void *         const contData) 
{  
  contDataStruct *largs = (contDataStruct*)contData;
  ObitErr *err = largs->err;
  ObitInfoList *myInput = ObitXMLGetCallArg (envP, paramArrayP, err);
  ObitInfoList *outList=NULL;
  ObitXML *outXML=NULL;
  ObitXMLValue *outXMLValue=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *OK = "OK", *failed="Failed";
  gboolean bfailed=TRUE;
  olong nlog;
  ObitImage *outImage=NULL;
  gchar *routine = "test";

  /* testing */
  ofloat *work=NULL;
  olong nwork, i, j;

  /* Initialize Return results */
  outList = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  dim[0] = strlen(failed);
  ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, failed);

  /* error? */
  if (err->error) goto finish;

  /* Do operation here */
  /* Create test image */
  outImage = ObitImageFromFileInfo("out", myInput, err);
  ObitImageOpen (outImage, OBIT_IO_ReadWrite, err);
  ObitImageClose (outImage, err);
  if (err->error) goto finish;

  /* Do something that takes a while */
  fprintf (stderr,"Start on port %d\n",port);
  nwork = 10000;
  work = g_malloc0(nwork*sizeof(ofloat));
  for (j=0; j<10000; j++) {
    for (i=0; i<nwork; i++) {
      work[i] = i*0.001;
      work[i] = cos(work[i]);
    }
  }

  /* DEBUG */
  ObitInfoListPrint (myInput, stdout);

  /* Test Logging */
  Obit_log_error(err, OBIT_InfoErr, "%s: Started on port %d",routine,port);
  Obit_log_error(err, OBIT_InfoErr, "%s: some logging message nwork %d",routine,nwork);
  Obit_log_error(err, OBIT_InfoErr, "%s: another logging message",routine);
  Obit_log_error(err, OBIT_InfoErr, "%s: yet another logging message",routine);
  /*Obit_log_error(err, OBIT_Error,   "%s: What the Hell",routine);*/
  Obit_log_error(err, OBIT_InfoErr, "%s: Finished",routine);

  /* Check if OK */
  if (!err->error) {
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut (outList, "Status", OBIT_string, dim, OK);
    bfailed = FALSE;
  }

  /* Done */
 finish:
  ObitErrLog(err); /* show any error messages on err */
  /* Move Logging messages to outList */
  nlog = myLog->number;
  if (bfailed) { /* DEBUG */
    ObitInfoListPrint (myLog, stdout);
  }    
  dim[0] = 1;
  ObitInfoListAlwaysPut (outList, "NumMessage",  OBIT_long, dim, &nlog);
  outList = ObitInfoListCopyData (myLog, outList);
  myLog = ObitInfoListUnref(myLog);   myLog = newObitInfoList();

  /* Create return value with status, message log and any return parameters */
  if (outXML) outXML = ObitXMLUnref(outXML);
  outXML = ObitXMLServerReturn (outList, err);
  outXMLValue = ObitXMLGetValue(outXML);
  
  /* Cleanup */
  myInput = ObitInfoListUnref(myInput);
  outList = ObitInfoListUnref(outList);
  if (outXML) outXML = ObitXMLUnref(outXML);
 
  return outXMLValue;
} /* end test */

