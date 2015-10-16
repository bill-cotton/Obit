/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014-2015                                          */
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

#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitUV.h"
#include "ObitData.h"
#include "ObitTableUtil.h"
#include "ObitTableSY.h"
#include "ObitUVSoln.h"
#include "ObitGainCal.h"
#include "ObitThread.h"

/*---------------Private structures----------------*/
/* Median filtering threaded function argument */
typedef struct {
  /* ObitThread with restartable queue */
  ObitThread *thread;
  /* thread number, <0 -> no threading   */
  olong        ithread;
  /* total number of values in arrays  */
  olong        ntime;
  /*  first (0-rel) value to smooth */
  olong lo;
  /*  highest (0-rel) value to smooth */
  olong hi;
  /* Controls smoothing 
     0 -> 1 = pure boxcar -> pure MWF (alpha of the 
     center data samples are averaged).
  */
  ofloat alpha;   
  /* smoothing time in unints of times */  
  ofloat smoTime;
  /* Clipping level in sigma, <=0 = no clip */  
  ofloat clip;
  /*  Array of values to smooth */
  ofloat *vals;
  /*  Smoothed Array of values */
  ofloat *out;
  /*  Array of weights (0 or 1) */
  ofloat *weight;
  /*  Array of output weights (0 or 1) */
  ofloat *outWt;
  /* times of data samples, <-1000 => ignore */
  ofloat* times;   
  /* Work array at least the size of hi-lo */
  ofloat* work;
} MednFuncArg;

/* internal prototypes */
/* Get inputs */
ObitInfoList* SYGainIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SYGainOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Write history */
void SYGainHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Copy selected data */
void SYCopy (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err);
/* Do Any Clipping */
void SYGainClip (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Smooth */
void SYGainSmooth (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Convert SY to SN table */
void SYGainSY2SN (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Clipping routines */
void clipPSum (ObitTableSY *SYTable, ObitUVSel *sel,  olong iif,
	       ofloat alpha, ofloat stpsum, ofloat mxpsum, olong sub,	
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       olong nThread, MednFuncArg** args, ObitErr* err);
void clipPDif (ObitTableSY *SYTable, ObitUVSel *sel, olong iif, 
	       ofloat alpha, ofloat stpdif, ofloat mxpdif, olong sub, 
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       olong nThread, MednFuncArg** args, ObitErr*  err);
void clipGain (ObitTableSY *SYTable, ObitUVSel *sel, olong iif, 
	       ofloat alpha, ofloat strate, ofloat mxrate, olong sub, 
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       olong nThread, MednFuncArg** args, ObitErr* err);
/*  Smooth SY table */
void static
SYSmooth (ObitTableSY *SYTab, gchar* smoFunc, gchar* smoType, ofloat alpha, 
	  ofloat *smoParm, olong iif, olong sub, 
	  olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, 
	  olong nThread, MednFuncArg** args, ObitErr* err);

/* Generic smoothing */
void static 
SYsmoIt (gchar* smmeth, ofloat width, ofloat alpha, ofloat clip, 
	 ofloat* x, ofloat *t, ofloat *w, olong n, 
	 ofloat* xs, ofloat* ws, ofloat *wrk1, ofloat *wrk2, 
	 gboolean doBlank, olong nThread, MednFuncArg** args);
/* count times in SY table */
olong SYCountTime (ObitTableSY *SYTab, olong isub, ObitErr* err);
/* Thread argument routines */
static olong MakeSmoFuncArgs (ObitThread *thread, olong lenwrk, 
			      MednFuncArg ***ThreadArgs);
static void KillSmoFuncArgs (olong nargs, MednFuncArg **ThreadArgs);
static gpointer ThreadBoxSmo (gpointer arg);
static gpointer ThreadMWFSmo (gpointer arg);
static gpointer ThreadGausSmo (gpointer arg);
/** Private: qsort ofloat comparison */
static int compare_ofloat  (const void* arg1,  const void* arg2);

/* Program globals */
gchar *pgmName = "SYGain";       /* Program name */
gchar *infile  = "SYGain.in" ;   /* File with program inputs */
gchar *outfile = "SYGain.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL ; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
olong inSYVer, outSYVer, outSYVer;   /* Input and output SY/SY tables */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task which smooths and filters Solution(SY) tables.             */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = SYGainIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Copy selected table */
  SYCopy (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Clipping */
  SYGainClip (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* smooth solutions */
  SYGainSmooth (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Convert SY to SN table */
  SYGainSY2SN (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  SYGainHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SYGainIn (int argc, char **argv, ObitErr *err)
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
  gchar *strTemp;
  oint    itemp, i, j, k;
  ObitInfoList* list=NULL;
  gchar *routine = "SYGainIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

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

    } else if (strcmp(arg, "-pgmNumber") == 0) { /*Program number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "pgmNumber", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-AIPSuser") == 0) { /* AIPS user number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "AIPSuser", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inSeq") == 0) { /* AIPS sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-inDisk") == 0) { /* input disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "inDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inName") == 0) { /* AIPS inName*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inClass") == 0) { /* AIPS inClass*/
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-inFile") == 0) { /*inFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inFile", OBIT_string, dim, strTemp);
      
   } else { /* unknown argument */
      Usage();
    }
    if (err->error) Obit_traceback_val (err, routine, "GetInput", list);
  } /* end parsing input arguments */
  
  /* Read defaults if no file specified */
  if (!init) ObitParserParse (infile, list, err);

  /* Extract basic information to program globals */
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

  /* Initialize output */
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end SYGainIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SYGain -input file -output ofile [args]\n");
    fprintf(stderr, "Clip and smooth an SY to SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SYGain.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def SYGain.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2\n");
    fprintf(stderr, "  -inFile input FITS uvdata file\n");
    fprintf(stderr, "  -inName input AIPS uvdata file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
/*  Note: Other parameters may be passed through the input text file      */
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
/*     inFile    Str [?]    input FITS image file name [def "Image.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListPut (out, "pgmNumber", OBIT_oint, dim, &itemp, err);
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

  /* Default type "FITS" */
  strTemp = "FITS";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "DataType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "SYGain.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SYGainName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Sources selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Sources", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Subarray */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "subA", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  ObitInfoList *out = newObitInfoList();
  /*  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
      ofloat ftemp;
      gchar *routine = "defaultOutputs";"*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;

  /* add parser items - nothing */
  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Digest inputs                                                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  oint         doSYGain;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doSYGain = -1;
  ObitInfoListGetTest(myInput, "doSYGain",  &type, dim, &doSYGain);
  doCalSelect = doCalSelect || (doSYGain>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  /* Initialize Threading */
  ObitThreadInit (myInput);

} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data, builds selector                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitUV with input data                                           */
/*----------------------------------------------------------------------- */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *inData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  olong        doCalib, nvis=1000;
  gchar        *dataParms[] = {  /* Parameters to select data */
    "Sources", "timeRange",  "subA", "Antennas", "FredID", "souCode", "Qual", 
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);
    
  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 
  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated, with selection and OK */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for SYGain                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SYGainHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "FreqID", "souCode", "Qual", "timeRange",  "subA", "Antennas", 
    "SWUse", "SYVer", "SYOut", "solnOut", "smoFunc", "smoParm", "alpha",
    "calInt", "refAnt", "clipSmo", "clipParm",
    NULL};
  gchar *routine = "SYGainHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

   /* Do history */
  outHistory = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadWrite, err);

  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  outHistory = ObitHistoryUnref(outHistory); /* cleanup */
 
} /* end SYGainHistory  */

/**
 * Routine to copy the input table to the output.
 * Can substitute values from one SW into another, uses "SWUse" on 
 * myInput
 *  \param     myInput   Input parameters on InfoList 
 *  \param     inData    ObitUV with tables  
 *  \param     err       Obit Error stack                
 */
void SYCopy (ObitInfoList* myInput, ObitUV* inData,  ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableSY *inTab=NULL, *outTab=NULL;
  ObitTableSYRow *row=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong   itemp, isuba,  fqid, npol, nif, iif;
  olong SYIn, SYOut, highVer;
  gboolean dropIt;
  ofloat tr[2];
  olong  loop, outrow, isub, *SWUse=NULL;
  gchar  *routine = "SYGainCopy";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  g_assert(ObitUVIsA(inData));
 
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SY");
  itemp = 0;
  ObitInfoListGetTest(myInput, "SYVer", &type, dim, &itemp);
  SYIn = itemp;
  if (SYIn==0) SYIn = highVer;
  /* Save to inputs */
  dim[0] = dim[1] = dim[2] = 1; itemp = SYIn;
  ObitInfoListAlwaysPut(myInput, "SYVer", OBIT_long, dim, &itemp);
  inSYVer = itemp;

  ObitInfoListGetTest(myInput, "SYOut", &type, dim, &itemp);
  SYOut = itemp;
  if (SYOut==0) SYOut = highVer+1;
  /* Save to inputs */
  dim[0] = dim[1] = dim[2] = 1; itemp = SYOut;
  ObitInfoListAlwaysPut(myInput, "SYOut", OBIT_long, dim, &itemp);
  outSYVer = itemp;

  inTab = newObitTableSYValue (inData->name, (ObitData*)inData, &SYIn, 
			       OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  /* Open input table */
  retCode = ObitTableSYOpen (inTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);

  outTab = newObitTableSYValue (inData->name, (ObitData*)inData, &SYOut, 
				OBIT_IO_ReadWrite, inTab->nIF, inTab->nPol, 
				err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Clear any existing rows */
  ObitTableClearRows((ObitTable*)outTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open output table */
  retCode = ObitTableSYOpen (outTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  outTab->nAnt = inTab->nAnt;  /* Save number of antennas */
  ((ObitTableDesc*)outTab->myIO->myDesc)->sort[0] = 
    ((ObitTableDesc*)inTab->myIO->myDesc)->sort[0];
  ((ObitTableDesc*)outTab->myIO->myDesc)->sort[1] = 
    ((ObitTableDesc*)inTab->myIO->myDesc)->sort[1];
  
  /* Get selection criteria */
  isuba  = inData->mySel->SubA;
  fqid   = inData->mySel->FreqID;
  nif    = inTab->nIF;
  npol   = inTab->nPol;

  /* Get substitution table */
  ObitInfoListGetP(myInput, "SWUse",  &type, dim, (gpointer)&SWUse);

  /* Timerange */
  tr[0] = inData->mySel->timeRange[0];
  tr[1] = inData->mySel->timeRange[1];
  if ((tr[0]==0.0) && (tr[1]==0.0)) {
    tr[0] = -1.0e20;
    tr[1] =  1.0e20;
  }

  /* Create Row */
  row = newObitTableSYRow (outTab);
  /* Attach row to output buffer */
  ObitTableSYSetRow (outTab, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);

  /* Tell */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Copy SY %d to SY %d", inSYVer, SYOut);
  ObitErrLog(err); 

  /* Loop through table */
  for (loop=1; loop<=inTab->myDesc->nrow; loop++) { /* loop 20 */

    retCode = ObitTableSYReadRow (inTab, loop, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outTab->name);
    if (row->status<0) continue;  /* Skip deselected records */

    /*  Drop this one? */
    dropIt = (row->SubA!=isuba) && (isuba>0);                    /* by subarray */
    dropIt = dropIt || ((row->FreqID!=fqid) && (fqid>0));        /* by FQ id */
    dropIt = dropIt || ((row->Time<tr[0]) || (row->Time>tr[1])); /* by time */
    if (dropIt) continue; /* Drop */

    if (!ObitUVSelWantAnt (inData->mySel,  row->antennaNo))  continue;  /* Check Antenna */
    if (!ObitUVSelWantSour (inData->mySel, row->SourID))     continue;  /* Check Source */

    /* Substitute if necessary */
    if (SWUse) {
      for (iif=0; iif<nif; iif++) {
	if (SWUse[iif]>0) {
	  isub = SWUse[iif]-1;  /* 0 rel value */
	  row->PwrDif1[iif] = row->PwrDif1[isub];
	  row->PwrSum1[iif] = row->PwrSum1[isub];
	  row->Gain1[iif]   = row->Gain1[isub];
	  /* two poln? */
	  if (npol>1) {
	    row->PwrDif2[iif] = row->PwrDif2[isub];
	    row->PwrSum2[iif] = row->PwrSum2[isub];
	    row->Gain2[iif]   = row->Gain2[isub];
	  }
	} /* end substitute */
      } /* end if loop */
    }  /* End substitute */

    /* Write record to output */
    outrow = -1;
    retCode = ObitTableSYWriteRow (outTab, outrow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  } /* End loop over table */

  /* Close tables */
  retCode = ObitTableSYClose (inTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inTab->name);
  retCode = ObitTableSYClose (outTab, err);
  if (err->error) Obit_traceback_msg (err, routine, outTab->name);
  row = ObitTableSYRowUnref (row);  /* Cleanup */
} /* end SYCopy */

void SYGainClip (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*  Do any Clipping comparing with smoothed data                          */
/*  Data selection based on selector on inData                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with tables                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitTableSY *SYTable =  NULL;
  ObitThread *thread=NULL;
  MednFuncArg **args=NULL;
  olong SYver, highVer, nThreads=1;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat clipSmo[5]  = {0.0, 0.0, 0.0, 0.0, 0.0};
  ofloat alpha=0.5, stgain, stpdif, stpsum;
  ofloat clipParm[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  ofloat mxgain, mxpdif, mxpsum;
  gboolean dogain, dopdif, dopsum;
  olong itemp, numSub, iSub, subA, bif=1, eif=0, i;
  olong mxtim, *wrkrec=NULL;
  ofloat *wrktim=NULL, *work1=NULL, *work2=NULL, *work3=NULL, *work4=NULL,
    *work5=NULL, *work6=NULL, *work7=NULL, *work8=NULL;
  gchar *routine = "SYGainClip";

  /* Get SY table to clip */
  itemp = 0;
  ObitInfoListGetTest(myInput, "SYOut", &type, dim, &itemp);
  SYver = itemp;
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SY");
  if (SYver==0) SYver = highVer;
  SYTable = newObitTableSYValue (inData->name, (ObitData*)inData, &SYver, 
				 OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get clipping parameters */
  ObitInfoListGetTest(myInput, "alpha",    &type, dim, &alpha);
  ObitInfoListGetTest(myInput, "clipSmo",  &type, dim, clipSmo);
  stpdif  = clipSmo[0] / 24.0;
  stpsum  = clipSmo[1] / 24.0;
  stgain  = clipSmo[2] / 24.0;
  ObitInfoListGetTest(myInput, "clipParm", &type, dim, clipParm);
  mxpdif  = fabs (clipParm[0]);
  mxpsum  = fabs (clipParm[1]);
  mxgain  = fabs (clipParm[2]);
  ObitInfoListGetTest(myInput, "BIF",      &type, dim, &bif);
  if (bif==0) bif = 1;
  if (inData->myDesc->jlocif>=0) 
    bif = MAX (1, MIN (bif, inData->myDesc->inaxes[inData->myDesc->jlocif]));
  else bif = 1;
  ObitInfoListGetTest(myInput, "EIF",      &type, dim, &eif);
  if (inData->myDesc->jlocif>=0) {
    if (eif<=0) eif = inData->myDesc->inaxes[inData->myDesc->jlocif];
    eif = MAX (1, MIN (eif, inData->myDesc->inaxes[inData->myDesc->jlocif]));
  } else eif = 1;
  ObitInfoListGetTest(myInput, "subA",     &type, dim, &subA);

  /* What's to be clipped? */
  dopdif = mxpdif >=  1.0e-10;
  dopsum = mxpsum >=  1.0e-10;
  dogain = mxgain >=  1.0e-18;

  /* Anything to do? */
  if (!dopdif && !dopsum && !dogain)
    {ObitTableSYUnref(SYTable); return;}

  /* Zero smoothing time actually very large */
  if (mxpdif < 1.0e-10) mxpdif = 1.0e20;
  if (mxpsum < 1.0e-10) mxpsum = 1.0e20;
  if (mxgain < 1.0e-18) mxgain = 1.0e20;

  /* Number of subarrays = number AN tables */
  numSub = ObitTableListGetHigh (inData->tableList, "AIPS AN");

  /* Sort to antenna time order */
  ObitTableUtilSort2f ((ObitTable*)SYTable, "ANTENNA", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Open table */
  ObitTableSYOpen (SYTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create work arrays - make big enough to swallow whole table */
  mxtim = SYTable->myDesc->nrow;
  wrkrec = g_malloc(mxtim*sizeof(olong));
  wrktim = g_malloc(mxtim*sizeof(ofloat));
  work1  = g_malloc(mxtim*sizeof(ofloat));
  work2  = g_malloc(mxtim*sizeof(ofloat));
  work3  = g_malloc(mxtim*sizeof(ofloat));
  work4  = g_malloc(mxtim*sizeof(ofloat));
  work5  = g_malloc(mxtim*sizeof(ofloat));
  work6  = g_malloc(mxtim*sizeof(ofloat));
  work7  = g_malloc(mxtim*sizeof(ofloat));
  work8  = g_malloc(mxtim*sizeof(ofloat));

  /* Tell */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Clipping SY %d ",  SYver);
  ObitErrLog(err); 

  /* Create threading arguments  */
  thread = newObitThread();
  nThreads = MakeSmoFuncArgs (thread, mxtim, &args);

  /* Loop over subarrays */
  for (iSub=1; iSub<=numSub; iSub++) {
    /* This one wanted? */
    if ((subA > 0)  &&  (iSub != subA)) continue;
    
    /* Clip Power Dif */
    if (dopdif) {
      Obit_log_error(err, OBIT_InfoErr, "Clipping PDif ");
      ObitErrLog(err); 
      /* Loop over IF */
      for (i= bif; i<=eif; i++) { /* loop 100 */
	/*fprintf (stderr,"DEBUG IF=%d\n", i);*/
	clipPDif (SYTable, inData->mySel, i, alpha, stpdif, mxpdif, iSub, 
		 mxtim, wrkrec, wrktim, work1, work2, work3,
		 work4, work5, work6, work7, work8,
		 nThreads, args, err);
	if (err->error) goto cleanup;
      } /* end loop  L100: */;
    }

    /* Clip Power Sum */
    if (dopsum) {
      Obit_log_error(err, OBIT_InfoErr, "Clipping PSum ");
      ObitErrLog(err); 
      /* Loop over IF */
      for (i= bif; i<=eif; i++) { /* loop 100 */
	/*fprintf (stderr,"DEBUG IF=%d\n", i);*/
	clipPSum (SYTable, inData->mySel, i, alpha, stpsum, mxpsum, iSub, 
		  mxtim, wrkrec, wrktim, work1, work2, work3,
		  work4, work5, work6, work7, work8,
		  nThreads, args, err);
	if (err->error) goto cleanup;
      } /* end loop */
    }
    /* Clip gain */
    if (dogain) {
      Obit_log_error(err, OBIT_InfoErr, "Clipping Gain");
      ObitErrLog(err); 
      /* Loop over IF */
      for (i= bif; i<=eif; i++) { /* loop 100 */
	clipGain(SYTable, inData->mySel, i, alpha, stgain, mxgain, iSub, 
		 mxtim, wrkrec, wrktim, work1, work2, work3, 
		 work4, work5, work6, work7, work8,
		 nThreads, args, err);
	if (err->error) goto cleanup;
      } 
   } /* end loop  L100: */;
  } /* end loop over subarrays */
    
 cleanup:
  /* Close SY Table */
  ObitTableSYClose (SYTable, err);
  SYTable = ObitTableSYUnref(SYTable);
  KillSmoFuncArgs (nThreads, args);
  thread = ObitThreadUnref(thread);
  if (wrkrec) g_free (wrkrec);
  if (wrktim) g_free (wrktim);
  if (work1)  g_free (work1);
  if (work2)  g_free (work2);
  if (work3)  g_free (work3);
  if (work4)  g_free (work4);
  if (work5)  g_free (work5);
  if (work6)  g_free (work6);
  if (work7)  g_free (work7);
  if (work8)  g_free (work8);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
} /* end SYGainClip */

void SYGainSY2SN (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*  Convert SY to SN table                                                */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with tables                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitTable *outCal=NULL;
  gboolean Tr=TRUE, Fl=FALSE;
  olong solnOut=1, SYOut=2, highVer;
  ofloat calInt=0.1;
  gchar *routine = "SYGainSY2SN";

  if (err->error) return;  /* Previous error? */

  /* Input control */
  ObitInfoListGetTest(myInput, "SYOut",   &type, dim, &SYOut);
  ObitInfoListGetTest(myInput, "solnOut", &type, dim, &solnOut);
  ObitInfoListGetTest(myInput, "calInt",  &type, dim, &calInt);
  if (calInt<=0.0) calInt = 0.1;  /* 6 sec default */
  /* Convert min to seconds */
  calInt *= 60.0;

  /* Set control parameters on inData->info */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (inData->info, "doGain",  OBIT_bool,  dim, &Fl);
  ObitInfoListAlwaysPut (inData->info, "doSens",  OBIT_bool,  dim, &Fl);
  ObitInfoListAlwaysPut (inData->info, "doOpac",  OBIT_bool,  dim, &Fl);
  ObitInfoListAlwaysPut (inData->info, "doSmoo",  OBIT_bool,  dim, &Fl);
  ObitInfoListAlwaysPut (inData->info, "doSwPwr", OBIT_bool,  dim, &Tr);
  ObitInfoListAlwaysPut (inData->info, "calInt",  OBIT_float, dim, &calInt);
  ObitInfoListAlwaysPut (inData->info, "calVer",  OBIT_int,   dim, &solnOut);
  ObitInfoListAlwaysPut (inData->info, "SYVer",   OBIT_int,   dim, &SYOut);

  /* Delete old SN Table */
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  if ((solnOut>0) && (solnOut<=highVer))
    ObitUVZapTable (inData, "AIPS SN", solnOut, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Tell */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Convert SY %d to SN %d", SYOut, solnOut);
  ObitErrLog(err); 

  /* Convert */
  outCal = ObitGainCalCalc (inData, Tr, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Cleanup */
  outCal = ObitTableUnref(outCal);

} /* end SYGainSY2SN */

void SYGainSmooth (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
/*----------------------------------------------------------------------- */
/*   Smooth                                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV with tables                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
{
  ObitTableSY *SYTable =  NULL;
  olong SYver, highVer, nwork, nThreads=1;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong itemp, numSub, iSub, subA, numif, numpol, numant, iif, mxtime, ntmp;
  gchar smoFunc[5], smoType[5]="DIF ";
  ofloat alpha=0.5, smoParm[5], *work1=NULL, *work2=NULL;
  ObitThread *thread=NULL;
  gboolean doBlank=TRUE;
  MednFuncArg **args=NULL;
  gchar *routine = "SYGainSmooth";

  /* Get control parameters */
  ObitInfoListGetTest(myInput, "alpha",    &type, dim, &alpha);
  ObitInfoListGetTest(myInput, "smoFunc",  &type, dim, smoFunc);
  ObitInfoListGetTest(myInput, "smoType",  &type, dim, smoType);
  ObitInfoListGetTest(myInput, "smoParm",  &type, dim, smoParm);

  /* Convert times hr to days */
  smoParm[0] /= 24.0;
  smoParm[1] /= 24.0;
  smoParm[2] /= 24.0;

  /* Get SY table to rereference */
  itemp = 0;
  ObitInfoListGetTest(myInput, "SYOut", &type, dim, &itemp);
  SYver = itemp;
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SY");
  if (SYver==0) SYver = highVer;
  outSYVer = SYver;
  SYTable  = newObitTableSYValue (inData->name, (ObitData*)inData, &SYver, 
				  OBIT_IO_ReadWrite, 0, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get descriptive info */
  numif  = SYTable->nIF;
  numpol = SYTable->nPol;
  numant = SYTable->nAnt;

   /* Sort to time-antenna  order */
  ObitTableUtilSort2f ((ObitTable*)SYTable,"TIME  " , 1, FALSE, "ANTENNA", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
 
  /* Count number of times appearing in table */
  iSub = 0;
  mxtime = SYCountTime(SYTable, iSub, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
  mxtime += 100;  /* Fudge a bit on the number of times */

  /* Primary work array */
  work1 = g_malloc(10*mxtime*sizeof(ofloat));
  
  /* How many secondary work arrays */
  if (!strncmp(smoFunc, "GAUS", 4)) ntmp = 3;
  else if (!strncmp(smoFunc, "MWF", 3)) ntmp = 4;
  else ntmp = 2;
  work2 = g_malloc(ntmp*mxtime*sizeof(ofloat));

  /* Sort to antenna-time order */
  ObitTableUtilSort2f ((ObitTable*)SYTable, "ANTENNA NO.", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) goto cleanup;
 
  /* Open table */
  ObitTableSYOpen (SYTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* desired subarray */
  subA = 0;
  numSub = (olong)ObitTableListGetHigh (inData->tableList, "AIPS AN");
  ObitInfoListGetTest(myInput, "subA",  &type, dim, &subA);

  /* Tell */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Smoothing SY %d ",  SYver);
  ObitErrLog(err); 

  /* Create threading arguments if MWF */
  if (!strncmp (smoFunc, "MWF",3)) nwork = mxtime;
  else                             nwork = 0;
  thread = newObitThread();
  nThreads = MakeSmoFuncArgs (thread, nwork, &args);


  /* Loop over subarrays */
  for (iSub=1; iSub<=numSub; iSub++) {
    /* This one wanted? */
    if ((subA > 0)  &&  (iSub != subA)) continue;
    
    for (iif=0; iif<numif; iif++) { /* loop over IFs*/
      /* Smooth */
      SYSmooth (SYTable, smoFunc, smoType, alpha, smoParm, iif, iSub, 
		mxtime, work1, work2, doBlank, nThreads, args, err);
      if (err->error) goto cleanup;
    } /* end IF loop */
  } /* end loop over subarrays */

  /* Sort to time-antenna  order */
  ObitTableUtilSort2f ((ObitTable*)SYTable,"TIME  " , 1, FALSE, "ANTENNA", 
		       1, FALSE, err);
  if (err->error) goto cleanup;

  /* Close output table */
cleanup: 
  ObitTableSYClose (SYTable, err);
  /* Cleanup */
  if (work1) g_free(work1);
  if (work2) g_free(work2);
  if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
  SYTable = ObitTableSYUnref(SYTable);
  KillSmoFuncArgs (nThreads, args);
  thread = ObitThreadUnref(thread);
} /* end SYGainSmooth */

/*----------------------------------------------------------------------- */
/*  Clip Power difference                                                 */
/*   Input:                                                               */
/*      SYTable  Table to clip                                            */
/*      iif      IF to operate on, 1-rel                                  */
/*      sel      Selector telling which data wanted                       */
/*      stpdif   Smoothing time (day) for amplitude                       */
/*      alpha    0 -> 1 = pure boxcar -> pure MWF (ALPHA of the data      */
/*               samples are discarded and the rest averaged).            */
/*      mxpdif   Max. amplitude excursion from smoothed value             */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*      nThreads number of threads to use for MWF                         */
/*      args     thread argument array                                    */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
 void clipPDif (ObitTableSY *SYTable, ObitUVSel *sel, olong iif, 
		ofloat alpha, ofloat stpdif, ofloat mxpdif, olong sub, 
		olong maxtim, olong *wrkrec, ofloat *wrktim,
		ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
		ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
		olong nThreads, MednFuncArg** args, ObitErr* err)   {
   ObitTableSYRow *row=NULL;
   ofloat fblank =  ObitMagicF();
   olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
   olong  isyrno=0;
   gboolean bad, bad2, want;
   odouble timoff=0.0;
   gchar *routine = "clipPDif";
   
   /* Create Row */
   row = newObitTableSYRow (SYTable);
   /* Attach row to output buffer */
   ObitTableSYSetRow (SYTable, row, err);
   fstrec = 0;  /* Record number read in table */
   
   /* Loop over antenna */
   for (loopa= 1; loopa<=SYTable->nAnt; loopa++) {
     ant = loopa;
     /* Want this antenna? */
     if (!ObitUVSelWantAnt (sel, ant)) continue;
     
     /* Set pointers, counters */
     numtim = 0;
     nleft = SYTable->myDesc->nrow - fstrec;  /* How many rows? */
     
     /* Loop in time, reading */
     for (loopr= 1; loopr<=nleft; loopr++) {
       isyrno = fstrec + loopr;
       ObitTableSYReadRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
       if (row->status<0) continue;  /* Skip deselected records */
       
       /* Finished antenna? */
       if (row->antennaNo < ant) continue;
       if (row->antennaNo > ant) break;
       
       /* See if wanted. */
       want = ObitUVSelWantSour (sel, row->SourID);
       
       /* Check subarray */
       want = want  &&  (row->SubA == sub);
       
       /* Not all antennas wanted */
       want = want  &&  (row->antennaNo == ant);
       if (want) {
	 /* See if flagged value */
	 bad = (row->PwrDif1[iif-1] == fblank);
	 bad2 = (SYTable->nPol <= 1)  ||  (row->PwrDif2[iif-1] == fblank);
	 if (numtim < maxtim) {
	   if (numtim == 0) timoff = row->Time;
	   wrktim[numtim] = row->Time - timoff;
	   wrkrec[numtim] = isyrno;
	   if (bad) {
	     work2[numtim] = fblank;
	     work4[numtim] = 0.0;
	   } else {
	     work2[numtim] = row->PwrDif1[iif-1];
	     work4[numtim] = 1.0;
	   } 
	   if (bad2) {
	     work3[numtim] = fblank;
	     work5[numtim] = 0.0;
	   } else {
	     work3[numtim] = row->PwrDif2[iif-1];
	     work5[numtim] = 1.0;
	   } 
	   numtim++;
	 }
       } 
     } /* end reading loop */;
     save = isyrno - 1; /* How far did we get? */
     
     if (numtim <= 0) goto endAnt;  /* Catch anything? */
     
     /* Smooth/clip as requested */
     /*ObitUVSolnSmooMWF (stpdif, 1.0, wrktim, work2, work4, numtim, work1, 
       work6, work7, work8, FALSE);*/
     SYsmoIt ("MWF", stpdif, alpha, mxpdif, wrktim, work2, work4, numtim, 
	      work1, work6, work7, work8, FALSE, nThreads, args);
     /* Second Poln? */
     if (SYTable->nPol > 1) {
       /*ObitUVSolnSmooMWF (stpdif, 1.0, wrktim, work3, work5, numtim, work2, 
	 work4, work7, work8, FALSE);*/
       SYsmoIt ("MWF", stpdif, alpha, mxpdif, wrktim, work3, work5, numtim, 
	      work2, work4, work7, work8, FALSE, nThreads, args);
     } 
     
     /* Apply Clip */
     for (itime=0; itime<numtim; itime++) { /* loop 200 */
       /* either poln blanked? */
       if (work6[itime]<=0.0) continue;
       if ((SYTable->nPol > 1) &&  (work4[itime]<=0.0)) continue;
       isyrno = wrkrec[itime];
       ObitTableSYReadRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
       
       /* Blank */
       if (work6[itime]<=0.0) row->PwrDif1[iif-1] = fblank;
       if (work4[itime]<=0.0) row->PwrDif2[iif-1] = fblank;

       /* Rewrite record */
       ObitTableSYWriteRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
     } /* end loop clipping  */;
   endAnt: fstrec = save;
   } /* end loop  of antenna loop */;
   
   /* Cleanup */
   row = ObitTableSYRowUnref (row);  
   
 } /* end clipPDif */
 
/*----------------------------------------------------------------------- */
/*  Clip Power Sum                                                        */
/*   Input:                                                               */
/*      SYTable  Table to clip                                            */
/*      iif      IF to operate on, 1-rel                                  */
/*      sel      Selector telling which data wanted                       */
/*      stpdif   Smoothing time (day) for amplitude                       */
/*      alpha    0 -> 1 = pure boxcar -> pure MWF (ALPHA of the data      */
/*               samples are discarded and the rest averaged).            */
/*      mxpdif   Max. amplitude excursion from smoothed value             */
/*      dopdif   Clip by Amplitude?                                       */
/*      stpsum   Smoothing time (day) for phases                          */
/*      mxpsum   Max. excursion from smoothed value                       */
/*      doph     Clip by Phase?                                           */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*      nThreads number of threads to use for MWF                         */
/*      args     thread argument array                                    */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
  void clipPSum (ObitTableSY *SYTable, ObitUVSel *sel, olong iif, 
		 ofloat alpha, ofloat stpsum, ofloat mxpsum, olong sub, 
		 olong maxtim, olong *wrkrec, ofloat *wrktim,
		 ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
		 ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
		 olong nThreads, MednFuncArg** args, ObitErr* err)   {
    ObitTableSYRow *row=NULL;
    ofloat fblank =  ObitMagicF();
    olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
    olong  isyrno=0;
    gboolean   bad, bad2, want;
    odouble timoff=0.0;
    gchar *routine = "clipPSum";
    
    /* Create Row */
    row = newObitTableSYRow (SYTable);
    /* Attach row to output buffer */
    ObitTableSYSetRow (SYTable, row, err);
    fstrec = 0;  /* Record number read in table */
    
    /* Loop over antenna */
    for (loopa= 1; loopa<=SYTable->nAnt; loopa++) { 
      ant = loopa;
      /* Want this antenna? */
      if (!ObitUVSelWantAnt (sel, ant)) continue;
      
      /* Set pointers, counters */
      numtim = 0;
      nleft = SYTable->myDesc->nrow - fstrec;  /* How many rows? */
      
      /* Loop in time, reading */
      for (loopr= 1; loopr<=nleft; loopr++) { 
	isyrno = fstrec + loopr;
	ObitTableSYReadRow (SYTable, isyrno, row, err);
	if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
	if (row->status<0) continue;  /* Skip deselected records */
	
	/* Finished antenna? */
	if (row->antennaNo < ant) continue;
	if (row->antennaNo > ant) break;
	
	/* See if wanted. */
	want = ObitUVSelWantSour (sel, row->SourID);
	
	/* Check subarray */
	want = want  &&  (row->SubA == sub);
	
	/* Not all antennas wanted */
	want = want  &&  (row->antennaNo == ant);
	if (want) {
	  /* See if flagged value */
	  bad = row->PwrSum1[iif-1] == fblank;
	  bad2 = (SYTable->nPol <= 1)  ||  (row->PwrSum2[iif-1] == fblank);
	  if (numtim < maxtim) {
	    if (numtim == 0) timoff = row->Time;
	    wrktim[numtim] = row->Time - timoff;
	    wrkrec[numtim] = isyrno;
	    if (bad) {
	      work2[numtim] = fblank;
	      work4[numtim] = 0.0;
	    } else {
	      work2[numtim] = row->PwrSum1[iif-1];
	      work4[numtim] = 1.0;
	    } 
	    if (bad2) {
	      work3[numtim] = fblank;
	      work5[numtim] = 0.0;
	    } else {
	      work3[numtim] = row->PwrSum2[iif-1];
	      work5[numtim] = 1.0;
	    } 
	    numtim++;
	  } 
	} 
      } /* end loop reading */;
      save = isyrno - 1; /* How far did we get? */
      
      if (numtim <= 0) goto endAnt;  /* Catch anything? */
      
      /* Smooth as requested */
      /*ObitUVSolnSmooMWF (stpsum, 1.0, wrktim, work2, work4, numtim, work1, 
			 work6, work7, work8, FALSE);*/
      SYsmoIt ("MWF", stpsum, alpha, mxpsum, wrktim, work2, work4, numtim, 
	       work1, work6, work7, work8, FALSE, nThreads, args);
      /* Second Poln? */
      if (SYTable->nPol > 1) {
	/*ObitUVSolnSmooMWF (stpsum, 1.0, wrktim, work3, work5, numtim, work2, 
	  work4, work7, work8, FALSE);*/
	SYsmoIt ("MWF", stpsum, alpha, mxpsum, wrktim, work3, work5, numtim, 
	      work2, work4, work7, work8, FALSE, nThreads, args);
      } 
      
      /* Apply Clip */
     for (itime=0; itime<numtim; itime++) { /* loop 200 */
       /* either poln blanked? */
       if (work6[itime]<=0.0) continue;
       if ((SYTable->nPol > 1) &&  (work4[itime]<=0.0)) continue;
       isyrno = wrkrec[itime];
       ObitTableSYReadRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
       
       /* Blank */
       if (work6[itime]<=0.0) row->PwrSum1[iif-1] = fblank;
       if (work4[itime]<=0.0) row->PwrSum2[iif-1] = fblank;

       /* Rewrite record */
       ObitTableSYWriteRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
     } /* end loop clipping  */;
    endAnt: fstrec = save;
    } /* end loop  of antenna loop */
    
    /* Cleanup */
    row = ObitTableSYRowUnref (row);  
  } /* end clipPSum */

/*----------------------------------------------------------------------- */
/*  Clip Gain                                                             */
/*   Input:                                                               */
/*      SYTable  Table to clip                                            */
/*      iif      IF to operate on, 1-rel                                  */
/*      sel      Selector telling which data wanted                       */
/*      stgain   Smoothing time (day)                                     */
/*      alpha    0 -> 1 = pure boxcar -> pure MWF (ALPHA of the data      */
/*               samples are discarded and the rest averaged).            */
/*      mxgain   Max. excursion from smoothed value                       */
/*      sub      Subarray (1-rel)                                         */
/*      mxtim    Number of elements in wrkrec, wrktim, work1-8            */
/*      wrktim, wrkrec, work* work arrays                                 */
/*      nThreads number of threads to use for MWF                         */
/*      args     thread argument array                                    */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void clipGain (ObitTableSY *SYTable, ObitUVSel *sel, olong iif, 
	       ofloat alpha, ofloat stgain, ofloat mxgain, olong sub, 
	       olong maxtim, olong *wrkrec, ofloat *wrktim,
	       ofloat *work1, ofloat *work2, ofloat *work3, ofloat *work4, 
	       ofloat *work5, ofloat *work6, ofloat *work7, ofloat *work8, 
	       olong nThreads, MednFuncArg** args, ObitErr* err) {
  ObitTableSYRow *row=NULL;
  ofloat fblank =  ObitMagicF();
  olong   loopr, loopa, numtim, ant, fstrec, nleft, save, itime;
  olong  isyrno=0;
  gboolean   bad, bad2, want;
  odouble timoff=0.0;
  gchar *routine = "clipGain";

  /* Create Row */
  row = newObitTableSYRow (SYTable);
  /* Attach row to output buffer */
  ObitTableSYSetRow (SYTable, row, err);
  fstrec = 0;  /* Record number read in table */

  /* Loop over antenna */
  for (loopa= 1; loopa<=SYTable->nAnt; loopa++) {
    ant = loopa;
    /* Want this antenna? */
    if (!ObitUVSelWantAnt (sel, ant)) continue;

    /* Set pointers, counters */
    numtim = 0;
    nleft = SYTable->myDesc->nrow - fstrec;  /* How many rows? */

    /* Loop in time, reading */
    for (loopr= 1; loopr<=nleft; loopr++) { /* loop 100 */
      isyrno = fstrec + loopr;
      ObitTableSYReadRow (SYTable, isyrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
      if (row->status<0) continue;  /* Skip deselected records */
    
      /* Finished antenna? */
      if (row->antennaNo > ant) break;

      /* See if wanted. */
      want = ObitUVSelWantSour (sel, row->SourID);

      /* Check subarray */
      want = want  &&  (row->SubA == sub);

      /* Not all antennas wanted */
      want = want  &&  (row->antennaNo == ant);
      if (want) {

	/* See if flagged value */
	bad = row->Gain1[iif-1] == fblank;
	bad2 = (SYTable->nPol <= 1)  ||  (row->Gain2[iif-1] == fblank);

	/* add to work arrays */	
	if (numtim < maxtim) {
	  if (numtim == 0) timoff = row->Time;
	  wrktim[numtim] = row->Time - timoff;
	  wrkrec[numtim] = isyrno;
	  if (bad) {
	    work2[numtim] = fblank;
	    work4[numtim] = 0.0;
	  } else {
	    work2[numtim] = row->Gain1[iif-1];
	    work4[numtim] = 1.0;
	  } 
	  if (bad2) {
	    work3[numtim] = fblank;
	    work5[numtim] = 0.0;
	  } else {
	    work3[numtim] = row->Gain2[iif-1];
	    work5[numtim] = 1.0;
	  } 
	  numtim++;
	} 
      } 
    } /* end reading loop */;
    save = isyrno - 1; /* How far did we get? */

    if (numtim <= 0) goto endAnt;  /* Catch anything? */

    /* Smooth/clip as requested */
    /*ObitUVSolnSmooMWF (stgain, 1.0, wrktim, work2, work4, numtim, work1, 
      work6, work7, work8, FALSE);*/
    SYsmoIt ("MWF", stgain, alpha, mxgain, wrktim, work2, work4, numtim, 
	     work1, work6, work7, work8, FALSE, nThreads, args);
    /* Second Poln? */
    if (SYTable->nPol > 1) {
      /*ObitUVSolnSmooMWF (stpdif, 1.0, wrktim, work3, work5, numtim, work2, 
	work4, work7, work8, FALSE);*/
      SYsmoIt ("MWF", stgain, alpha, mxgain, wrktim, work3, work5, numtim, 
	       work2, work4, work7, work8, FALSE, nThreads, args);
    } 

    /* Apply Clip */
     for (itime=0; itime<numtim; itime++) { /* loop 200 */
       /* either poln blanked? */
       if (work6[itime]<=0.0) continue;
       if ((SYTable->nPol > 1) &&  (work4[itime]<=0.0)) continue;
       isyrno = wrkrec[itime];
       ObitTableSYReadRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
       
       /* Blank */
       if (work6[itime]<=0.0) row->Gain1[iif-1] = fblank;
       if (work4[itime]<=0.0) row->Gain2[iif-1] = fblank;

       /* Rewrite record */
       ObitTableSYWriteRow (SYTable, isyrno, row, err);
       if (err->error) Obit_traceback_msg (err, routine, SYTable->name);
     } /* end loop clipping  */
  endAnt: fstrec = save;
  } /* end loop  of antenna loop L600: */;

  /* Cleanup */
  row = ObitTableSYRowUnref (row);  

} /* end clipGain */

/**
 * Routine to smooth values in an open SY table.  
 * All poln present are smoothed but only one IF.
 * \param SYTab  SY table object; must be opened/closed externally
 * \param smoFunc  Smoothing function: 'MWF', 'GAUS', else BOX 
 * \param smoType  Type of data to smooth
 *                "DIF ", "ALL ", "    " = "ALL " [def "DIF "]
 * \param alpha    Alpha clip for MWF (0 -> box, 1 -> pure MWF) 
 * \param smoParm  Smoothing time in days for:
 *                 PwrDif, PwrSum, Gain
 *                 0=>fill in for blanked only. 
 * \param iif      Desired IF (0-rel)
 * \param sub      Desired subarray (1-rel)
 * \param nxt      Number of times allowed in wrk 
 * \param work1    Work buffer (nxt*16) 
 * \param work2    Work area >= (nxt*m)  (m=2 BOX, 3 GAUS, 4 MWF) 
 * \param doBlank  replace blanked values with interpolated values?
 * \param nThreads number of threads to use for MWF
 * \param args     thread argument array
 * \param err      Error/message stack, returns if error.
 */
void static
SYSmooth (ObitTableSY *SYTab, gchar* smoFunc, gchar* smoType, ofloat alpha, 
	  ofloat *smoParm, olong iif, olong sub, 
	  olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, 
	  olong nThreads, MednFuncArg** args, ObitErr* err) 
{
  olong   loopa, numtim, ant, numrec, nleft, isyrno=0, itime, 
    n1good, n2good, i, numif, numpol, numant;
  olong loopr, fstrec, save;
  ofloat    stdif=0.0, stsum=0.0, stgain, weight, fblank =  ObitMagicF();
  gboolean  need2, dodif, dosum=FALSE, dogain=FALSE;
  odouble   timoff=0.0;
  ObitIOCode retCode;
  ObitTableSYRow *row=NULL;
  gchar *routine = "SYSmooth";
  
  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSYIsA(SYTab));
  
  /* Get number of records in table */
  numrec = SYTab->myDesc->nrow;
  if (numrec <= 0) return;   /* bail if nothing */
  
  /* Get descriptive info */
  numif  = SYTab->nIF;
  numpol = SYTab->nPol;
  numant = SYTab->nAnt;
  
 
  /* Are there 2 polarizations? */
  need2 = numpol>1;

  /* Only dif ? */
  if (!strncmp(smoType, "DIF ",4)) {
    dodif   = TRUE;
    dosum   = FALSE;
    dogain  = FALSE;
    stdif   = smoParm[0];
    stsum   = 0.0;
    stgain  = 0.0;

    /* All? */
  } else if (!strncmp(smoType, "ALL ",4)) {
    dodif   = TRUE;
    dosum   = TRUE;
    dogain  = TRUE;
    stdif   = smoParm[0];
    stsum   = smoParm[1];
    stgain  = smoParm[2];
  }

  /* Create Row */
  row = newObitTableSYRow (SYTab);
  /* Attach row to output buffer */
  ObitTableSYSetRow (SYTab, row, err);
  fstrec = 0;  /* Record number read in table */
  
  /* Loop over antenna */
  for (loopa= 1; loopa<=numant; loopa++) {
    ant = loopa;
    /* Set pointers, counters */
    numtim = 0;
    nleft = SYTab->myDesc->nrow - fstrec;  /* How many rows? */
    n1good = 0;
    n2good = 0;
    /* Loop in time, reading */
    for (loopr=1; loopr<=nleft; loopr++) { /* loop 100 */
      isyrno = fstrec + loopr;
      retCode = ObitTableSYReadRow (SYTab, isyrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Finished antenna? */
      if (row->antennaNo < ant) continue; /* Shouldn't happen */
      if (row->antennaNo > ant) break;

      /* Want this record */
      if ((row->SubA == sub)  &&  (row->antennaNo == ant)) {

	/* Put in buffer */
	if (numtim >= nxt) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Exceed time limit of %d for %s", 
			 routine, nxt, SYTab->name);
	  row = ObitTableSYRowUnref(row); /* delete row object */
	  return;
	} 
	/* Work1 usage :
	   0 = dif pol 1
	   1 = sum pol 1
	   2 = gain pol 1
	   3 = Weight pol 1
	   4 = dif pol 2
	   5 = sum pol 2
	   6 = gain pol 2
	   7 = Weight pol 2
	   8 = Time(day) relative to first
	   9 = row number
	*/
	
	if (numtim == 0) timoff = row->Time;  /* First time */
	work1[8*nxt+numtim] = row->Time - timoff;
	work1[9*nxt+numtim] = (ofloat)isyrno;
	/*USE WEIGHT TO BLANK CRAZIES;*/
	weight = 1.0;
	/* First polarization */
	if (row->PwrDif1[iif]!=fblank) {
	  work1[0*nxt+numtim]  = row->PwrDif1[iif];
	  work1[1*nxt+numtim]  = row->PwrSum1[iif];
	  work1[2*nxt+numtim]  = row->Gain1[iif];
	  work1[3*nxt+numtim]  = weight;
	  n1good = n1good + 1;
	} else {
	  work1[0*nxt+numtim]  = fblank;
	  work1[1*nxt+numtim]  = fblank;
	  work1[2*nxt+numtim]  = fblank;
	  work1[3*nxt+numtim]  = fblank;
	}
	
	if (need2) {  /* Second polarization */
	  /*USE WEIGHT TO BLANK CRAZIES;*/
	  weight = 1.0;
	if (row->PwrDif2[iif]!=fblank) {
	    work1[4*nxt+numtim]  = row->PwrDif2[iif];
	    work1[5*nxt+numtim]  = row->PwrSum2[iif];
	    work1[6*nxt+numtim]  = row->Gain2[iif];
	    work1[7*nxt+numtim]  = weight;
	    n2good = n2good + 1;
	  } else {
	    work1[4*nxt+numtim]  = fblank;
	    work1[5*nxt+numtim]  = fblank;
	    work1[6*nxt+numtim]  = fblank;
	    work1[7*nxt+numtim]  = fblank;
	  } 
	} /* end second polarization */
      } /* end if want record */ 
      numtim++;   /* count times */
    } /* end loop */

    save = isyrno - 1; /* How far did we get? */
    if (numtim <= 0) goto endAnt;  /* Catch anything? */
    
    /* Smooth as requested */
    if (n1good > 0) { /* First polarization */
      if (dodif) {  /* Diff */
	SYsmoIt (smoFunc, stdif, alpha, -1.0, &work1[8*nxt], &work1[0*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[0*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!(dosum||dogain))) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
     }
      if (dosum) {  /* Sum */
	SYsmoIt (smoFunc, stsum, alpha, -1.0, &work1[8*nxt], &work1[1*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[1*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!dogain)) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
      }
      if (dogain) {  /* Gain */
	SYsmoIt (smoFunc, stsum, alpha, -1.0, &work1[8*nxt], &work1[2*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[2*nxt+i] = work2[i];
	/* Save deblanked weights  */
	if (doBlank) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
     }
    } /* end first polarization */
    
    if (n2good > 0) {  /* Second polarization */
      if (dodif) {  /* Diff */
	SYsmoIt (smoFunc, stdif, alpha, -1.0, &work1[8*nxt], &work1[4*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[4*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!(dosum||dogain))) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
      }
      if (dosum) {  /* Sum */
	SYsmoIt (smoFunc, stsum, alpha, -1.0, &work1[8*nxt], &work1[5*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[5*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!dogain)) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
     }
      if (dogain) {  /* Gain */
	SYsmoIt (smoFunc, stsum, alpha, -1.0, &work1[8*nxt], &work1[6*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank, nThreads, args);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[6*nxt+i] = work2[i];
	/* Save deblanked weights  */
	if (doBlank) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
     }
   } /* end second polarization */
    
    /* Replace with smoothed values */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isyrno = (olong)(work1[9*nxt+itime]+0.5);
      retCode = ObitTableSYReadRow (SYTab, isyrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Update */
      /* weights zero rather than fblank */
      if (work1[3*nxt+itime]==fblank) work1[3*nxt+itime] = 0.0;
      if (work1[3*nxt+itime]>0.0) {
	row->PwrDif1[iif]   = work1[0*nxt+itime];
	row->PwrSum1[iif]   = work1[1*nxt+itime];
	row->Gain1[iif]     = work1[2*nxt+itime];
      } else {  /* Datum bad */
	row->PwrDif1[iif]   = fblank;
	row->PwrSum1[iif]   = fblank;
	row->Gain1[iif]     = fblank;
      }
      if (need2) {
	/* weights zero rather than fblank */
	if (work1[7*nxt+itime]==fblank) work1[7*nxt+itime] = 0.0;
	if (work1[7*nxt+itime] > 0.0) {
	  row->PwrDif2[iif]   = work1[4*nxt+itime];
	  row->PwrSum2[iif]   = work1[5*nxt+itime];
	  row->Gain2[iif]     = work1[6*nxt+itime];
	} else {  /* Datum bad */
	  row->PwrDif2[iif]   = fblank;
	  row->PwrSum2[iif]   = fblank;
	  row->Gain2[iif]     = fblank;
	}
      }
      
      /* Rewrite record */
      retCode = ObitTableSYWriteRow (SYTab, isyrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
    } /* end loop rewriting smoothed solutions L200: */;
    /* First SY number of next antenna */
    
    /* End of antenna loop */
  endAnt: fstrec = save;
  } /* end loop over antennas*/;

  row = ObitTableSYRowUnref(row); /* delete row object */
} /* end of routine SYSmooth */ 
/**
  * Routine to call appropriate smoothing routine.  Magic value blanking  
 * is supported.  
 * Routine adopted from the AIPSish 
 * 31DEC02/APL/PGM/NOTST/SNSMO.FOR/SNSMSM  
 * \param smmeth  Method 'BOX','MWF', 'GAUS', unknown = 'BOX' 
 * \param width   Smoothing time (days) 
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \param clip    Clipping level in sigma
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      [out] Smoothed values. 
 * \param ws      [out] Smoothed weights 
 * \param yor     Scratch 
 * \param wor     Scratch 
 * \param doBlank replace blanked values with interpolated values.
 * \param nThreads number of threads to use for MWF
 * \param args     thread argument array
 */
void static
SYsmoIt (gchar* smmeth, ofloat width, ofloat alpha, ofloat clip, 
	 ofloat* x, ofloat *y, ofloat *w, olong n, 
	 ofloat* ys, ofloat* ws, ofloat *wrk1, ofloat *wrk2, gboolean doBlank, 
	 olong nThreads, MednFuncArg** args) 
{
  ObitThreadFunc Func;
  olong iTh, nValPth;
  /* Any work to do? */
  if (n <= 0) return;

  if (!strncmp (smmeth, "BOX",3)) {
    Func = ThreadBoxSmo;
  } else if (!strncmp (smmeth, "MWF",3)) {
    Func = ThreadMWFSmo;
  } else if (!strncmp (smmeth, "Gaus",4)) {
    Func = ThreadGausSmo;
  } else {   /* Default = box */
    Func = ThreadBoxSmo;
  }
  
  /* Prepare thread arguments */
  nValPth = n/nThreads;
  for (iTh=0; iTh<nThreads; iTh++) {
    args[iTh]->ntime = n;
    args[iTh]->lo = iTh*nValPth;
    args[iTh]->hi = args[iTh]->lo + nValPth - 1;
    if (iTh==nThreads)  args[iTh]->hi = n-1;
    args[iTh]->alpha   = alpha;
    args[iTh]->clip    = clip;
    args[iTh]->smoTime = width;
    args[iTh]->times   = x;
    args[iTh]->vals    = y;
    args[iTh]->out     = ys;
    args[iTh]->weight  = w;
    args[iTh]->outWt   = ws;
  }
  /* Do operation */
  ObitThreadIterator (args[0]->thread, nThreads, 
		      Func, (gpointer**)args);
} /* end of routine SYsmoIt */ 

/**
 * Routine to determine number of times in an open SY table.  
 * \param SYTab    SY table object 
 * \param isub     Subarray number, 0=>1 
 * \param err      Error/message stack, returns if error.
 * \return number of times
 */
olong SYCountTime (ObitTableSY *SYTab, olong isub, ObitErr* err) 
{
  olong  loop, sub, count=0;
  ofloat lastTime;
  ObitTableSYRow *row=NULL;
  gchar *routine = "SYCountTime";

  /* Error checks */
  if (err->error) return count;  /* previous error? */
  g_assert(ObitTableSYIsA(SYTab));

  /* Subarray */
  sub = MAX (1, isub);
  lastTime = -1.0e20;
  count = 0;
  
  /* Open table */
  ObitTableSYOpen (SYTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, SYTab->name, count);

  /* Create Row */
  row = newObitTableSYRow (SYTab);
  /* Loop through table */
  for (loop=1; loop<=SYTab->myDesc->nrow; loop++) {

    ObitTableSYReadRow (SYTab, loop, row, err);
    if (err->error) break;
    if (row->status<0) continue;  /* Skip deselected record */

    /* Right subarray? */
    if ((row->SubA!=sub) && (row->SubA>0)) continue;

    /* Count times - only allow epsilon time difference */
    if (row->Time>(lastTime+0.0005*row->TimeI)) {
      lastTime = row->Time;
      count++;
    }
  } /* end loop  */

  ObitTableSYClose (SYTab, err);
  row = ObitTableSYRowUnref(row); /* delete row object */
  if (err->error) Obit_traceback_val (err, routine, SYTab->name, count);

  return count;
} /* end of routine timeCount */ 

/**
 * Make arguments for a Threaded ThreadFAFunc?
 * \param thread     ObitThread object to be used
 * \param lenwrk     Length of work array
 * \param ThreadArgs[out] Created array of FAFuncArg, 
 *                   delete with KillFAFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeSmoFuncArgs (ObitThread *thread, olong lenwrk, 
			     MednFuncArg ***ThreadArgs)

{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(MednFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(MednFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->clip     = -1.0;
    (*ThreadArgs)[i]->thread   = ObitThreadRef(thread);
    (*ThreadArgs)[i]->vals     = NULL;
    (*ThreadArgs)[i]->out      = NULL;
    (*ThreadArgs)[i]->weight   = NULL;
    (*ThreadArgs)[i]->outWt    = NULL;
    (*ThreadArgs)[i]->times    = NULL;
    if (lenwrk>0)
      (*ThreadArgs)[i]->work   = g_malloc0(lenwrk*sizeof(ofloat));
    else
      (*ThreadArgs)[i]->work   = NULL;
    (*ThreadArgs)[i]->ithread  = i;
  }

  return nThreads;
} /*  end  MakeSmoFuncArgs*/

/**
 * Delete arguments for ThreadFAFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of FAFuncArg
 */
static void KillSmoFuncArgs (olong nargs, MednFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->work)   g_free(ThreadArgs[i]->work);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillSmoFuncArgs */
/**
 * Determine alpha median value of a ofloat array
 * \param n       Number of points
 * \param value   Array of values
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \return alpha median value
 */
static ofloat MedianLevel (olong n, ofloat *value, ofloat alpha)
{
  ofloat out=0.0;
  ofloat fblank = ObitMagicF();
  ofloat beta, sum;
  olong i, i1, i2, count;

  if (n<=0) return out;

  /* Sort to ascending order */
  qsort ((void*)value, n, sizeof(ofloat), compare_ofloat);

  out = value[n/2];

  beta = MAX (0.05, MIN (0.95, 1.0-alpha)) / 2.0; /*  Average around median factor */

  /* Average around the center */
  i1 = MAX (0, (n/2)-(olong)(beta*n+0.5));
  i2 = MIN (n, (n/2)+(olong)(beta*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += value[i];
	count++;
      }
    }
    if (count>0) out = sum / count;
  }
   
  return out;
} /* end MedianLevel */

/**
 * Determine robust RMS value of a ofloat array about mean
 * Use center 90% of points, excluding at least one point from each end
 * \param n       Number of points, needs at least 4
 * \param value   Array of values assumed sorted
 * \param mean    Mean value of value
 * \return RMS value, fblank if cannot determine
 */
ofloat MedianSigma (olong n, ofloat *value, ofloat mean)
{
  ofloat fblank = ObitMagicF();
  ofloat out;
  ofloat sum;
  olong i, i1, i2, count;

  out = fblank;
  if (n<=4) return out;
  if (mean==fblank) return out;

  /* Get RMS around the center 90% */
  i1 = MAX (1,   (n/2)-(olong)(0.45*n+0.5));
  i2 = MIN (n-1, (n/2)+(olong)(0.45*n+0.5));

  if (i2>i1) {
    sum = 0.0;
    count = 0;
    for (i=i1; i<i2; i++) {
      if (value[i]!=fblank) {
	sum += (value[i]-mean)*(value[i]-mean);
	count++;
      }
    }
    if (count>1) out = sqrt(sum / (count-1));
  }
   
  return out;
} /* end MedianSigma */

/**
 * Thread alpha median window smooth an array with optional clip
 * Callable as thread
 * \param arg Pointer to MednFuncArg argument with elements:
 * \li ntime    total number of values in arrays
 * \li lo       first (0-rel) value in array to smooth
 * \li hi       highest (0-rel) value in array to smooth
 * \li alpha    control averaging of center values
 * \li clip     clipping level in sigma, <=0 => no clip
 * \li smoTime  Smoothing time
 * \li vals     Array to smooth
 * \li out      output smoothed array
 * \li weight   Array of input weights (0 or 1)
 * \li outWt    Array of output weights (0 or 1)
 * \li times    times of data samples, <-1000 => ignore
 * \li work     work array for sorting
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadMWFSmo (gpointer arg)
{
  /* Get arguments from structure */
  MednFuncArg *largs = (MednFuncArg*)arg;
  olong  ntime    = largs->ntime;
  olong  lo       = largs->lo;
  olong  hi       = largs->hi;
  ofloat alpha    = largs->alpha;
  ofloat clip     = largs->clip;
  ofloat smoTime  = largs->smoTime/2;
  ofloat *vals    = largs->vals;
  ofloat *out     = largs->out;
  ofloat *weight  = largs->weight;
  ofloat *outWt   = largs->outWt;
  ofloat* times   = largs->times;   
  ofloat* work    = largs->work;
 
  /* local */
  olong i, j, count;
  ofloat time, old, RMS, fblank = ObitMagicF();

  if (hi<lo) goto finish;
  /* Clipping? */
  if (clip>0.0) {  /* Clip */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      for (j=0; j<ntime; j++) {
	if ((fabs(time-times[j])<smoTime) && (weight[j]>0.0)) work[count++] = vals[j];
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      old = out[i];   /* Save initial value */
      if (count>3) {
	out[i]   = MedianLevel (count, work, alpha);
	RMS      = MedianSigma (count, work, out[i]);
	/* Clipping */
	if ((RMS!=fblank)&&(fabs(old-out[i])<RMS*clip)) outWt[i] = 1.0;
	else                                            outWt[i] = 0.0;
      } else {
	out[i]  = vals[i]*weight[i];
	outWt[i]= weight[i];
      }
      
    } /* End smoothing loop */
  } else {         /* no Clip */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      for (j=0; j<ntime; j++) {
	if ((fabs(time-times[j])<smoTime) && (weight[j]>0.0)) work[count++] = vals[j];
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      if (count>3) {
	out[i]   = MedianLevel (count, work, alpha);
	outWt[i] = 1.0;
      } else {
	out[i]  = vals[i]*weight[i];
	outWt[i]= weight[i];
      }
      
    } /* End smoothing loop */
  } /* end no clip */
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadMWFSmo */
/**
 * Thread box smooth an array with optional clip
 * Callable as thread
 * \param arg Pointer to MednFuncArg argument with elements:
 * \li ntime    total number of values in arrays
 * \li lo       first (0-rel) value in array to smooth
 * \li hi       highest (0-rel) value in array to smooth
 * \li clip     clipping level in sigma, <=0 => no clip
 * \li smoTime  Smoothing time
 * \li vals     Array to smooth
 * \li out      output smoothed array
 * \li weight   Array of input weights (0 or 1)
 * \li outWt    Array of output weights (0 or 1)
 * \li times    times of data samples, <-1000 => ignore
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadBoxSmo (gpointer arg)
{
  /* Get arguments from structure */
  MednFuncArg *largs = (MednFuncArg*)arg;
  olong  ntime    = largs->ntime;
  olong  lo       = largs->lo;
  olong  hi       = largs->hi;
  ofloat smoTime  = largs->smoTime/2;
  ofloat clip     = largs->clip;
  ofloat *vals    = largs->vals;
  ofloat *out     = largs->out;
  ofloat *weight  = largs->weight;
  ofloat *outWt   = largs->outWt;
  ofloat* times   = largs->times;   
 
  /* local */
  olong i, j, count;
  ofloat time, sum, sum2, RMS, old;

  if (hi<lo) goto finish;

  /* Clipping? */
  if (clip>0.0) {  /* Clip */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      sum   = sum2 = 0.0;
      for (j=0; j<ntime; j++) {
	if ((fabs(time-times[j])<smoTime) && (weight[j]>0.0)) {
	  sum  += vals[j];
	  sum2 += vals[j]*vals[j];
	  count++;
	}
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      old = out[i];
      if (count>3) {
	out[i]   = sum / count;
	outWt[i] = 1.0;
	RMS      = sqrt(MIN (0.01*out[i], (sum2 / count) - out[i]*out[i]));
	/* Clipping */
	if (fabs(old-out[i])<RMS*clip) outWt[i] = 1.0;
	else                           outWt[i] = 0.0;
      } else {
	out[i]   = vals[i]*weight[i];
	outWt[i] = weight[i];
      }
    } /* End smoothing loop */
  } else { /* no clip */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      sum   = sum2 = 0.0;
      for (j=0; j<ntime; j++) {
	if ((fabs(time-times[j])<smoTime) && (weight[j]>0.0)) {
	  sum  += vals[j];
	  sum2 += vals[j]*vals[j];
	  count++;
	}
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      if (count>3) {
	out[i]   = sum / count;
	outWt[i] = 1.0;
      } else {
	out[i]   = vals[i]*weight[i];
	outWt[i] = weight[i];
      }
    } /* End smoothing loop */
  } /* end no clipping */
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadBoxSmo */
/**
 * Thread Gaussian smooth an array with optional clip
 * Callable as thread
 * \param arg Pointer to MednFuncArg argument with elements:
 * \li ntime    total number of values in arrays
 * \li lo       first (0-rel) value in array to smooth
 * \li hi       highest (0-rel) value in array to smooth
 * \li clip     clipping level in sigma, <=0 => no clip
 * \li smoTime  Smoothing time, twice sigma of Gaussian.
 * \li vals     Array to smooth
 * \li out      output smoothed array
 * \li weight   Array of input weights (0 or 1)
 * \li outWt    Array of output weights (0 or 1)
 * \li times    times of data samples, <-1000 => ignore
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadGausSmo (gpointer arg)
{
  /* Get arguments from structure */
  MednFuncArg *largs = (MednFuncArg*)arg;
  olong  ntime    = largs->ntime;
  olong  lo       = largs->lo;
  olong  hi       = largs->hi;
  ofloat smoTime  = largs->smoTime;
  ofloat clip     = largs->clip;
  ofloat ivar     = 1.0 / (0.25*largs->smoTime*largs->smoTime);
  ofloat *vals    = largs->vals;
  ofloat *out     = largs->out;
  ofloat *weight  = largs->weight;
  ofloat *outWt   = largs->outWt;
  ofloat* times   = largs->times;   
 
  /* local */
  olong i, j, count;
  ofloat time, sum, sum2, sum3, RMS, g, dt, old;

  if (hi<lo) goto finish;

   /* Clipping? */
  if (clip>0.0) {  /* Clip */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      sum   = sum2 = sum3 = 0.0;
      for (j=0; j<ntime; j++) {
	dt = fabs(time-times[j]);
	if ((dt<smoTime) && (weight[j]>0.0)) {
	  g = exp(-ivar*dt*dt);
	  sum  += g*vals[j];
	  sum2 += vals[j];
	  sum3 += vals[j]*vals[j];
	  count++;
	}
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      old = out[i];
      if (count>3) {
	out[i]   = sum / count;
	outWt[i] = 1.0;
	sum2 /= count;
	RMS      = sqrt(MIN (0.01*out[i], (sum3 / count) - sum2*sum2));
	/* Clipping */
	if (fabs(old-out[i])<RMS*clip) outWt[i] = 1.0;
	else                           outWt[i] = 0.0;
       } else {
	out[i]   = vals[i]*weight[i];
	outWt[i] = weight[i];
      }
    } /* End smoothing loop */
  } else { /* no clipping */
    /* Loop over values to smooth */
    for (i=lo; i<=hi; i++) {
      time = times[i];
      /* Accumulate values to smooth */
      count = 0;
      sum   = 0.0;
      for (j=0; j<ntime; j++) {
	dt = fabs(time-times[j]);
	if ((dt<smoTime) && (weight[j]>0.0)) {
	  g = exp(-ivar*dt*dt);
	  sum += g*vals[j];
	  count++;
	}
	if (times[j]>(time+smoTime)) break;
      } /* end accumulating */
      if (count>3) {
	out[i]   = sum / count;
	outWt[i] = 1.0;
      } else {
	out[i]   = vals[i]*weight[i];
	outWt[i] = weight[i];
      }
    } /* End smoothing loop */
  } /* end no clipping */
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadGausSmo */
/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_ofloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2)      out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_ofloat */


