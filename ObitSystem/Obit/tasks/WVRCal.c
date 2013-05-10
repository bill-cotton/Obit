/* $Id$  */
/* Convert an ALMA WVR dataset to an SN table              */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011-2013                                          */
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

#include "ObitUV.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitWVRCoef.h"
#include "ObitTableSUUtil.h"
#include "ObitTableANUtil.h"
#include "ObitUVSoln.h"
#include "ObitAntennaList.h"
#include "ObitSourceList.h"
#include "ObitPrecess.h"
#include "ObitThread.h"
/* libALMAWVR stuff */
#ifdef HAVE_WVR  /* Only if libALMAWVR available */
#include "almawvr/almaabs_c.h"
#endif /* HAVE_WVR */
  /* Speed of light */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
  /* CI = 1/speed of light */
#ifndef CI
#define CI 1.0 / VELIGHT
#endif

/*---------------Private structures----------------*/
/* Threaded function argument */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* thread number  */
  olong        ithread;
/* libALMAWVR stuff */
#ifdef HAVE_WVR  /* Only if libALMAWVR available */
  /* Input data */
  ALMAAbsInput wvrIn, *wvrInArr[1];
  /* Return results */
  ALMARes_Basic wvrOut, *wvrOutArr[1];
#endif /* HAVE_WVR */
} WVRFuncArg;

/* internal prototypes */
/* Get inputs */
ObitInfoList* WVRCalIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void WVRCalOut (ObitInfoList* outList, ObitErr *err);
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
/* Create output uvdata */
ObitUV* getOutputData (ObitInfoList *myInput, ObitErr *err);
/* Get scan averaged sensitivities of WVR  */
void WVRCalScanCal (ObitUV* inData, ObitUV* outData, ObitWVRCoef *wvrcoef,
		    ObitErr* err);
/* Convert WVR data to SN  */
void WVRCalConvert (ObitUV* inData, ObitUV* outData, ObitWVRCoef *wvrcoef,
		    ObitErr* err);
/* Write history */
void WVRCalHistory (ObitInfoList* myInput, ObitUV* inData, 
		    ObitUV* outData, ObitErr* err);
/* Get SourceList   */
ObitSourceList* GetSourceList (ObitUV* inData, ObitErr* err);
/* Get AntennaList   */
ObitAntennaList* GetAntennaList (ObitUV* inUV, olong suba, ObitErr* err);
/* Scan average next scan and determine sensivities */
gboolean AvgWVR1 (ObitUV* inUV, ofloat solInt, olong numAnt, 
		  ObitAntennaList *AList, ObitSourceList *SList, olong *sid, 
		  ObitWVRCoef *wvrcoef,
		  olong *count, odouble *sumWt, odouble *sumWtT, odouble *sumEl,
		  olong nTh, WVRFuncArg **ThreadArgs, ObitErr* err);
/* Convert WVR data to delays, average */
gboolean AvgWVR2 (ObitUV* inUV, ofloat solInt, olong numAnt, 
		  ObitAntennaList *AList, ObitSourceList *SList,
		  odouble *timec, ofloat *timei, olong *sid, olong *fqid, 
		  olong *refAnt, gboolean *gotAnt, ofloat *delay, ofloat *wt,  
		  ObitWVRCoef *wvrcoef,
		  olong *count, odouble *sumWt, odouble *sumWtT, odouble *sumEl,
		  ObitErr* err);
/* Make threading arguments */
static olong MakeWVRFuncArgs (ObitThread *thread, 
			      WVRFuncArg ***ThreadArgs);
/* Kill threading arguments */
static void KillWVRFuncArgs (olong nargs, WVRFuncArg **ThreadArgs);
/* Thread function */
static gpointer ThreadWVR (gpointer arg);
/** Time to String */
static void T2String (ofloat time, gchar *msgBuf);

/* Program globals */
gchar *pgmName = "WVRCal";       /* Program name */
gchar *infile  = "WVRCal.in" ;   /* File with program inputs */
gchar *outfile = "WVRCal.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* Convert an ALMA WVR dataset to an SN table                             */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *outData=NULL;
  ObitWVRCoef  *wvrcoef=NULL;
  ObitErr      *err= NULL;
  /* Selection controls for uv data */
  gchar *UVParms1[] = {"timeRange", "flagVer", "Antennas", "calSour", "Qual",
		       "calCode", "solInt1", "solnVer", "refAnt", "doFlip",
		      NULL };
  gchar *UVParms1Rename[] = {"timeRange", "flagVer", "Antennas", "Sources", "Qual",
			     "souCode", "solInt", "solnVer", "refAnt", "doFlip",
		      NULL };
  gchar *UVParms2[] = {"timeRange", "flagVer", "Antennas", "Sources", "Qual",
		       "souCode", "solInt2", "solnVer",  "refAnt", "doFlip",
		      NULL };
  gchar *UVParms2Rename[] = {"timeRange", "flagVer", "Antennas", "Sources", "Qual",
			     "souCode", "solInt", "solnVer",  "refAnt", "doFlip",
			     NULL };
   /* Startup - parse command line */
  err = newObitErr();
  myInput = WVRCalIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto cleanup;

  /* Get output uvdata */
  outData = getOutputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create wvrcoef structure */
  wvrcoef = newObitWVRCoef();

  /* Copy control info */
  ObitInfoListCopyListRename (myInput, inData->info,  UVParms1, UVParms1Rename);
  ObitInfoListCopyListRename (myInput, outData->info, UVParms1, UVParms1Rename);

  /* Scan average sensitivities */
  WVRCalScanCal (inData, outData, wvrcoef, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto cleanup;

  /* Copy control info */
  ObitInfoListCopyListRename (myInput, inData->info,  UVParms2, UVParms2Rename);
  ObitInfoListCopyListRename (myInput, outData->info, UVParms2, UVParms2Rename);

  /* Convert WVR to SN */
  WVRCalConvert (inData, outData, wvrcoef, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto cleanup;

  /* History */
  WVRCalHistory (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto cleanup;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto cleanup;
  
  /* cleanup */
cleanup:
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  wvrcoef   = ObitWVRCoefUnref(wvrcoef);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* WVRCalIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "WVRCalIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

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
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end WVRCalIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: WVRCal -input file -output ofile [args]\n");
    fprintf(stderr, "WVRCal Obit task to image/CLEAN data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def WVRCal.in\n");
    fprintf(stderr, "  -output output result file, def WVRCal.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
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
/*     inFile    Str [?]    input FITS uv file name [no def]              */
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
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
  strTemp = "WVRCal.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "WVRCalName";
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

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

 /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "flagVer", OBIT_oint, dim, &itemp, err);
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
  /*ObitInfoType type;
    gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    ofloat ftemp;*/
  gchar *routine = "digestInputs";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Initialize Threading */
  ObitThreadInit (myInput);

} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data                                                        */
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
  ObitIOType   IOType;
  olong         Aseq, disk, cno, nvis;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7],  *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  inData = newObitUV("input UV data");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    IOType = OBIT_IO_AIPS;  /* Save file type */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    /* input AIPS name */
    if (ObitInfoListGetP(myInput, "inName", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aname, strTemp, 13);
    } else { /* Didn't find */
      strncpy (Aname, "No Name ", 13);
    } 
    Aname[12] = 0;

    /* input AIPS class */
    if  (ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    }

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
    /* define object */
    nvis = 1;
    ObitUVSetAIPS (inData, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }

    g_snprintf (Aname, 12, "%s", inFile);  /* Use as root of scratch images */
    IOType = OBIT_IO_FITS;  /* Save file type */
   
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* define object */
    nvis = 1;
    ObitUVSetFITS (inData, nvis, disk, inFile,  err); 
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }

 /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* getOutputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129];
  gchar     *routine = "getOutputData";

  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* outName given? */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

      
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "WVRCal" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "UVSub", 7);

    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

    /* Output file really should exist */
    Obit_retval_if_fail(exist, err, outUV,
			"%s: Output data does not exist", routine);
    
    /* define object */
    nvis = 1000;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to outDisk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* define object */
    nvis = 1000;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS UV data %s on disk %d", outFile, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  ObitErrLog(err); /* Show messages */
  return outUV;
} /* end getOutputUV */

/*----------------------------------------------------------------------- */
/*  Write History for WVRCal                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to copy history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void WVRCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		    ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "timeRange" "Sources",  "Qual", "souCode", "flagVer", 
    "outFile",  "outDisk",  "outName", "outClass", "outSeq",
    "refAnt", "solInt1", "solInt2", "solnVer", "nThread", "doFlip",
    NULL};
  gchar *routine = "WVRCalHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
  } else { /* simply copy history */
    ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end WVRCalHistory  */

/**
 *  Scan Average selected WVR data to sensitivities
 * \param   inUV    ObitUV with WVR data
 * Control parameters are on the info member.
 * \li "solInt"   OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \param   outUV   ObitUV with output SN table
 * \param   wvrcoef Structure with fitted sensitivities
 * \param   err     Obit Error stack
 */
void WVRCalScanCal (ObitUV* inUV, ObitUV* outUV, ObitWVRCoef *wvrcoef,
		    ObitErr* err)
{
  ObitAntennaList *AList=NULL;
  ObitSourceList  *SList=NULL;
  WVRFuncArg **ThreadArgs=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong ifreq, sid, nTh, nextVisBuf, nVisPIO, *count=NULL;
  oint numAnt, suba=1, prtlv;
  ofloat solInt;
  odouble *sumWt=NULL, *sumWtT=NULL, *sumEl=NULL;
  gboolean done, doCalSelect;
  ObitIOCode retCode;
  gchar *routine = "WVRCalScanCal";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));
  
  /* Get Solution interval */
  solInt = 0.0;
  ObitInfoListGetTest(inUV->info, "solInt", &type, dim, &solInt);
  solInt /= 1440.0;  /* Convert to days */

  /* One visibility per access */
  nVisPIO = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, & nVisPIO);
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, & nVisPIO);
  
  /* open UV data  */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine,outUV->name);
  
  /* Make sure some data */
  Obit_return_if_fail((inUV->myDesc->nvis>0), err, 
		      "%s: NO data in %s", routine, inUV->name);

  /* Look like WVR data? Check frequency and number of channels */
  ifreq = inUV->myDesc->jlocf;
  Obit_return_if_fail(((inUV->myDesc->crval[ifreq]>=170.0e9) &&
		       (inUV->myDesc->crval[ifreq]<=190.0e9) &&
		       (inUV->myDesc->inaxes[ifreq]==4)), err, 
		      "%s: NO WVR data in %s", routine, inUV->name);
  /* Check compatability of inUV and outUV */
  ifreq = inUV->myDesc->jlocf;
  Obit_return_if_fail((inUV->myDesc->JDObs==outUV->myDesc->JDObs), err, 
		      "%s: inUV and outUV different dates", routine);
  
  /* Update frequency tables on outUV */
  if (!inUV->myDesc->freqArr) ObitUVGetFreq (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  /* Need array information */
  if (!inUV->myDesc->numAnt)   ObitUVGetSubA (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  
  /* Source list */
  SList = GetSourceList (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  if (outUV->mySel->sources) sid = outUV->mySel->sources[0];
  else                       sid = 1;

  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(outUV->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_return_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);

  /* Antenna list */
  AList = GetAntennaList(inUV, suba, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  
  /* Create arrays */
  numAnt   = inUV->myDesc->numAnt[suba-1];
  count    = g_malloc0(numAnt*sizeof(olong));
  sumWt    = g_malloc0(numAnt*sizeof(odouble));
  sumWtT   = g_malloc0(4*numAnt*sizeof(odouble));
  sumEl    = g_malloc0(numAnt*sizeof(odouble));
  
  /* Get parameters from inUV */
  prtlv = 0;
  ObitInfoListGetTest(inUV->info, "prtLv", &type, dim, &prtlv);
  
  /* Initialize threading */
  nTh = MakeWVRFuncArgs (inUV->thread, &ThreadArgs);

  /* Read first WVR record */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) ObitUVReadSelect (inUV, inUV->buffer, err);
  else ObitUVRead (inUV, inUV->buffer, err);
  if (err->error) goto cleanup;
  
  /* Loop until done */
  done = FALSE;
  nextVisBuf = -1;
  while (!done) {
    /* Read and average next solution interval */
    done =  AvgWVR1 (inUV, solInt, numAnt, AList, SList, &sid,
		    wvrcoef,
		    count, sumWt, sumWtT, sumEl, 
		    nTh, ThreadArgs, err);
    if (err->error) goto cleanup;

    /* Done? */
    if (done) break; /* still have OK data */
    
    /* Messages */
    if (err->prtLv>1) ObitErrLog(err);
    
  } /* end loop processing data */
  
  /* Close uv data */
  ObitUVClose (inUV, err);
  ObitUVClose (outUV, err);
  if (err->error) goto cleanup;
  
  /* DEBUG - print structure */
  if (err->prtLv>=5) {
    ObitWVRCoefPrint (wvrcoef, stderr);
  }

  goto cleanup; /* Cleanup */
 cleanup: 
  KillWVRFuncArgs (nTh, ThreadArgs);
  AList   = ObitAntennaListUnref(AList);
  SList   = ObitSourceListUnref(SList);
  if (count)    g_free(count);
  if (sumWt)    g_free(sumWt);
  if (sumWtT)   g_free(sumWtT);
  if (sumEl)    g_free(sumEl);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  return;
} /* end WVRCalScanCal  */

/**
 *  Convert selected WVR data to SN table
 * \param   inUV    ObitUV with WVR data
 * Control parameters are on the info member.
 * \li "solnVer"  OBIT_int   (1,1,1) Solution (SN) table to write; 0=> create new.
 * \li "solInt"   OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \li "refAnt"   OBIT_int (1,1,1) reference antenna, 0=>pick
 * \li "doFlip"   OBIT_boo (1,1,1) Flip sign of correction?
 * \param   outUV   ObitUV with output SN table
 * \param   wvrcoef Structure with fitted sensitivities
 * \param   err     Obit Error stack
 */
void WVRCalConvert (ObitUV* inUV, ObitUV* outUV, ObitWVRCoef *wvrcoef,
		    ObitErr* err)
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *row=NULL;
  ObitUVSel *sel = NULL;
  ObitAntennaList *AList=NULL;
  ObitSourceList  *SList=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, iAnt, SNver, ifreq, sid, fqid=0;
  olong nextVisBuf, iSNRow, numFreq, cntGood=0, cntPoss=0, cntBad=0;
  olong itemp, nVisPIO, *count=NULL, refAnt=0, refant;
  oint numPol, numIF, numAnt, suba=1, prtlv;
  ofloat solInt, *delay=NULL, *wt=NULL;
  ofloat timei=0.0, phase, fblank = ObitMagicF();
  odouble *sumWt=NULL, *sumWtT=NULL, *sumEl=NULL;
  odouble timec=0.0;
  gboolean done, good, oldSN, *gotAnt, doCalSelect, empty, doFlip=FALSE;
  ObitIOCode retCode;
  gchar *tname;
  /*ofloat T[4];  DEBUG */
  gchar *routine = "WVRCalConvert";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));
  
  /* Get Solution interval */
  solInt = 0.0;
  ObitInfoListGetTest(inUV->info, "solInt", &type, dim, &solInt);
  solInt /= 1440.0;  /* Convert to days */

  /* reference antenna */
  ObitInfoListGetTest(inUV->info, "refAnt", &type, dim, &refAnt);

  /* Flip sign? */
  ObitInfoListGetTest(inUV->info, "doFlip", &type, dim, &doFlip);

  /* One visibility per access */
  nVisPIO = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, & nVisPIO);
  ObitInfoListAlwaysPut(inUV->info, "nVisPIO", OBIT_long, dim, & nVisPIO);
  
  /* open UV data  */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine,outUV->name);
  
  /* Make sure some data */
  Obit_return_if_fail((inUV->myDesc->nvis>0), err, 
		      "%s: NO data in %s", routine, inUV->name);

  /* Look like WVR data? Chech frequency and number of channels */
  ifreq = inUV->myDesc->jlocf;
  Obit_return_if_fail(((inUV->myDesc->crval[ifreq]>=170.0e9) &&
		       (inUV->myDesc->crval[ifreq]<=190.0e9) &&
		       (inUV->myDesc->inaxes[ifreq]==4)), err, 
		      "%s: NO WVR data in %s", routine, inUV->name);
  /* Check compatability of inUV and outUV */
  ifreq = inUV->myDesc->jlocf;
  Obit_return_if_fail((inUV->myDesc->JDObs==outUV->myDesc->JDObs), err, 
		      "%s: inUV and outUV different dates", routine);
  
  /* Update frequency tables on outUV */
  if (!inUV->myDesc->freqArr) ObitUVGetFreq (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  /* Need array information */
  if (!inUV->myDesc->numAnt)   ObitUVGetSubA (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  
  /* Source list */
  SList = GetSourceList (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  if (outUV->mySel->sources) sid = outUV->mySel->sources[0];
  else                       sid = 1;

  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(outUV->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_return_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);
  
  /* Antenna list */
  AList = GetAntennaList(inUV, suba, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Create output - version requested? */
  itemp = 0;
  ObitInfoListGetTest(inUV->info, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  oldSN = SNver > 0;  /* Does SN table already exist? */
  /* 0=> make new */
  itemp = ObitTableListGetHigh (outUV->tableList, "AIPS SN") + 1;
  if (SNver<=0) SNver = itemp;
  tname = g_strconcat ("SN Calibration for: ", outUV->name, NULL);
  if (inUV->myDesc->jlocs>=0)
    numPol = MIN (2, outUV->myDesc->inaxes[outUV->myDesc->jlocs]);
  else numPol = 1;
  if (outUV->myDesc->jlocif>=0)
    numIF  = outUV->myDesc->inaxes[outUV->myDesc->jlocif];
  else numIF  = 1;
  outSoln = newObitTableSNValue(tname, (ObitData*)outUV, &SNver, OBIT_IO_ReadWrite, 
				numPol, numIF, err);
  g_free (tname);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  if (outUV->myDesc->jlocf>=0)
    numFreq = outUV->myDesc->inaxes[outUV->myDesc->jlocf];
  else numFreq = 1;

  /* If SN table previously existed, deselect values about to be redetermined. 
     get information from selector on inUV */
  sel = outUV->mySel;
  if (oldSN) ObitUVSolnDeselSN (outSoln, sel->SubA, sel->FreqID, 
				sel->numberAntList, sel->ants, 
				sel->numberSourcesList, sel->sources, 
				sel->timeRange, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
 
  /* Create arrays */
  numAnt   = inUV->myDesc->numAnt[suba-1];
  gotAnt   = g_malloc0(numAnt*sizeof(gboolean));
  delay    = g_malloc0(numAnt*sizeof(ofloat));
  wt       = g_malloc0(numAnt*sizeof(ofloat));
  count    = g_malloc0(numAnt*sizeof(olong));
  sumWt    = g_malloc0(numAnt*sizeof(odouble));
  sumWtT   = g_malloc0(4*numAnt*sizeof(odouble));
  sumEl    = g_malloc0(numAnt*sizeof(odouble));
  
  /* Get parameters from inUV */
  prtlv = 0;
  ObitInfoListGetTest(inUV->info, "prtLv", &type, dim, &prtlv);
  
  /* Open output table */
  retCode = ObitTableSNOpen (outSoln, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  /* Anything already there? */
  empty = outSoln->myDesc->nrow==0;
  if (empty) {  /* Init if empty */
    outSoln->numAnt = numAnt;  /* Number of antennas */
    outSoln->mGMod  = 1.0;     /* initial mean gain modulus */
  }
  
  /* Create Row */
  row = newObitTableSNRow (outSoln);
  
  /* Attach row to output buffer */
  ObitTableSNSetRow (outSoln, row, err);
  if (err->error) goto cleanup;
  
  /* Initialize solution row */
  row->Time   = 0.0; 
  row->TimeI  = 0.0; 
  row->SourID = sid; 
  row->antNo  = 0; 
  row->SubA   = 0; 
  row->FreqID = 0; 
  row->IFR    = 0.0; 
  row->NodeNo = 0; 
  row->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    row->Real1[i]   = 0.0; 
    row->Imag1[i]   = 0.0; 
    row->Delay1[i]  = 0.0; 
    row->Rate1[i]   = 0.0; 
    row->Weight1[i] = 0.0; 
    row->RefAnt1[i] = 0; 
  }
  if (numPol>1) {
    row->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      row->Real2[i]   = 0.0; 
      row->Imag2[i]   = 0.0; 
      row->Delay2[i]  = 0.0; 
      row->Rate2[i]   = 0.0; 
      row->Weight2[i] = 0.0; 
      row->RefAnt2[i] = 0; 
    }
  }

  /* Read first WVR record */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) ObitUVReadSelect (inUV, inUV->buffer, err);
  else ObitUVRead (inUV, inUV->buffer, err);
  if (err->error) goto cleanup;
  
  /* Loop until done */
  done = FALSE;
  nextVisBuf = -1;
  while (!done) {
    refant = refAnt;
    /* Read and average next solution interval */
    done =  AvgWVR2 (inUV, solInt, numAnt, AList, SList, 
		     &timec, &timei, &sid, &fqid, 
		     &refant, gotAnt, delay, wt, wvrcoef,
		     count, sumWt, sumWtT, sumEl, 
		     err);
    if (err->error) goto cleanup;

    /* Done? */
    if (done) break; /* still have OK data */
    
    /* Messages */
    if (err->prtLv>1) ObitErrLog(err);
    
    /* How many possible solutions */
    for (iAnt= 0; iAnt<numAnt; iAnt++) if (gotAnt[iAnt]) cntPoss += numIF*numPol;
    
    /* Write solutions to SN table */
    /* Common values */
    row->Time   = timec; 
    row->TimeI  = timei; 
    if (sid>0) row->SourID = sid; 
    row->SubA   = suba; 
    row->FreqID = fqid; 
    
    /* Loop over antennas */
    iSNRow = -1;
    for (iAnt= 0; iAnt<numAnt; iAnt++) {
      if (gotAnt[iAnt]) {
	good = FALSE; /* Until proven */
	row->antNo  = iAnt+1; 
	if (doFlip)
	  row->MBDelay1 = -delay[iAnt];
	else
	  row->MBDelay1 = +delay[iAnt];
	for (i=0; i<numIF; i++) {
	  if (delay[iAnt]!=fblank) {
	    phase = 2*G_PI*outUV->myDesc->freqIF[i]*row->MBDelay1;
	    row->Real1[i]   = cos(phase); 
	    row->Imag1[i]   = sin(phase); 
	  } else {
	    row->Real1[i]   = fblank; 
	    row->Imag1[i]   = fblank; 
	  }
	  row->Delay1[i]  = row->MBDelay1;
	  row->Weight1[i] = wt[iAnt]; 
	  row->RefAnt1[i] = refant; 
	  if (wt[iAnt]>0.0) {good = TRUE; cntGood++;}
	  if (wt[iAnt]<=0.0) cntBad++;    /* DEBUG */
	}
	if (numPol>1) {
	  row->MBDelay2 = row->MBDelay1;
	  for (i=0; i<numIF; i++) {
	    row->Real2[i]   = row->Real1[i]; 
	    row->Imag2[i]   = row->Imag1[i]; 
	    row->Delay2[i]  = row->MBDelay2;
	    row->Weight2[i] = wt[iAnt]; 
	    row->RefAnt2[i] = refant;
	    if (wt[iAnt]>0.0) {good = TRUE; cntGood++;}
	    if (wt[iAnt]<=0.0) cntBad++;    /* DEBUG */
	  }
	}
	retCode = ObitTableSNWriteRow (outSoln, iSNRow, row, err);
	if (err->error) goto cleanup;
      } /* end if gotant */
    }
  } /* end loop processing data */
  
  /* Close uv data */
  ObitUVClose (inUV, err);
  ObitUVClose (outUV, err);
  if (err->error) goto cleanup;
  
  /* If suba>1 or not empty at start mark as unsorted */
  if ((suba == 1) && empty) {  
    /* IF subarray 1, new table the sort is time, antenna */
    outSoln->myDesc->sort[0] = outSoln->TimeCol+1;
    outSoln->myDesc->sort[1] = outSoln->antNoCol+1; 
  } else { /* otherwise unsorted */
    outSoln->myDesc->sort[0] = 0;
    outSoln->myDesc->sort[1] = 0;
  }

  /* Close output table */
  retCode = ObitTableSNClose (outSoln, err);
  if (err->error) goto cleanup;
  
  /* Give results */
  if (err->prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, " Write %d records to SN version %d",
		   cntGood, outSoln->tabVer);
  }
  /* Give success rate
  if (err->prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr, " %d of %d possible solutions found",
		   cntGood, cntPoss); 
  }*/

  /* DEBUG
  T[0] = 195.1;  T[1] = 103.5; T[2] = 51.7; T[3] = 39.6;
  fprintf(stderr, "test %g\n",
	  ObitWVRCoefCal (wvrcoef, 0.25, 3, 2, T)); */

  goto cleanup; /* Cleanup */
 cleanup: 
  row     = ObitTableSNRowUnref(row);
  outSoln = ObitTableSNUnref(outSoln);
  AList   = ObitAntennaListUnref(AList);
  SList   = ObitSourceListUnref(SList);
  if (delay)    g_free(delay);
  if (wt)       g_free(wt);
  if (gotAnt)   g_free(gotAnt);
  if (count)    g_free(count);
  if (sumWt)    g_free(sumWt);
  if (sumWtT)   g_free(sumWtT);
  if (sumEl)    g_free(sumEl);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  return;
} /* end WVRCalConvert  */

/**
 *  Get Source List
 * \param   inUV       ObitUV with SU Table
 * \param   err        Obit Error stack
 * \return SourceList, should be Unrefed when done
 */
ObitSourceList* GetSourceList (ObitUV* inUV, ObitErr* err)
{
  olong iver = 1;
  ObitTableSU *SUTable=NULL;
  ObitSourceList  *SList=NULL;
  gchar *routine = "GetSourceList";

  if (err->error) return SList;

  SUTable = newObitTableSUValue ("SUTable", (ObitData*)inUV, &iver, 
			       OBIT_IO_ReadOnly, 0, err);
  if (SUTable) {
    SList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, SUTable->name, SList);
    SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
  } else {  /* Use position from header */
    SList = ObitSourceListCreate ("SList", 1);
    SList->SUlist[0]->equinox = inUV->myDesc->equinox;
    SList->SUlist[0]->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    SList->SUlist[0]->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, SList->SUlist[0]);
  }
  return SList;
} /* end GetSourceList */

/**
 *  Get Antenna List
 * \param   inUV       ObitUV with AN Table
 * \param   suba       Subarray
 * \param   err        Obit Error stack
 * \return AntennaList, should be Unrefed when done
 */
ObitAntennaList* GetAntennaList (ObitUV* inUV, olong suba, ObitErr* err)
{
  olong iver = 1;
  ObitTableAN *ANTable=NULL;
  ObitAntennaList  *AList=NULL;
  gchar *routine = "GetAntennaList";

  if (err->error) return AList;

  iver = MAX (1, suba);
  ANTable = newObitTableANValue ("AList", (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  Obit_retval_if_fail((ObitTableANIsA(ANTable)), err, AList,
		      "%s: Antenna table not found", routine);
  AList = ObitTableANGetList (ANTable, err);
  ANTable = ObitTableANUnref(ANTable);   /* Done with table */
  if (err->error) Obit_traceback_val (err, routine, ANTable->name, AList);

  return AList;
} /* end GetAntennaList */

/**
 *  Scan average selected sources and calculate sensitivities
 * \param   inUV     ObitUV with WVR data
 * \param   solInt   Desired solution interval in days
 * \param   numAnt   Number of antennas
 * \param   AList    Antenna list
 * \param   SList    Source list
 * \param   sid      [in/out] source ID
 * \param   wvrcoef  [out] Structure saving fitted sensitivities
 * \param   count    Work array
 * \param   sumWt    Work array
 * \param   sumWtT   Work array
 * \param   sumEl      Work array
 * \param   nTh        Number of threads
 * \param   ThreadArgs Threading function arguments 
 * \param   err        Obit Error stack
 * \return TRUE if all data processed, could be some OK in output 
 */
gboolean AvgWVR1 (ObitUV* inUV, ofloat solInt, olong numAnt, 
		  ObitAntennaList *AList, ObitSourceList *SList, olong *sid,
		  ObitWVRCoef *wvrcoef,
		  olong *count, odouble *sumWt, odouble *sumWtT, odouble *sumEl,
		  olong nTh, WVRFuncArg **ThreadArgs, ObitErr* err)
{
  gboolean done=TRUE;
#ifdef HAVE_WVR  /* Only if libALMAWVR available */
  olong iAnt, jAnt, iTh, jTh, tsid;
  odouble endTime, begTime;
  ofloat elev, time, weight, base, tr[2], dTdL[4], TObs[4], wt;
  gboolean more, doCalSelect, OK, first;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar tStr[25];
  gchar *routine = "AvgWVR1";
  
  /* error checks */
  if (err->error) return done;
  
  /* Init */
  begTime =  1.0e20;
  endTime = -1.0e20;
  for (iAnt=0; iAnt<numAnt; iAnt++) {
    count[iAnt]      = 0;
    sumWt[iAnt]      = 0.0;
    sumEl[iAnt]      = 0.0;
    sumWtT[iAnt*4+0] = 0.0;
    sumWtT[iAnt*4+1] = 0.0;
    sumWtT[iAnt*4+2] = 0.0;
    sumWtT[iAnt*4+3] = 0.0;
  }
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, 
		      &doCalSelect);
  
  if (inUV->myDesc->ilocsu>=0) 
    *sid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]+0.5);
  more = first = TRUE;
  while (more) {
    /* Info */
    time = inUV->buffer[inUV->myDesc->iloct];
    begTime = MIN (begTime, (odouble)time);
    endTime = MAX (endTime, (odouble)time);
    if (inUV->myDesc->ilocsu>=0) 
      tsid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]+0.5);
    else tsid = *sid;
    /* finished? */
    more    = (endTime<(begTime+solInt)) && (tsid==(*sid));
    if (!more) break;
    first = FALSE;
    
    base = inUV->buffer[inUV->myDesc->ilocb]; /* Baseline */
    iAnt = (base / 256.0) + 0.001;
    if (inUV->myDesc->ilocsu>=0) 
      *sid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]+0.5);
    
    /* Get source elevation */
    elev = ObitAntennaListElev (AList, iAnt, time, 
				SList->SUlist[(*sid)-1]);
    
    /* Accumulate */
    TObs[0] = inUV->buffer[inUV->myDesc->nrparm+0];
    TObs[1] = inUV->buffer[inUV->myDesc->nrparm+3];
    TObs[2] = inUV->buffer[inUV->myDesc->nrparm+6];
    TObs[3] = inUV->buffer[inUV->myDesc->nrparm+9];
    weight  = inUV->buffer[inUV->myDesc->nrparm+2];
    /* Select good data */
    if ((TObs[0]<100.0) || (TObs[0]>300.0) ||
	(TObs[1]<20.0)  || (TObs[1]>200.0) ||
	(TObs[2]<10.0)  || (TObs[2]>200.0) ||
	(TObs[3]<10.0)  || (TObs[3]>200.0)) weight = 0.0;
    /* bad ALMA position? */
    if (elev<0.0) weight = 0.0;
    if (weight>0) {
      count[iAnt-1]++;
      sumWt[iAnt-1]        += weight;
      sumWtT[(iAnt-1)*4+0] += weight*TObs[0];
      sumWtT[(iAnt-1)*4+1] += weight*TObs[1];
      sumWtT[(iAnt-1)*4+2] += weight*TObs[2];
      sumWtT[(iAnt-1)*4+3] += weight*TObs[3];
      sumEl[iAnt-1]        += weight*elev;
    }

    /* Read next WVR record */
    if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else retCode =  ObitUVRead (inUV, inUV->buffer, err);
    /* Done? */
    done = (retCode==OBIT_IO_EOF) || (inUV->myDesc->firstVis>=inUV->myDesc->nvis);
    if (done) break;
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, done);
  } /* end loop over data */

  /* Package for threads */
  tr[0] = (ofloat)begTime;
  tr[1] = (ofloat)endTime;
  for (iAnt=0; iAnt<numAnt; ) {
    iTh = 0;
    while (iTh<nTh) {
      if (iAnt>=numAnt) break;
      /* Any data? */
      if ((count[iAnt]>0) && (sumWt[iAnt]>0.0)) {
	/* Data to libALMAWVR struct */
	ThreadArgs[iTh]->wvrIn.el      = (double)(sumEl[iAnt] / sumWt[iAnt]);
	ThreadArgs[iTh]->wvrIn.TObs[0] = (double)(sumWtT[iAnt*4+0] / sumWt[iAnt]);
	ThreadArgs[iTh]->wvrIn.TObs[1] = (double)(sumWtT[iAnt*4+1] / sumWt[iAnt]);
	ThreadArgs[iTh]->wvrIn.TObs[2] = (double)(sumWtT[iAnt*4+2] / sumWt[iAnt]);
	ThreadArgs[iTh]->wvrIn.TObs[3] = (double)(sumWtT[iAnt*4+3] / sumWt[iAnt]);
	time = 0.5*(begTime+endTime);
	ThreadArgs[iTh]->wvrIn.time    = (double)((time)*86400.0);  /* Not used */
	ThreadArgs[iTh]->wvrIn.state   = 0;                           /* Not used */
	ThreadArgs[iTh]->wvrIn.antno   = (size_t)(iAnt+1);            /* Not used */
	ThreadArgs[iTh]->wvrIn.source  = (size_t)(*sid);              /* Not used */
	iTh++;
      }
      iAnt++;    
    } /* end loop filling arguments */

    /* any good data? */
    if (iTh<=0) continue;

    /* If only one, indicate to thread */
    if (iTh==1) ThreadArgs[0]->ithread = -1;
    else        ThreadArgs[0]->ithread = 1;

    /* DEBUG
    fprintf (stderr, "nTh  %d ants %d %d\n",iTh, 
	     ThreadArgs[0]->wvrIn.antno, ThreadArgs[1]->wvrIn.antno); */
    
    /* Do operation */
    OK = ObitThreadIterator (inUV->thread, iTh, 
			     (ObitThreadFunc)ThreadWVR,
			     (gpointer**)ThreadArgs);
    
    /* Check for problems */
    Obit_retval_if_fail((OK), err, done,
			"%s: WVR Solution failed", routine);
    
    for (jTh =0; jTh<iTh; jTh++) {
      jAnt = (olong)(ThreadArgs[jTh]->wvrIn.antno - 1);
      if (ThreadArgs[jTh]->wvrOut.c_err>0.0) 
	wt = ThreadArgs[jTh]->wvrOut.c/ThreadArgs[jTh]->wvrOut.c_err;
      else
	wt = 0.0;
      dTdL[0] = (ofloat)ThreadArgs[jTh]->wvrOut.dTdL[0];
      dTdL[1] = (ofloat)ThreadArgs[jTh]->wvrOut.dTdL[1];
      dTdL[2] = (ofloat)ThreadArgs[jTh]->wvrOut.dTdL[2];
      dTdL[3] = (ofloat)ThreadArgs[jTh]->wvrOut.dTdL[3];

      /* Save entry if a valid number */
      if (!isnan(ThreadArgs[jTh]->wvrOut.c)) {
	ObitWVRCoefAdd (wvrcoef, tr, *sid, jAnt+1, dTdL, wt);

	/* Diagnostic */
	if (err->prtLv>=3) {
	  T2String ((ofloat)(ThreadArgs[jTh]->wvrIn.time/86400.0), tStr);
	  elev = ThreadArgs[jTh]->wvrIn.el * RAD2DG;
	  Obit_log_error(err, OBIT_InfoErr, "%s ant %3d sou %3d data %6.1f %6.1f %6.1f %6.1f el %5.1f path %8.6f (%8.6f)",
			 tStr, jAnt+1, *sid,
			 ThreadArgs[jTh]->wvrIn.TObs[0], ThreadArgs[jTh]->wvrIn.TObs[1],
			 ThreadArgs[jTh]->wvrIn.TObs[2], ThreadArgs[jTh]->wvrIn.TObs[3],
			 elev,
			 ThreadArgs[jTh]->wvrOut.c, ThreadArgs[jTh]->wvrOut.c_err);
	}
      } else {  /* Solution went bad */
	  T2String ((ofloat)(ThreadArgs[jTh]->wvrIn.time/86400.0), tStr);
	  Obit_log_error(err, OBIT_InfoErr, "%s ant %3d sou %3d Failed", tStr, jAnt+1, *sid);
      }
    } /* end loop retrieving values */
  } /* end loop over antenna */ 

#else  /* No libALMAWVR - stubb */
  gchar *routine = "AvgWVR1";
  Obit_log_error(err, OBIT_Error, 
		 "%s: libALMAWVR not available - cannot do fit", 
		     routine);
#endif /* HAVE_WVR */

  return done;
} /* end  AvgWVR1 */

/**
 *  Convert selected WVR data to delay and average for solInt
 * \param   inUV     ObitUV with WVR data
 * \param   solInt   Desired solution interval in days
 * \param   numAnt   Number of antennas
 * \param   AList    Antenna list
 * \param   SList    Source list
 * \param   timec    [out] Center time (day)
 * \param   timei    [out] Time interval (day)
 * \param   sid      [out] source ID
 * \param   fqid     [out] freq ID
 * \param   refAnt   [out] reference antenna
 * \param   gotAnt   [out] Flag if have data for antennas
 * \param   delay    [out] Delay (sec) per antenna, fblank if invalid
 * \param   wt       [out] Weights for delays
 * \param   wvrcoef  [out] Structure saving fitted sensitivities
 * \param   count    Work array
 * \param   sumWt    Work array
 * \param   sumWtT   Work array
 * \param   sumEl      Work array
 * \param   err        Obit Error stack
 * \return TRUE if all data processed, could be some OK in output 
 */
gboolean AvgWVR2 (ObitUV* inUV, ofloat solInt, olong numAnt, 
		  ObitAntennaList *AList, ObitSourceList *SList,
		  odouble *timec, ofloat *timei, olong *sid, olong *fqid, 
		  olong *refAnt, gboolean *gotAnt, ofloat *delay, ofloat *wt, 
		  ObitWVRCoef *wvrcoef,
		  olong *count, odouble *sumWt, odouble *sumWtT, odouble *sumEl,
		  ObitErr* err)
{
  gboolean done=TRUE;
  olong iAnt, tsid;
  odouble endTime, begTime;
  ofloat elev, time, weight, base, refDelay, TObs[4];
  ofloat fblank = ObitMagicF();
  gboolean more, doCalSelect, first;
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *routine = "AvgWVR2";
  
  /* error checks */
  if (err->error) return done;
  
  /* Init */
  *fqid  = -1;
  *timec = -1.0e20;
  *timei = 0.0;
  begTime =  1.0e20;
  endTime = -1.0e20;
  for (iAnt=0; iAnt<numAnt; iAnt++) {
    gotAnt[iAnt]     = FALSE;
    wt[iAnt]         = 0.0;
    count[iAnt]      = 0;
    sumWt[iAnt]      = 0.0;
    sumEl[iAnt]      = 0.0;
    sumWtT[iAnt*4+0] = 0.0;
    sumWtT[iAnt*4+1] = 0.0;
    sumWtT[iAnt*4+2] = 0.0;
    sumWtT[iAnt*4+3] = 0.0;
  }
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, 
		      &doCalSelect);
  
  if (inUV->myDesc->ilocsu>=0) 
    *sid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]+0.5);
  more = first = TRUE;
  while (more) {
    /* Info */
    time = inUV->buffer[inUV->myDesc->iloct];
    begTime = MIN (begTime, (odouble)time);
    endTime = MAX (endTime, (odouble)time);
    if (inUV->myDesc->ilocsu>=0) 
      tsid = (olong)(inUV->buffer[inUV->myDesc->ilocsu]+0.5);
    else tsid = *sid;
    /* finished? */
    more    = (endTime<(begTime+solInt)) && (tsid==(*sid));
    if (!more) break;
    first = FALSE;
    
    base = inUV->buffer[inUV->myDesc->ilocb]; /* Baseline */
    iAnt = (base / 256.0) + 0.001;
    if (inUV->myDesc->ilocsu>=0) *sid = tsid;
    
    /* Get source elevation */
    elev = ObitAntennaListElev (AList, iAnt, time, 
				SList->SUlist[(*sid)-1]);
    
    /* Accumulate */
    TObs[0] = inUV->buffer[inUV->myDesc->nrparm+0];
    TObs[1] = inUV->buffer[inUV->myDesc->nrparm+3];
    TObs[2] = inUV->buffer[inUV->myDesc->nrparm+6];
    TObs[3] = inUV->buffer[inUV->myDesc->nrparm+9];
    weight  = inUV->buffer[inUV->myDesc->nrparm+2];
    /* bad ALMA position? */
    if (elev<0.0) weight = 0.0;
    /* Select good data */
    if ((TObs[0]<100.0) || (TObs[0]>300.0) ||
	(TObs[1]<20.0)  || (TObs[1]>200.0) ||
	(TObs[2]<10.0)  || (TObs[2]>200.0) ||
	(TObs[3]<10.0)  || (TObs[3]>200.0)) weight = 0.0;
    if (weight>0) {
      count[iAnt-1]++;
      sumWt[iAnt-1]        += weight;
      sumWtT[(iAnt-1)*4+0] += weight*TObs[0];
      sumWtT[(iAnt-1)*4+1] += weight*TObs[1];
      sumWtT[(iAnt-1)*4+2] += weight*TObs[2];
      sumWtT[(iAnt-1)*4+3] += weight*TObs[3];
      sumEl[iAnt-1]        += weight*elev;
    }

    /* Read next WVR record */
    if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else retCode =  ObitUVRead (inUV, inUV->buffer, err);
    /* Done? */
    done = (retCode==OBIT_IO_EOF) || (inUV->myDesc->firstVis>=inUV->myDesc->nvis);
    if (done) break;
    if ((retCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inUV->name, done);
  } /* end loop over data */

  /* Calculate corrections */
  *timec = 0.5 * (begTime + endTime);
  *timei = (ofloat)(endTime - begTime);
  for (iAnt=0; iAnt<numAnt; iAnt++) {
    if ((count[iAnt]>0) && (sumWt[iAnt]>0.0)) {
      wt[iAnt] = sumWt[iAnt] / count[iAnt];
      gotAnt[iAnt] = wt[iAnt]>0.0;
      TObs[0] = (ofloat)(sumWtT[iAnt*4+0] / sumWt[iAnt]);
      TObs[1] = (ofloat)(sumWtT[iAnt*4+1] / sumWt[iAnt]);
      TObs[2] = (ofloat)(sumWtT[iAnt*4+2] / sumWt[iAnt]);
      TObs[3] = (ofloat)(sumWtT[iAnt*4+3] / sumWt[iAnt]);
      time = (ofloat)(*timec);
      /* Get delay in mm as correction */
      /* DEBUG delay[iAnt] = -ObitWVRCoefCal (wvrcoef, time, -1, iAnt+1, TObs);*/
      delay[iAnt] = ObitWVRCoefCal (wvrcoef, time, -1, iAnt+1, TObs);
      /* to Seconds */
      delay[iAnt] *= 1.0e-3 * CI;
    } else {
      delay[iAnt] = fblank;
      wt[iAnt]    = 0.0;
    }
  } /* end loop saving results */
  
  /* rereference phases - first make sure reference antenna valid */
  if (!gotAnt[(*refAnt)-1]) {
    for (iAnt=0; iAnt<numAnt; iAnt++) {
      if (gotAnt[iAnt]) {*refAnt = iAnt+1; break;}
    }
  }
  refDelay = delay[(*refAnt)-1];
  /* refDelay = 0.0;  DEBUG */
  for (iAnt=0; iAnt<numAnt; iAnt++) {
    if (gotAnt[iAnt]) {
      /* Refer to reference antenna */
      if ((wt[iAnt]>0) && (delay[iAnt]!=fblank) && ((*refAnt)>0))
	 delay[iAnt] -= refDelay;
    }
  } /* end loop referencing */
 return done;
  } /* end  AvgWVR2 */

/**
 * Make arguments for a Threaded ThreadWVRFunc?
 * \param thread     ObitThread object to be used
 * \param ThreadArgs[out] Created array of WVRFuncArg, 
 *                   delete with KillWVRFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeWVRFuncArgs (ObitThread *thread, 
			      WVRFuncArg ***ThreadArgs)

{
  olong nThreads=1;

#ifdef HAVE_WVR  /* Only if libALMAWVR available */
  olong i;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(WVRFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(WVRFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (*ThreadArgs)[i]->ithread  = i;
    (*ThreadArgs)[i]->wvrInArr[0]  = &(*ThreadArgs)[i]->wvrIn;
    (*ThreadArgs)[i]->wvrOutArr[0] = &(*ThreadArgs)[i]->wvrOut;
  }

#endif /* HAVE_WVR */
  return nThreads;
} /*  end MakeInterpImageArgs */

/**
 * Delete arguments for ThreadWVRFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of WVRFuncArg
 */
static void KillWVRFuncArgs (olong nargs, WVRFuncArg **ThreadArgs)
{
#ifdef HAVE_WVR  /* Only if libALMAWVR available */
   olong i;

 if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
#endif /* HAVE_WVR */
} /*  end KillWVRFuncArgs */

/**
 * Thread calculate WVR coefficients
 * Callable as thread
 * \param arg Pointer to WVRFuncArg argument with elements:
 * \li wvrInArr   Input data array.
 * \li wvrOutArr  Input data array.
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadWVR (gpointer arg)
{
#ifdef HAVE_WVR  /* Only if libALMAWVR available */
  /* Get arguments from structure */
  WVRFuncArg *largs           = (WVRFuncArg*)arg;

  /* Calculate water */
  almaabs_ret (largs->wvrInArr[0], 1, largs->wvrOutArr[0]);

  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
#endif /* HAVE_WVR */
  return NULL;
  
} /*  end ThreadWVR */

/**
 * Convert a time as time in days to a printable string
 * \param    time  Beginning time, end time in days
 * \msgBuff  Human readable string as "dd/hh:mm:ss.s"
 *           must be allocated at least 13 characters
 */
static void T2String (ofloat time, gchar *msgBuf)
{
  ofloat rtemp, rt1;
  olong   id1, it1, it2;

  id1 = time;
  rtemp = 24.0 * (time - id1);
  it1 = rtemp;
  it1 = MIN (23, it1);
  rtemp = (rtemp - it1)*60.0;
  it2 = rtemp;
  it2 = MIN (59, it2);
  rt1 = (rtemp - it2)*60.0;
  g_snprintf (msgBuf, 30, "%2.2d/%2.2d:%2.2d:%5.2f",
	      id1, it1, it2, rt1);
} /* end of routine T2String   */ 
