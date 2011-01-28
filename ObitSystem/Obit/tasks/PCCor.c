/* $Id$  */
/* Convert Pulse Cal (PC) table into an SN table                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUVUtil.h"
#include "ObitUVGSolveWB.h"
#include "ObitTableSN.h"
#include "ObitTablePC.h"
#include "ObitTableFQUtil.h"
#include "ObitPCal.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* PCCorIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void PCCorOut (ObitInfoList* outList, ObitErr *err);
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
/* Resolve single band delay ambiguities */
void GetSBAmbig(ObitInfoList *myInput, ObitUV *inData, ObitErr *err);
/* Convert PC Table to SN table */
void PC2SN (ObitInfoList *myInput, ObitUV *inData, ObitErr *err);
/* Add dummy entries for missing antennas from the PC table */
void DummySN(ObitInfoList *myInput, ObitUV *inData, ObitErr *err);
/* Write history */
void PCCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Average Delays/phases in SN Table */
ofloat** AverageDelayPhase (ObitTableSN *SNTable, ofloat *avgTime, 
			    ObitErr* err);
/* Get Frequency info */
void GetFQInfo (ObitInfoList *myInput, ObitUV *inData, ObitErr *err);
/* Fit a delay to a set of phases */
ofloat PCFitDelay (olong n, odouble *Freq, ofloat *Phase, ofloat *residRMS);

/* Program globals */
gchar *pgmName = "PCCor";       /* Program name */
gchar *infile  = "PCCor.in" ;   /* File with program inputs */
gchar *outfile = "PCCor.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
ObitTablePC  *PCin = NULL;     /* Input PC Table */
ObitTableSN  *SNOut = NULL;    /* Output SN Table */
oint  numAnt=1;                /* Number of Antennas */
ofloat **SBDelay = NULL;       /* List of Singleband delays per antenna, then per IF/poln */
ofloat **MBDelay = NULL;       /* List of Multiband delays per antenna, then per poln */
odouble *CableCal = NULL;      /* List of cable cal delays per antenna in cal reference time */
gboolean *gotAnt = NULL;       /* List of flags showing if an antenna 
				  has a valid solution */
ofloat *SBFact = NULL;         /* Sideband factors  per IF */
ofloat *chBandw=NULL;          /* Channel bandwidth per IF */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task to convert an AIPS PC table to an AIPS SN table            */
/*----------------------------------------------------------------------- */
{
  oint         i, ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;
  
  /* Startup - parse command line */
  err = newObitErr();
  myInput = PCCorIn (argc, argv, err);
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
  
  /* Get frequency info to global */
  GetFQInfo (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Resolve single band delay ambiguities */
  GetSBAmbig(myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Convert PC Table to SN table */
  PC2SN (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Add dummy entries for missing antennas from the PC table */
  DummySN (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Write history */
  PCCorHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  if (gotAnt) g_free(gotAnt); gotAnt = NULL;
  if (SBFact) g_free(SBFact); SBFact = NULL;
  if (CableCal) g_free(CableCal);  CableCal = NULL;
  if (SBDelay) {
    for (i=0; i<numAnt; i++) if (SBDelay[i]) g_free(SBDelay[i]);
     g_free(SBDelay); SBDelay = NULL;
  }
  
  if (MBDelay) {
    for (i=0; i<numAnt; i++) if (MBDelay[i]) g_free(MBDelay[i]);
     g_free(MBDelay); MBDelay = NULL;
  }
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* PCCorIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "PCCorIn";

  /* error checks */
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
} /* end PCCorIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: PCCor -input file -output ofile [args]\n");
    fprintf(stderr, "PCCor Converts AIPS PC to AIPS SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def PCCor.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def PCCor.out\n");
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
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  gboolean btemp;
  ofloat farray[2];
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

  /* error checks */
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
  strTemp = "PCCor.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "PCCorName";
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

  /*  Apply calibration/selection?, def=False */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doCalSelect", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  >0 => apply gain calibration, 2=> cal. wt, def=no cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doCalib", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Gain table (CL/SN) table to apply, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "gainUse", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  If >0.5 apply bandpass cal, def = no BP cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doBand", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Bandpass table version, 0=highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "BPVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "flagVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Apply polarization calibration?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "doPol", OBIT_bool, dim, &btemp, err);
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
  oint         doCalib;
  gchar        calSour[17] = "                ";
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doPCCor",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Copy "calSour" to "Sources" */
  ObitInfoListGetTest(myInput, "calSour",  &type, dim, calSour);
  dim[1] = strlen(calSour); dim[1] =  dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (myInput, "Sources", OBIT_string, dim, calSour);

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
  olong         Aseq, disk, cno, nvis, nThreads;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "timeRange", "BIF", "EIF", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  inData = newObitUV("input UV data");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
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
    nvis = 1000;
    nThreads = 1;
    ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
    nvis *= nThreads;
    ObitUVSetAIPS (inData, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    ObitTrimTrail(inFile);  /* remove trailing blanks */
  
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* define object */
    nvis = 1000;
    ObitUVSetFITS (inData, nvis, disk, inFile,  err); 
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }

   /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

 /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for PCCor                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void PCCorHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "calSour", "timeRange",  "subA", "BIF", "EIF", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "PCVer", "solnVer", "doCabCor", "doZeroMB",
    NULL};
  gchar *routine = "PCCorHistory";

  /* error checks */
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
 
} /* end PCCorHistory  */

/*----------------------------------------------------------------------- */
/* Resolve delay ambiguities                                              */
/* Follows general method in AIPSish PCCOR.FOR PCCAL                      */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/* Also global MBDelay = antenna list of MB delays per poln               */
/* Also global SBDelay = antenna list of SB delays per IF/poln            */
/*----------------------------------------------------------------------- */
void GetSBAmbig(ObitInfoList *myInput, ObitUV *inData, ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, PCVer=1, suba=1, npol, nif, nant, nfreq, size, indx, iant;
  olong nt, itone, refAnt=0;
  oint numPol=0, numBand=0, numTones=0;
  const ObitUVGSolveClassInfo *solnClass;
  ObitPCal* PCal=NULL;
  ObitTableSN *SNTable=NULL;
  ObitUV  *scrData = NULL;
  ofloat ftemp, **avgDelayPhase=NULL, avgTime, fblank = ObitMagicF();
  ofloat *fdata=NULL, *PCal1=NULL, *PCal2=NULL, *PCalR1=NULL, *PCalR2=NULL;
  odouble *ftone=NULL;
  odouble CableCl, dPhase;
  odouble *Freq1=NULL, *Freq2=NULL, *FreqR1=NULL, *FreqR2=NULL;
  odouble *refDelay=NULL, *IFFreq, refFreq;
  ofloat faz, del, f, ft, residRMS; 
  ObitUVGSolve *solver=NULL;
  gchar        *solverParms[] = {  /* Calibration parameters */
    "solInt", "solnVer", "solType", "solMode", "avgPol", "avgIF", "doMGM", "elevMGM",
    "refAnt", "refAnts", "doTwo", "ampScalar", "minSNR",  "minNo", "prtLv",
    NULL};
  gchar *routine = "GetSBAmbig";
  
  /* Reference antenna */
  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);
  if (refAnt<=0) {
    refAnt = 1;
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    ObitInfoListAlwaysPut (myInput, "refAnt", OBIT_long, dim, &refAnt);
  }
  
  /* Create input PC table object */
  ObitInfoListGetTest(myInput, "PCVer", &type, dim, &PCVer);
  ObitInfoListGetTest(myInput, "subA",  &type, dim, &suba);
  suba = MAX (1, suba);
  PCin = newObitTablePCValue (inData->name, (ObitData*)inData, &PCVer, 
			      OBIT_IO_ReadOnly, numPol, numBand, numTones, err);
  /* Open input table */
  ObitTablePCOpen (PCin, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  numPol   = PCin->numPol;
  numBand  = PCin->numBand;
  numTones = PCin->numTones;
  ObitTablePCClose (PCin, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Check compatability with UV data */
  npol  = MIN (2, inData->myDesc->inaxes[inData->myDesc->jlocs]);
  nif   = inData->myDesc->inaxes[inData->myDesc->jlocif];
  nfreq = inData->myDesc->inaxes[inData->myDesc->jlocf];
  nant  = inData->myDesc->numAnt[suba-1];/* actually highest antenna number */
  Obit_return_if_fail (((numPol==npol) && (numBand==nif)), err, 
		       "%s PC incomptatble with UV: pol %d != %d IF %d != %d", 
		       routine, numPol, npol, numBand, nif);
  /* Multiband delays and phases */
  numAnt = nant;
  MBDelay = g_malloc0(nant*sizeof(ofloat*));
  for (i=0; i<nant; i++) MBDelay[i] = g_malloc0(npol*sizeof(ofloat));
  
  /* Singleband delays */
  SBDelay = g_malloc0(nant*sizeof(ofloat*));
  for (i=0; i<nant; i++) SBDelay[i] = g_malloc0(nif*npol*sizeof(ofloat));
  
  /* Cable cal delays */
  CableCal = g_malloc0(nant*sizeof(odouble));
  
  /* Frequency sideband factors */
  SBFact = g_malloc0(nif*sizeof(ofloat));
  
  /* If there is only a single tone per band per pol then we're done */
  if (numTones<=1) return;
  
  /* Copy/select/calibrate to scratch file */
  scrData = newObitUVScratch (inData,err);
  scrData = ObitUVCopy (inData, scrData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Index scrData */
  dim[0] = dim[1] = 1;
  ftemp = 15.0;  /* Max scan time 15 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxScan", OBIT_float, dim, &ftemp);
  ftemp = 1.0; /* Max allowable gap 1 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxGap", OBIT_float, dim, &ftemp);
  ObitUVUtilIndex (scrData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create Calibration solution object */
  solver = (ObitUVGSolve*)ObitUVGSolveWBCreate("Gain solver");

  /* Copy calibration control to solver */
  ObitInfoListCopyList (myInput, solver->info, solverParms);

  /* Do solution */
  solnClass = (ObitUVGSolveClassInfo*)solver->ClassInfo;
  SNTable   = solnClass->ObitUVGSolveCal (solver, scrData, scrData, inData->mySel, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Average Delays/phases, note: solutions expressed as corrections */
  avgDelayPhase = AverageDelayPhase (SNTable, &avgTime, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Diagnostics */
  if (err->prtLv>=2) {
    for (iant=0; iant<nant; iant++) {
      for (i=0; i<nif; i++) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Cal Data Ant=%d, IF=%d, Poln 1 delay=%8.2f nsec phase=%8.2f deg,",
		       iant+1, i+1, avgDelayPhase[iant][i*2]*1.0e9, 
		       avgDelayPhase[iant][i*2+1]*57.296);
	if (npol>1) 
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Cal Data Ant=%d, IF=%d, Poln 2 delay=%8.2f nsec phase=%8.2f deg,",
			 iant+1, i+1, avgDelayPhase[iant][(i+nif)*2]*1.0e9, 
			 avgDelayPhase[iant][(i+nif)*2+1]*57.296);
      }
    }
  } /* end diagnostics */

  /* Get PCal phases at avgTime */
  PCal = ObitPCalCreate ("PCal interp", PCin, inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Work arrays */
  fdata = g_malloc0((numTones*nif+5)*sizeof(ofloat)); 
  ftone = g_malloc0((numTones*nif+5)*sizeof(odouble));
  size = nif * numTones;
  Freq1  = g_malloc0(size*sizeof(odouble));
  PCal1  = g_malloc0(size*sizeof(ofloat));
  FreqR1 = g_malloc0(size*sizeof(odouble));
  PCalR1 = g_malloc0(size*sizeof(ofloat));
  if (numPol>1) {
    Freq2  = g_malloc0(size*sizeof(odouble));
    PCal2  = g_malloc0(size*sizeof(ofloat));
    FreqR2 = g_malloc0(size*sizeof(odouble));
    PCalR2 = g_malloc0(size*sizeof(ofloat));
  }

  /* Reference antenna delay, entry per IF/poln */
  refDelay = g_malloc0(npol*nif*sizeof(odouble));

  /* Reference Frequency of first IF */
  IFFreq  = inData->myDesc->freqIF;
  refFreq = IFFreq[0];
  
  /* Get Antenna SB delays from PC - fblank if missing or one tone */
  for (iant=0; iant<nant; iant++) {
    /* Interpolate values */
    ObitPCalReport (PCal, avgTime, iant+1, suba, &CableCl, 
		    Freq1, PCal1, Freq2, PCal2, err);
    CableCal[iant] = CableCl;    /* Save cable cal at avgTime */
    /* Loop over IFs */
    for (i=0; i<nif; i++) {
      indx = i*numTones;
      if ((PCal1[indx]!=fblank) && (PCal1[indx+numTones-1]!=fblank) &&
	  (Freq1[indx]!=0.0) && (Freq1[indx+numTones-1]!=0.0) &&
	  (numTones>1)) {
	if (Freq1[indx] > Freq1[indx+numTones-1]) SBFact[i] = -1.0;
	else                                      SBFact[i] =  1.0;
	dPhase = SBFact[i]*(PCal1[indx] - PCal1[indx+numTones-1])/(2*G_PI);
	/* Phase difference in range +/- half turn */
	dPhase -= (olong)dPhase;
	if (dPhase> 0.5) dPhase -= 1.0;
	if (dPhase<-0.5) dPhase += 1.0;
	SBDelay[iant][i] = dPhase / (Freq1[indx] - Freq1[indx+numTones-1]);
      } else SBDelay[iant][i] = fblank;  /* No data or only one tone */
      /* Diagnostics */
      if (err->prtLv>=2) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Ant=%d, IF=%d, Poln 1 fit=%f PC=%f nsec,",
		       iant+1, i+1, avgDelayPhase[iant][i*2]*1.0e9, 
		       SBDelay[iant][i]*1.0e9);
      }
      /* Second poln */
      if (npol>1) {
	if ((PCal2[indx]!=fblank) && (PCal2[indx+numTones-1]!=fblank) &&
	    (Freq2[indx]!=0.0) && (Freq2[indx+numTones-1]!=0.0) &&
	    (numTones>1)) {
	  dPhase = SBFact[i]*(PCal2[indx] - PCal2[indx+numTones-1])/(2*G_PI);
	  /* Phase difference in range +/- half turn */
	  dPhase -= (olong)dPhase;
	  if (dPhase> 0.5) dPhase -= 1.0;
	  if (dPhase<-0.5) dPhase += 1.0;
	  SBDelay[iant][i+nif] = dPhase / (Freq2[indx] - Freq2[indx+numTones-1]);
	} else SBDelay[iant][i+nif] = fblank;  /* No data or only one tone  */
	/* Diagnostics */
	if (err->prtLv>=2) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Ant=%d, IF=%d, Poln 2 fit=%f PC=%f nsec",
			 iant+1, i+1, avgDelayPhase[iant][(i+nif)*2]*1.0e9, 
			 SBDelay[iant][i+nif]*1.0e9);
	}
      }
    }
  } /* End antnna loop getting antenna PC delays */
  
  /* Any messages */
  ObitErrLog(err);


  /* Interpolate values for reference antenna */
  ObitPCalReport (PCal, avgTime, refAnt, suba, &CableCl, 
		  FreqR1, PCalR1, FreqR2, PCalR2, err);

  /* Loop over antennas getting MB delay, resolving ambiguities */
  for (iant=0; iant<nant; iant++) {
    if ((iant+1)==refAnt) continue;  /* Not for reference antenna */
    /* Interpolate values */
    ObitPCalReport (PCal, avgTime, iant+1, suba, &CableCl, 
		    Freq1, PCal1, Freq2, PCal2, err);
    /* Collect phases and tone frequencies to determine MBdelay */
    nt = 0;
    /* Poln 1 */
    /* Loop over IFs */
    for (i=0; i<nif; i++) {
      /* Fitted from vis data: */
      del = -avgDelayPhase[iant][i*2];   /* Expressed as corrections */
      faz = -avgDelayPhase[iant][i*2+1];
      /* Correct for cal data fitting to correct ambiguity */
      f = faz + 2*G_PI * del * ((Freq1[i*numTones]-IFFreq[i]));
      for (itone=0; itone<numTones; itone++) {
	indx = i*numTones + itone;
	if ((PCal1[indx]!=fblank) && (Freq1[indx]!=0.0)) {
	  ftone[nt] = Freq1[indx];
	  if (itone>0) ft = f + 2*G_PI * del * (Freq1[indx]-Freq1[i*numTones]);
	  else ft = f;
	  /* Subtract PC phase of tone wrt ref ant */
	  ft = ft - PCal1[indx] + PCalR1[indx];
	  fdata[nt] = ft;
	  nt++;
	}
      } /* end loop over tones */
	
      indx = i*numTones;
      /* save difference between SB  calibrator and PC delays, refer to reference antenna */
      if ((SBDelay[iant][i]!=fblank) && (SBDelay[refAnt-1][i]!=fblank)) {
	SBDelay[iant][i] -= -del + SBDelay[refAnt-1][i];
      } else {
	SBDelay[iant][i] = del;
      }
    } /* end loop over IFs */

    /* Fit MB delay */
    MBDelay[iant][0] = -PCFitDelay (nt, ftone, fdata, &residRMS);

    /* Diagnostics */
    if (err->prtLv>=2) {
      Obit_log_error(err, OBIT_InfoErr, 
		    "Ant=%d, Poln 1 MBfit=%8.2f nsec RMS residPC=%8.5f turns, phase %8.2f deg",
		     iant+1, MBDelay[iant][0]*1.0e9, residRMS, fdata[0]*57.296);
    }
  
    /* Second Poln */
    if (npol>1) {
      nt = 0;
      /* Loop over IFs */
      for (i=0; i<nif; i++) {
	/* Fitted from vis data: */
	del = -avgDelayPhase[iant][(i+nif)*2];
	faz = -avgDelayPhase[iant][(i+nif)*2+1];
	/* Correct for cal data fitting to correct ambiguity */
	f = faz + 2*G_PI * del * ((Freq2[i*numTones]-IFFreq[i]));
	for (itone=0; itone<numTones; itone++) {
	  indx = i*numTones + itone;
	  if ((PCal2[indx]!=fblank) && (Freq2[indx]!=0.0)) {
	    ftone[nt] = Freq2[indx];
	    if (itone>0) ft = f + 2*G_PI * del * (Freq2[indx]-Freq2[i*numTones]);
	    else ft = f;
	    /* Subtract PC phase of tone wrt ref ant */
	    ft = ft - PCal2[indx] + PCalR2[indx];
	    fdata[nt] = ft;
	    nt++;
	  }
	} /* end loop over tones */
	
	indx = i*numTones;
        /* save difference between SB calibrator and PC delays, refer to reference antenna */
	if ((SBDelay[iant][i+nif]!=fblank) && (SBDelay[refAnt-1][i+nif]!=fblank)) {
	  SBDelay[iant][i+nif] -= -del + SBDelay[refAnt-1][i+nif];
	} else {
	  SBDelay[iant][i+nif] = del;
	}
	
      } /* end loop over IFs */
      /* Fit MB delay */
      MBDelay[iant][1] = -PCFitDelay (nt, ftone, fdata, &residRMS);

      /* Diagnostics */
      if (err->prtLv>=2) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Ant=%d, Poln 2 MBfit=%8.2f nsec RMS residPC=%8.5f turns, phase %8.2f deg",
		       iant+1, MBDelay[iant][1]*1.0e9, residRMS, fdata[0]*57.296);
      }
    } /* end second poln */
  }  /* End loop over antennas getting MB delay, resolving ambiguities */

  /* Cleanup */
  PCal    = ObitPCalUnref(PCal);
  SNTable = ObitTableSNUnref(SNTable);
  if (refDelay) g_free(refDelay);
  if (Freq1)    g_free(Freq1);
  if (Freq2)    g_free(Freq2);
  if (FreqR1)   g_free(FreqR1);
  if (FreqR2)   g_free(FreqR2);
  if (PCal1)    g_free(PCal1);
  if (PCal2)    g_free(PCal2);
  if (PCalR1)   g_free(PCalR1);
  if (PCalR2)   g_free(PCalR2);
  if (fdata)    g_free(fdata);
  if (ftone)    g_free(ftone);
  if (avgDelayPhase) {
    for (i=0; i<nant; i++) if (avgDelayPhase[i]) g_free(avgDelayPhase[i]);
    g_free(avgDelayPhase);
  }
     
} /* end GetSBAmbig */

/*----------------------------------------------------------------------- */
/* Convert PC Table to SN table                                           */
/* Follows general method in AIPSish PCCOR.FOR PCTAG                      */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void PC2SN (ObitInfoList *myInput, ObitUV *inData, ObitErr *err)
{
  olong i, PCVer=1, suba=1, refAnt, iAnt, indx;
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *SNrow=NULL;
  ObitTablePC *PCin=NULL;
  ObitTablePCRow *PCrow=NULL;
  olong SNver, iSNRow, iPCRow, itemp, n2pi;
  oint numPol, numIF, numBand=0, numTones=0;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  odouble dPhase, delay, addDelay1, addDelay2, DFreq, *IFFreq;
  ofloat temp, phz, phase, fblank = ObitMagicF();
  odouble refFreq, CabCor;
  gboolean doCabCor=FALSE, doZeroMB=FALSE;
  gchar *tname;
  gchar *routine = "PC2SN";
 
  ObitInfoListGetTest(myInput, "doZeroMB", &type, dim, &doZeroMB);
  ObitInfoListGetTest(myInput, "doCabCor", &type, dim, &doCabCor);
  ObitInfoListGetTest(myInput, "refAnt",   &type, dim, &refAnt);
  if (refAnt<=0) {
    refAnt = 1;
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    ObitInfoListAlwaysPut (myInput, "refAnt", OBIT_long, dim, &refAnt);
  }

  /* Create input PC table object */
  ObitInfoListGetTest(myInput, "PCVer", &type, dim, &PCVer);
  ObitInfoListGetTest(myInput, "subA",  &type, dim, &suba);
  PCin = newObitTablePCValue (inData->name, (ObitData*)inData, &PCVer, 
			      OBIT_IO_ReadOnly, numPol, numBand, numTones, err);
  /* Open input table */
  retCode = ObitTablePCOpen (PCin, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  numPol   = PCin->numPol;
  numBand  = PCin->numBand;
  numTones = PCin->numTones;

  /* Create output - version requested? */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  /* 0=> make new */
  itemp = ObitTableListGetHigh (inData->tableList, "AIPS SN") + 1;
  if (SNver<=0) SNver = itemp;
  tname = g_strconcat ("SN Calibration for: ", inData->name, NULL);
  if (inData->myDesc->jlocs>=0)
    numPol = MIN (2, inData->myDesc->inaxes[inData->myDesc->jlocs]);
  else numPol = 1;
  if (inData->myDesc->jlocif>=0)
    numIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else numIF  = 1;
  outSoln = newObitTableSNValue(tname, (ObitData*)inData, &SNver, OBIT_IO_ReadWrite, 
				numPol, numIF, err);
  g_free (tname);
  /* Open output table */
  retCode = ObitTableSNOpen (outSoln, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Header stuff */
  outSoln->numAnt = numAnt;
  
  /* Save actual SN table used */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (myInput, "solnVer", OBIT_long, dim, &SNver);

  /* Create Rows */
  SNrow = newObitTableSNRow (outSoln);
  PCrow = newObitTablePCRow (PCin);
  
  /* Attach output row to output buffer */
  ObitTableSNSetRow (outSoln, SNrow, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Reference Frerquency of first IF */
  IFFreq  = inData->myDesc->freqIF;
  refFreq = IFFreq[0];
  
  /* Initialize solution row */
  SNrow->Time   = 0.0; 
  SNrow->TimeI  = 0.0; 
  SNrow->SourID = 0; 
  SNrow->antNo  = 0; 
  SNrow->SubA   = 0; 
  SNrow->FreqID = 0; 
  SNrow->IFR    = 0.0; 
  SNrow->NodeNo = 0; 
  SNrow->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    SNrow->Real1[i]   = 0.0; 
    SNrow->Imag1[i]   = 0.0; 
    SNrow->Delay1[i]  = 0.0; 
    SNrow->Rate1[i]   = 0.0; 
    SNrow->Weight1[i] = 1.0; 
    SNrow->RefAnt1[i] = refAnt; 
  }
  if (numPol>1) {
    SNrow->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      SNrow->Real2[i]   = 0.0; 
      SNrow->Imag2[i]   = 0.0; 
      SNrow->Delay2[i]  = 0.0; 
      SNrow->Rate2[i]   = 0.0; 
      SNrow->Weight2[i] = 1.0; 
      SNrow->RefAnt2[i] = refAnt; 
    }
  }
  
  /* Keep track of which antennas have values */
  gotAnt  = g_malloc0(numAnt*sizeof(gboolean));
  for (i=0; i<numAnt; i++) gotAnt [i]= FALSE;
  
  /* Loop over PC table converting */
  for (iPCRow=1; iPCRow<=PCin->myDesc->nrow; iPCRow++) {
    ObitTablePCReadRow (PCin, iPCRow, PCrow, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    if (PCrow->status==-1) continue;
    
    /* Set Solution */
    iAnt          = PCrow->antennaNo-1;
    SNrow->Time   = PCrow->Time;
    SNrow->TimeI  = PCrow->TimeI;
    SNrow->SourID = PCrow->SourID; 
    SNrow->antNo  = PCrow->antennaNo; 
    SNrow->SubA   = PCrow->Array; 
    SNrow->FreqID = PCrow->FreqID; 

    /* Cable corrections */
    if (doCabCor) CabCor = CableCal[iAnt] - PCrow->CableCal;
    else          CabCor = 0.0;

    /* Multiband delays */
    addDelay1 = MBDelay[iAnt][0] + CabCor;
    if (doZeroMB) SNrow->MBDelay1 = 0.0;
    else SNrow->MBDelay1 = addDelay1; 
    if (SNrow->antNo==refAnt) SNrow->MBDelay1 = 0.;  /* Reference antenna = 0 */
    if (numPol>1) {
      addDelay2 = MBDelay[iAnt][1] + CabCor;
      if (doZeroMB) SNrow->MBDelay2 = 0.0;
      else SNrow->MBDelay2 = addDelay2; 
      if (SNrow->antNo==refAnt) SNrow->MBDelay2 = 0.;  /* Reference antenna = 0 */
    }

    /* IF Loop Poln 1 */
    for (i=0; i<numIF; i++) {
      indx = i*numTones;
      /* Good? */
      if ((PCrow->PCReal1[indx]!=fblank) && (PCrow->PCImag1[indx]!=fblank)) {
	/* Single band delay? */
	if (numTones>1) {
	  dPhase = SBFact[i]*(atan2(PCrow->PCImag1[indx],PCrow->PCReal1[indx]) -
			      atan2(PCrow->PCImag1[indx+numTones-1],PCrow->PCReal1[indx+numTones-1])) /
	    (2*G_PI);
	  /* Phase difference in range +/- half turn */
	  dPhase -= (olong)dPhase;
	  if (dPhase> 0.5) dPhase -= 1.0;
	  if (dPhase<-0.5) dPhase += 1.0;
	  DFreq = (PCrow->PCFreq1[indx] - PCrow->PCFreq1[indx+numTones-1]);
	  delay = dPhase / DFreq;

	  /* Resolve ambiguity */
	  temp = (delay - SBDelay[iAnt][i]) * DFreq;
	  if (temp>0.0) n2pi = (olong)(temp + 0.5);
	  else  n2pi = (olong)(temp - 0.5);
	  n2pi = (olong)(temp);
	  delay -= n2pi / DFreq;

	  SNrow->Delay1[i]  = delay - SBDelay[iAnt][i] - SBDelay[refAnt-1][i];
	  if (SNrow->antNo==refAnt)  SNrow->Delay1[i] = 0.;  /* Reference antenna = 0 */

	} else { /* Only one tone */
	  SNrow->Delay1[i]  = - SBDelay[iAnt][i];
	}
	/* Use phase of first tone */
	phz   = SBFact[i]*atan2(PCrow->PCImag1[indx], PCrow->PCReal1[indx]);
	phase = phz - 2*G_PI * (delay * (PCrow->PCFreq1[indx]- IFFreq[i]) +
						   (addDelay1 * (PCrow->PCFreq1[indx] - refFreq)));
	/* DEBUG if (SNrow->antNo==refAnt) phase = 0.;  Reference antenna = 0 */
	phase = fmodf(phase, 2*G_PI);
	SNrow->Real1[i]   = cos(phase);
	SNrow->Imag1[i]   = -sin(phase);
	SNrow->Weight1[i] = 1.0; 
	SNrow->Rate1[i]   =  PCrow->PCRate1[i]; 
	/* Diagnostics */
	if (err->prtLv>=3) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "time=%f Ant=%d, IF=%d, Poln 1 delay=%8.2f nsec phase %8.2f deg",
			 PCrow->Time, iAnt+1, i+1, SNrow->Delay1[i]*1.0e9, phase*57.296);
	  ObitErrLog(err);
	}
      } else { /* bad */
	SNrow->Real1[i]   = fblank;
	SNrow->Imag1[i]   = fblank;
	SNrow->Weight1[i] = 0.0; 
      }
    }
    /* IF Loop Poln 2 */
    if (numPol>1) {
      for (i=0; i<numIF; i++) {
	/* Good? */
	indx = i*numTones;
	if ((PCrow->PCReal2[indx]!=fblank) && (PCrow->PCImag2[indx]!=fblank)) {
	  if (numTones>1) {
	    SNrow->Delay2[i]  = 0.0; 
	    dPhase = SBFact[i]*(atan2(PCrow->PCImag2[indx],PCrow->PCReal2[indx]) -
				atan2(PCrow->PCImag2[indx+numTones-1],PCrow->PCReal2[indx+numTones-1])) /
	      (2*G_PI);
	    /* Phase difference in range +/- half turn */
	    dPhase -= (olong)dPhase;
	    if (dPhase> 0.5) dPhase -= 1.0;
	    if (dPhase<-0.5) dPhase += 1.0;
	    DFreq = (PCrow->PCFreq2[indx] - PCrow->PCFreq2[indx+numTones-1]);
	    delay = dPhase / DFreq;

	    /* Resolve ambiguity */
	    temp = (delay - SBDelay[iAnt][i+numIF]) * DFreq;
	    if (temp>0.0) n2pi = (olong)(temp + 0.5);
	    else  n2pi = (olong)(temp - 0.5);
	    n2pi = (olong)(temp);
	    delay -= n2pi / DFreq;
	    
	    SNrow->Delay2[i]  = delay - SBDelay[iAnt][i+numIF] - SBDelay[refAnt-1][i+numIF];
	    if (SNrow->antNo==refAnt) SNrow->Delay2[i] = 0.;  /* Reference antenna = 0 */

	  } else {  /* Single tone */
	    SNrow->Delay2[i]  = - SBDelay[iAnt][i+numIF];
	  }
	  /* Use phase of first tone */
	  phz = SBFact[i]*atan2(PCrow->PCImag2[indx], PCrow->PCReal2[indx]);
	  phase = phz - 2*G_PI * (delay * (PCrow->PCFreq2[indx]- IFFreq[i]) +
							(addDelay2 * (PCrow->PCFreq2[indx] - refFreq)));
	  /* DEBUG if (SNrow->antNo==refAnt) phase = 0.;   Reference antenna = 0 */
	  phase = fmodf(phase, 2*G_PI);
	  SNrow->Real2[i]   = cos(phase);
	  SNrow->Imag2[i]   = -sin(phase);
	  SNrow->Weight2[i] = 1.0; 
	  SNrow->Rate2[i]   =  PCrow->PCRate2[i]; 
	  /* Diagnostics */
	  if (err->prtLv>=3) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "time=%f Ant=%d, IF=%d, Poln 2 delay=%8.2f nsec phase %8.2f deg",
			   PCrow->Time, iAnt+1, i+1, SNrow->Delay2[i]*1.0e9, phase*57.296);
	    ObitErrLog(err);
	  }
	} else { /* bad */
	  SNrow->Real2[i]   = fblank;
	  SNrow->Imag2[i]   = fblank;
	  SNrow->Weight2[i] = 0.0; 
	}
      } /* end IF loop */
    }
    /* If good mark antenna */
    gotAnt[PCrow->antennaNo-1] = TRUE;
    
    /* Write SN table */
    iSNRow = -1;
    retCode = ObitTableSNWriteRow (outSoln, iSNRow, SNrow, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end loop over table */
  
  /* Close */
  ObitTableSNClose (outSoln, err);
  ObitTablePCClose (PCin, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Cleanup */
  outSoln = ObitTableSNUnref(outSoln);
  SNrow   = ObitTableSNRowUnref(SNrow);
  PCin    = ObitTablePCUnref(PCin);
  PCrow   = ObitTablePCRowUnref(PCrow);

} /* end PC2SN */

/*----------------------------------------------------------------------- */
/* Add dummy entries for missing antennas from the PC table               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void DummySN(ObitInfoList *myInput, ObitUV *inData, ObitErr *err) 
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *SNrow=NULL;
  olong i, SNver, iSNRow, refAnt;
  oint numPol, numIF;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode;
  gboolean done=FALSE;
  gchar *tname;
  gchar *routine = "DummySN";
 
  if (err->error) return;

  /* Need to do anything? */
  for (i=0; i<numAnt; i++) done = done && gotAnt[i];
  if (done) return;

  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);
  ObitInfoListGetTest(myInput, "solnVer", &type, dim, &SNver);
  tname = g_strconcat ("SN Calibration for: ", inData->name, NULL);
  if (inData->myDesc->jlocs>=0)
    numPol = MIN (2, inData->myDesc->inaxes[inData->myDesc->jlocs]);
  else numPol = 1;
  if (inData->myDesc->jlocif>=0)
    numIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else numIF  = 1;
  outSoln = newObitTableSNValue(tname, (ObitData*)inData, &SNver, OBIT_IO_ReadWrite, 
				numPol, numIF, err);
  g_free (tname);
  /* Open output table */
  retCode = ObitTableSNOpen (outSoln, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Create Row */
  SNrow = newObitTableSNRow (outSoln);
  /* Attach output row to output buffer */
  ObitTableSNSetRow (outSoln, SNrow, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Initialize solution row */
  SNrow->Time   = 0.0; 
  SNrow->TimeI  = 0.0; 
  SNrow->SourID = -1; 
  SNrow->antNo  = 0; 
  SNrow->SubA   = 0; 
  SNrow->FreqID = 0; 
  SNrow->IFR    = 0.0; 
  SNrow->NodeNo = 0; 
  SNrow->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    SNrow->Real1[i]   = 1.0; 
    SNrow->Imag1[i]   = 0.0; 
    SNrow->Delay1[i]  = 0.0; 
    SNrow->Rate1[i]   = 0.0; 
    SNrow->Weight1[i] = 1.0; 
    SNrow->RefAnt1[i] = refAnt; 
  }
  if (numPol>1) {
    SNrow->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      SNrow->Real2[i]   = 1.0; 
      SNrow->Imag2[i]   = 0.0; 
      SNrow->Delay2[i]  = 0.0; 
      SNrow->Rate2[i]   = 0.0; 
      SNrow->Weight2[i] = 1.0; 
      SNrow->RefAnt2[i] = refAnt; 
    }
  }
  
  /* Loop over antennas */
  for (i=0; i<numAnt; i++) {
    if (gotAnt[i]) continue;  /* Need this? */
    SNrow->antNo = i+1;
    /* Write SN table */
    iSNRow = -1;
    retCode = ObitTableSNWriteRow (outSoln, iSNRow, SNrow, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end loop over antenna */

  /* Close */
  ObitTableSNClose (outSoln, err);

  /* Cleanup */
  outSoln = ObitTableSNUnref(outSoln);
  SNrow   = ObitTableSNRowUnref(SNrow);
} /* end DummySN */

/*----------------------------------------------------------------------- */
/* Average Delay/phase                                                    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*      nant      Number of antennas (actually max)                       */
/*      nif       Number of IFs                                           */
/*      npol      Number of Poln (parallel only)                          */
/*   Output:                                                              */
/*      avgTime Average time (day) of SN Table                            */
/*      err     Obit Error stack                                          */
/*   Returns averaged delay array,per antenna each delay+phase per IF/poln*/
/*      should be g_freeed when done, may be NULL on failure              */
/*----------------------------------------------------------------------- */
ofloat** AverageDelayPhase(ObitTableSN *SNTable, ofloat *avgTime, ObitErr* err)
{
  ofloat **out=NULL;
  ofloat **wt=NULL, **sumR=NULL, **sumI=NULL, sumTime=0.0;
  ofloat fblank = ObitMagicF();
  ObitTableSNRow *SNRow=NULL;
  olong iRow, nant, nif, npol, i, j, size, iant, cntTime=0;
  ObitIOCode retCode;
  gchar *routine="AverageDelay";

  if (err->error) return out;
    
  /* Open input table */
  retCode = ObitTableSNOpen (SNTable, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, SNTable->name, out);
  /* Row for I/O */
  SNRow  = newObitTableSNRow (SNTable);
  
  /* Create output, counts */
  nant  = SNTable->numAnt;
  nif   = SNTable->numIF;
  npol  = SNTable->numPol;
  out   = g_malloc0(nant*sizeof(ofloat*));
  wt    = g_malloc0(nant*sizeof(ofloat*));
  sumR  = g_malloc0(nant*sizeof(ofloat*));
  sumI  = g_malloc0(nant*sizeof(ofloat*));
  size  =  nif*npol;
  for (i=0; i<nant; i++) {
    out[i]  = g_malloc0(2*size*sizeof(ofloat));
    wt[i]   = g_malloc0(size*sizeof(ofloat));
    sumR[i] = g_malloc0(size*sizeof(ofloat));
    sumI[i] = g_malloc0(size*sizeof(ofloat));
    for (j=0; j<size; j++) {
      out[i][j]  = 0.0;
      wt[i][j]   = 0.0;
      sumR[i][j] = 0.0;
      sumI[i][j] = 0.0;
    }
  }
 
  /* Loop over table copying selected data */
  for (iRow=1; iRow<=SNTable->myDesc->nrow; iRow++) {
    retCode = ObitTableSNReadRow (SNTable, iRow, SNRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, SNTable->name, out);
    if (SNRow->status==-1) continue;
    /* Time */
    sumTime += SNRow->Time;
    cntTime++;
    iant = SNRow->antNo - 1;     /* zero rel antenna number */
    iant = MAX(0, MIN(iant, (nant-1)));
    /* Loop over IFs */
    for (i=0; i<nif; i++) {
      if ((SNRow->Weight1[i]>0.0) && (SNRow->Real1[i]!=fblank)) {
	out[iant][i*2] += SNRow->Delay1[i]*SNRow->Weight1[i];
	sumR[iant][i]  += SNRow->Real1[i]*SNRow->Weight1[i];
	sumI[iant][i]  += SNRow->Imag1[i]*SNRow->Weight1[i];
	wt[iant][i]    += SNRow->Weight1[i];
      }
      /* two poln? */
      if ((npol>1) && (SNRow->Weight2[i]>0.0) && (SNRow->Real2[i]!=fblank)) {
	out[iant][(i+nif)*2] += SNRow->Delay2[i]*SNRow->Weight2[i];
	sumR[iant][i+nif]    += SNRow->Real2[i]*SNRow->Weight2[i];
	sumI[iant][i+nif]    += SNRow->Imag2[i]*SNRow->Weight2[i];
	wt[iant][i+nif]      += SNRow->Weight2[i];
      }
    }
  } /* end loop over table */

  /* Close table */
  retCode = ObitTableSNClose (SNTable, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, SNTable->name, out);
  
  /* release row object */
  SNRow  = ObitTableSNRowUnref(SNRow);
  
  /* Average */
  if (cntTime>0) *avgTime = sumTime/cntTime;
  else  *avgTime = 0.0;
  for (iant=0; iant<nant; iant++) {
    for (i=0; i<nif; i++) {
      if (wt[iant][i]>0.0) {
	  out[iant][i*2]   /= wt[iant][i];
	  out[iant][i*2+1] = atan2(sumI[iant][i], sumR[iant][i]);
      }
      if ((npol>1) && (wt[iant][i+nif]>0.0)) {
	out[iant][(i+nif)*2] /= wt[iant][i+nif];
	out[iant][(i+nif)*2+1] = atan2(sumI[iant][i+nif], sumR[iant][i+nif]);
      }
    }
  }
  
  /* Cleanup */
  if (wt) {
    for (iant=0; iant<nant; iant++) 
      if (wt[iant]) g_free(wt[iant]);
    g_free(wt);
  }
  return out;
} /* end AverageDelayPhase */

/*----------------------------------------------------------------------- */
/* Get Frequency information                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV                                            */
/*   Output:                                                              */
/*      err     Obit Error stack                                          */
/* output in common chBandw                                               */
/*----------------------------------------------------------------------- */
void GetFQInfo (ObitInfoList *myInput, ObitUV *inData, ObitErr *err)
{
  ObitTableFQ *FQTab=NULL;
  olong ver;
  oint  numIF, *sideBand=NULL, fqid;
  odouble *freqOff=NULL;
  gchar *routine = "GetFQInfo";

  if (inData->myDesc->jlocif>=0)
    numIF  = inData->myDesc->inaxes[inData->myDesc->jlocif];
  else numIF  = 1;

  /* Create output 
     SBFact = g_malloc0(numIF*sizeof(ofloat));*/
  /* Default 
     for (i=0; i<numIF; i++) SBFact[i] = 1.0;*/

  ver = 1;
  FQTab = newObitTableFQValue("FQTab", (ObitData*)inData, &ver, OBIT_IO_ReadOnly, 
			      numIF, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  if (FQTab==NULL) return; /* Bummer */

  fqid = 1; 

  ObitTableFQGetInfo (FQTab, fqid, &numIF, &freqOff, &sideBand, &chBandw, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Save sideband info 
     for (i=0; i<numIF; i++) SBFact[i] *= sideBand[i]; */

  /* Cleanup */
  FQTab = ObitTableFQUnref(FQTab);
  if (sideBand) g_free(sideBand);
  if (freqOff)  g_free(freqOff);
  
} /* end GetFQInfo */

/**
 *  Fit a delay to a set of phases
 *  Use direct parameter search.
 * \param n     Number of samples
 * \param Freq  Frequencies (Hz)
 * \param Phase Phases (radians) possibly fblanked
 * \param residRMS [out] residual RMS
 * \return Delay in sec.
 */
ofloat PCFitDelay (olong n, odouble *Freq, ofloat *Phase, ofloat *residRMS) 
{
  ofloat Delay = 0.0;
  olong i, j, m, besti;
  ofloat trial=0.0, bestDelay, model, dif, delta, fblank = ObitMagicF();
  ofloat *tPhase=NULL;
  odouble sum, bestSum, *tFreq=NULL, minDelta;

  if (n<=0) return Delay;

  /* Count good values */
  *residRMS = 0.0;
  m = 0;
  for (j=0; j<n; j++) if (Phase[j]!=fblank) m++;
  if (m<=0) return Delay;
  
  /* Make temporary arrays */
  tPhase = g_malloc0(m*sizeof(ofloat));
  tFreq  = g_malloc0(m*sizeof(odouble));
  i = 0;
  for (j=0; j<n; j++) {
    if (Phase[j]!=fblank) {
      tPhase[i] = Phase[j] / (2*G_PI);  /* Phase in turns */
      tFreq[i]  = Freq[j];
      if (i>0) {  /* Relative to first */
	tPhase[i] -= tPhase[0];
	tFreq[i]  -= tFreq[0];
	/* No more than half turn */
	tPhase[j] -= (olong)tPhase[j];
	if (tPhase[j]> 0.5) tPhase[j] -= 1.0;
	if (tPhase[j]<-0.5) tPhase[j] += 1.0;
      }
      i++;
    }
  }

  /* Zero first frequency/phase */
  tFreq[0]  = 0.0;
  tPhase[0] = 0.0;

  /* get min spacing for ambiguity */
  minDelta = 1.0e20;
  for (j=1; j<m; j++) minDelta = MIN (minDelta, tFreq[j]-tFreq[j-1]);
  minDelta = MAX(1.0e6, minDelta);  /* In case - it cannot be less than this */

  /* First test loop over delay */
  delta =  0.5e-3/minDelta;   /* Closest ambiguity */
  bestSum = 1.0e20; bestDelay = 0.0;
  for (i=0; i<2001; i++) {
    trial = (i-1000)*delta;
    
    /* Get RMS residual for this trial delay */
    sum = 0.0;
    for (j=1; j<m; j++) {
      model = trial * tFreq[j];
      dif   = model - tPhase[j];
      dif  -= (olong)dif;
      /* No more than half turn */
      if (dif> 0.5) dif -= 1.0;
      if (dif<-0.5) dif += 1.0;
      sum += dif*dif;
    }

    /* Is this the best so far? */
    if (sum<bestSum) {
      bestSum   = sum;
      bestDelay = trial;
      besti     = i;
    }
  } /* end first delay loop */

  Delay = bestDelay;  /* Save best */

  /* Second, fine test loop over delay */
  delta =  0.5e-5/minDelta;  
  bestSum = 1.0e20; bestDelay = 0.0;
  for (i=0; i<201; i++) {
    trial = Delay + (i-100)*delta;

    /* Get RMS residual for this trial delay */
    sum = 0.0;
    for (j=1; j<m; j++) {
      model = trial * tFreq[j];
      dif   = model - tPhase[j];
      dif  -= (olong)dif;
      /* No more than half turn */
      if (dif> 0.5) dif -= 1.0;
      if (dif<-0.5) dif += 1.0;
      sum += dif*dif;
    }

    /* Is this the best so far? */
    if (sum<bestSum) {
      bestSum   = sum;
      bestDelay = trial;
      besti     = i;
    }
  } /* end second delay loop */

  Delay = bestDelay;  /* Save best */

  /* RMS Resid */
  *residRMS = sqrt (sum/m);

  /* Debug 
  for (j=0; j<n; j++) {
    tPhase[j] -= Delay* tFreq[j];
    if (tPhase[j]> 0.5) tPhase[j] -= 1.0;
    if (tPhase[j]<-0.5) tPhase[j] += 1.0;
    } */

  /* Cleanup */
  if (tPhase) g_free(tPhase);
  if (tFreq) g_free(tFreq);

  return Delay;
} /* end PCFitDelay */

