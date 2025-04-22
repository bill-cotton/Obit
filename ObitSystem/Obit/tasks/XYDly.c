/* $Id$  */
/* X-Y delay/phase calibration                                        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableSN.h"
#include "ObitTableCL.h"
#include "ObitSinCos.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* XYDlyIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void XYDlyOut (ObitInfoList* outList, ObitErr *err);
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
/* Do solutions */
void XYDelayCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		ObitErr* err);
/* Write history */
void XYDlyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Actual fitting */
ObitTableSN* ObitUVXYDelayCal (ObitUV *inUV, ObitUV *outUV, ObitErr *err);
/* Fit phase/delay to complex spectrum */
void XYFitDelay (olong numFreq, olong BChan, ofloat dFreq,
		 ofloat *xpol1, ofloat *xpol2, olong fitType, 
		 ofloat *delay, ofloat *phase, ofloat *snr, ObitErr *err);
/* test delay */
ofloat XYtestDelay (olong numFreq, ofloat *xpol1, ofloat *xpol2, ofloat delay);
/* Get average IF freq */
odouble* GetAvgFreq (olong numIF, olong numFreq, ofloat *xpol1, ofloat *xpol2, 
		     ObitUV *inUV, ObitErr *err);
/* Get source lookup table */
olong* GetLookupSU (ObitSourceList *SList, ObitInfoList *info);
/* Setup EVPA sin/cos arrays */
void SetupEVPA (olong numIF, olong numFreq, ObitUV *inUV, 
		ofloat EVPA, ofloat RM, ofloat *EVPAFS, ofloat *EVPAFC);
/* Calculate sign flips */
void GetFlip (olong numIF, olong numFreq, ObitUV *inUV, ofloat *EVPAFS, ofloat *EVPAFC, 
	      ofloat chisin, ofloat chicos, ofloat *flip);
void sincosf(float x, float *sin, float *cos);
/* Human readable time string */
void day2dhms(ofloat time, gchar *timeString);

/* Program globals */
gchar *pgmName = "XYDly";      /* Program name */
gchar *infile  = "XYDly.in" ;  /* File with program inputs */
gchar *outfile = "XYDly.out";  /* File to contain program outputs */
olong  pgmNumber;               /* Program number (like POPS no.) */
olong  AIPSuser;                /* AIPS user number number  */
olong  nAIPS=0;                 /* Number of AIPS directories */
gchar **AIPSdirs=NULL;          /* List of AIPS data directories */
olong  nFITS=0;                 /* Number of FITS directories */
gchar **FITSdirs=NULL;          /* List of FITS data directories */
ObitInfoList *myInput  = NULL;  /* Input parameter list */
ObitInfoList *myOutput = NULL;  /* Output parameter list */
olong  souNo=0;                 /* Single source number */
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    X-Y (linear feed) phase/delay calibration                           */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData =NULL, *scrData=NULL, *avgData1=NULL, *avgData2=NULL;
  ObitErr      *err=NULL;
  ofloat       timeAvg;
  olong        chAvg;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select/process data */
    "Sources", "timeRange", "EVPA", "RM", "SNSoln", "solInt", "refAnt", "fitType", 
    NULL};

  /* Startup - parse command line / input */
  err = newObitErr();
  myInput = XYDlyIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Copy/select/calibrate to scratch file */
  scrData = newObitUVScratch (inData,err);
  scrData = ObitUVCopy (inData, scrData, err);
  /* Get input parameters from myInput, copy to scrData */
  ObitInfoListCopyList (myInput, scrData->info, dataParms);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Average in Frequency */
  chAvg = 8; /* Default 8 */
  ObitInfoListGetTest (myInput, "chAvg", &type, dim, &chAvg); 
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "NumChAvg", OBIT_long, dim, &chAvg);
  avgData1 = ObitUVUtilAvgF (scrData, TRUE, avgData1, err);
  /* Get input parameters from myInput, copy to avgData1 */
  ObitInfoListCopyList (myInput, avgData1->info, dataParms);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* Average in time */
  timeAvg = 2.0;   /* Default 2 min */
  ObitInfoListGetTest (myInput, "chAvg", &type, dim, &chAvg); 
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "timeAvg", OBIT_float, dim, &timeAvg);
  avgData2 = ObitUVUtilAvgT (avgData1, TRUE, avgData2, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* Do solutions, convert to SN table */
  XYDelayCal(myInput, avgData2, inData, err);
  if (err->error) {ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Write history */
  XYDlyHistory (myInput, inData, err); 
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  scrData   = ObitUnref(scrData);
  avgData1  = ObitUnref(avgData1);
  avgData2  = ObitUnref(avgData2);
 
  /* Shutdown Obit */
 exit: 
  if (myOutput) ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* XYDlyIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "XYDlyIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

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
} /* end XYDlyIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: XYDly -input file -output ofile [args]\n");
    fprintf(stderr, "XYDly Obit task XY phase/delay for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def XYDly.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def XYDly.out\n");
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
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*     Alpha      Flt(1)    default spectral index (0)                    */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat ftemp, farray[3];
  gboolean btemp;
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
  strTemp = "XYDly.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "XYDlyName";
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

  /*  Apply calibration/selection?, def=true */
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

  /* default Spectral index */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Alpha", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* BChan, EChan */
  dim[0] = dim[1] = dim[2] = 1;
  itemp = 1;
  ObitInfoListAlwaysPut (out, "BChan", OBIT_long, dim, &itemp);
  itemp = 0;
  ObitInfoListAlwaysPut (out, "EChan", OBIT_long, dim, &itemp);
 
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
  /*ObitInfoType type;*/
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

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
  olong        nvis, doCalib;
  gboolean     doCalSelect;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select/process data */
    "Sources", "timeRange", "FreqID", "BChan", "EChan",   "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "doPol", "PDVer", "keepLin", "UVR_Full", "WtUV", "SNSoln", "solInt",
    "timeAvg", "chAvg", "EVPA", "RM", "refAnt", "fitType", 
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
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
 
  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

   /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

 /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for XYDly                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void XYDlyHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA", "UVR_Full", "WtUV", 
    "selBand", "selFreq", "FreqID", "BChan", "EChan", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ", "BPVer", "flagVer", 
    "doPol", "PDVer", "keepLin", "Antennas",  "fitType", 
    "SNSoln", "solInt", "timeAvg", "chAvg", "EVPA", "RM", 
    "minSNR", "prtLv", "nThreads",
    NULL};
  gchar *routine = "XYDlyHistory";

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
 
} /* end XYDlyHistory  */

/**
 * Determine cross pol delay and phase
 * \param myInput Input parameters on InfoList    
 * \param avgUV   Averaged data to be used to determine calibration
 * \param inUV    UV data onto which the SN table to be written
 * \param err     ObitErr stack for reporting problems.
 */
void  XYDelayCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		 ObitErr* err)
{
  ObitTableSN *SNTable=NULL;
  gchar        *dataParms[] = {  /* Parameters to control calibration */
    "Sources", "UVR_Full", "WtUV", "minSNR", "timeAvg", "chAvg", 
    "SNSoln", "solInt", "EVPA", "RM", "refAnt", "fitType", 
     NULL};
  gchar *routine = "XYDelayCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(avgData));


  /* Any data in avgData? */
  if (avgData->myDesc->nvis<=0) {
    Obit_log_error(err, OBIT_Error, "NO averaged/calibrated data");
    return;
  }

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, avgData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Fit - write SN table */
  SNTable = ObitUVXYDelayCal (avgData, inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Cleanup */
  SNTable = ObitTableSNUnref(SNTable);
} /* end XYDelayCal  */

/**
 * Determine cross polarized delays and phase for an UV dataset divided by
 * a source model.
 * All selected spectra are averaged and the fitting is done to these spectra.
 * Corrections are added to the second polarization in the output SN table.
 * \param inUV   Input UV data. 
 * Control parameters are on the info member.
 * \li "Sources" OBIT_str  (16,10,1) Selected Sources
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "UVR_Full" OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (lamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "SNSoln"  OBIT_long  (1,1,1)  Output SN table, 0=>new
 * \li "solInt"  OBIT_float  (1,1,1) Solution interval (min), def [long]
 * \li "EVPA"    OBIT_float  (10,1,1) EVPA (deg) at ref. freq. of Sources
 * \li "RM"      OBIT_float  (10,1,1) RM (ras/m^2)  of Sources
 * \li "fitType" OBIT_long  (1,1,1) fit Type def=0
 *                                   0=Both, 1=XY,2=YX
 * \param outUV  UV with which the output  SN is to be associated
 * Control parameters are on the info member.
 * \li "BChan"   OBIT_int   (1,1,1) BChan used to generate inUV (default 1)
 * \param err    Error/message stack, returns if error.
 * \return Pointer to the newly created SN object which is associated with outUV.
 */
ObitTableSN* ObitUVXYDelayCal (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *row=NULL;
  ObitTableAN *ANTab=NULL;
  ObitAntennaList *ANList=NULL;
  ObitSource *Source=NULL;
  ObitSourceList *SList=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  gboolean doCalSelect=TRUE;
  ObitIOCode retCode;
  ObitInfoType type;
  olong i, j, k, kndx, ii, iAnt, SNver, suba, refAnt, BChan, SNSoln, suID, lastSU=-10;
  olong iSNRow, numFreq, jndx, indx, numPol, numIF, numPCal, numAnt;
  olong incs, incf, incif, js, jf, jif, numVis, numOrb, iANver, cntTime;
  olong fitType, *lookupSU = NULL;
  ofloat uvrang[2], wtuv, weight, wt, dFreq, snrmin, bl, chi, chisin, chicos;
  ofloat p1, p2, cp1, cp2, tr, ti, solInt, phase0, *flip=NULL;
  ofloat *vis, *u, *v, *base, *time, *delay=NULL, *phase=NULL, *snr=NULL;
  ofloat *antwt=NULL, *EVPA=NULL, *EVPAFS=NULL, *EVPAFC=NULL, *RM=NULL;
  ofloat *xpol1=NULL, *xpol2=NULL;
  ofloat sumTime, startTime, endTime=0.0, lastTime=0.0;
  ofloat amp1, amp2;
  odouble *avgIFFreq=NULL;
  gboolean empty, done, dump, multiSU=FALSE;
  gchar timeString[25], *tname;
  gchar *routine = "ObitUVXYDelayCal";
  /* error checks */
  if (err->error) return outSoln;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Make sure you can solve for delay */
  Obit_retval_if_fail(((inUV->myDesc->inaxes[inUV->myDesc->jlocf]>1)),
		       err, outSoln,
		      "%s: MUST have multiple frequencies to solve for Delay", 
		      routine);

  /* Make sure 4 Stokes correlations are available */
  Obit_retval_if_fail(((inUV->myDesc->inaxes[inUV->myDesc->jlocs]>=4)),
		       err, outSoln,
		      "%s: MUST have cross polarized data", 
		      routine);

  
  /* Selection of input? should have been done before */
  access = OBIT_IO_ReadCal;

  /* open UV data  */
  retCode = ObitUVOpen (inUV, access, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  /* Update frequency tables on inUV */
  if (!inUV->myDesc->freqArr) ObitUVGetFreq (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);

  /* Need array information */
  if (!outUV->myDesc->numAnt)   ObitUVGetSubA (outUV, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
  
   /* open output UV data  */
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);

  /* Better be linear feeds? */
  Obit_retval_if_fail((inUV->myDesc->crval[inUV->myDesc->jlocs]<-4.0),
		      err, outSoln, "%s: MUST be linear feed data", routine);

  /* Antenna List for linear feeds */
  numOrb  = 0;  numPCal = 0; numIF = 0; iANver=1;
  ANTab = 
    newObitTableANValue (inUV->name, (ObitData*)inUV, &iANver, OBIT_IO_ReadOnly, 
			 numIF, numOrb, numPCal, err);
  ANList = ObitTableANGetList (ANTab, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);

  /* Get parameters from inUV */
  snrmin = 5.0;
  ObitInfoListGetTest(inUV->info, "minSNR", &type, dim, &snrmin);
  uvrang[0] = 0.0; uvrang[1] = 1.0e15;
  if (!ObitInfoListGetTest(inUV->info, "UVR_Full", &type, dim, &uvrang[0])
      || (uvrang[1]<=uvrang[0])) {
    /* If no explicit uv range given, use everything */
    uvrang[0] = 0.0;
    uvrang[1] = 1.0e15;
  }  /* end derive uv range */
  wtuv = 1.0;
  ObitInfoListGetTest(inUV->info, "WtUV", &type, dim, &wtuv);
  /* rlp = RLPhase/57.296;  R-L phase in radians 
  qpol = cos(rlp); upol = sin(rlp);   Fake poln for lin. feeds */
  refAnt = 1;
  ObitInfoListGetTest(inUV->info, "refAnt", &type, dim, &refAnt);
  fitType = 0;
  ObitInfoListGetTest(inUV->info, "fitType", &type, dim, &fitType);
  BChan = 1;
  ObitInfoListGetTest(outUV->info, "BChan", &type, dim, &BChan);
  BChan -= 1; /* Zero rel */
  ObitInfoListGetP(inUV->info, "RM", &type, dim, (gpointer)&RM);
  ObitInfoListGetP(inUV->info, "EVPA", &type, dim, (gpointer)&EVPA);
 
  SNSoln = 0; /* Output SN table version */
  ObitInfoListGetTest(outUV->info, "SNSoln", &type, dim, &SNSoln);
  solInt = 1440.0*1.0e5;   /* Very long default solution interval */
  ObitInfoListGetTest(outUV->info, "solInt", &type, dim, &solInt);
  solInt /= 1440.0; /* to days */
  if (solInt<=0.0) solInt = 1.0e5;   /* Very long default solution interval */
 
  /* Output SN table  */
  if (SNSoln>0) SNver = SNSoln;
  else SNver = 1 + ObitTableListGetHigh (outUV->tableList, "AIPS SN");
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
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
  if (inUV->myDesc->jlocf>=0)
    numFreq  = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  else numFreq  = 1;
  dFreq = inUV->myDesc->cdelt[inUV->myDesc->jlocf];

  /* Close - this should update header */
  retCode = ObitUVClose (outUV, err);
  if (err->error) goto cleanup;

  if (err->prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, "Writing cross poln gains in SN table %d",SNver);
  }

  /* Increments */
  incs  = inUV->myDesc->incs;
  incf  = inUV->myDesc->incf;
  incif = inUV->myDesc->incif;
 
  /* Create arrays - pad */
  xpol1  =  g_malloc0(numFreq*(numIF+2)*3*sizeof(ofloat));
  xpol2  =  g_malloc0(numFreq*(numIF+2)*3*sizeof(ofloat));
  delay  =  g_malloc0((numIF+2)*sizeof(ofloat));
  phase  =  g_malloc0((numIF+2)*sizeof(ofloat));
  snr    =  g_malloc0((numIF+2)*sizeof(ofloat));
  EVPAFS =  g_malloc0((numFreq*(numIF+2))*sizeof(ofloat));
  EVPAFC =  g_malloc0((numFreq*(numIF+2))*sizeof(ofloat));
  flip   =  g_malloc0((numFreq*(numIF+2))*sizeof(ofloat));
  for (j=0; j<numFreq*numIF; j++) flip[j] = 1.0;  /* Default no flip */
  /* zero accumulators */
  for (j=0; j<numFreq*(numIF+2)*3; j++) xpol1[j] = xpol2[j] = 0.0;
  for (j=0; j<numIF+2; j++) delay[j] = phase[j] = snr[j] = 0.0;
  
  multiSU = (inUV->myDesc->ilocsu>=0); /* MultiSource? */
  if (!multiSU) {
    suID = -1;
    Source= ObitUVGetSource (Source, inUV, suID, err);
    /* Source lookup table */
    lookupSU   = g_malloc0(2*sizeof(olong));
    lookupSU[0]= 0;
    /* Setup EVPA arrays */
    SetupEVPA(numIF, numFreq, inUV, EVPA[0], RM[0], EVPAFS, EVPAFC);
  } else { /* end single source */
    /* Multi source, need sourceList, lookup table */
    SList = ObitUVGetSourceList (inUV, err); /* Get full Source List */
    lookupSU = GetLookupSU (SList, inUV->info);
  } /* end multiSource */
  if (err->error) Obit_traceback_val (err, routine, inUV->name,outSoln);

  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(inUV->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_retval_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err, outSoln,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);
  
  numAnt = outUV->myDesc->numAnt[suba-1];
  
  /* wavelength at reference frequency 
     lambda0 = VELIGHT/inUV->myDesc->freq;*/

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
  
  /* Initialize solution row - add corrections in second poln */
  row->Time   = 0.5; row->TimeI  = 2.0; row->SourID = 0; row->antNo  = 0; row->SubA   = 0; 
  row->FreqID = 1; row->IFR    = 0.0; row->NodeNo = 0; row->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    row->Real1[i]   = 1.0; row->Imag1[i]   = 0.0; 
    row->Delay1[i]  = 0.0; row->Rate1[i]   = 0.0; 
    row->RefAnt1[i] = refAnt;
    if (snr[i]>=snrmin) {
      row->Weight1[i] = snr[i]; 
    } else {
      /* No way to flag Xpol only here */
      row->Weight1[i] = 1.0; 
    } 
  }
  if (numPol>1) {
    row->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      row->Rate2[i]   = 0.0; 
      row->RefAnt2[i] = refAnt; 
      row->Real2[i]   = 1.0; row->Imag2[i]   = 0.0;
      row->Delay2[i]  = 0.0; row->Weight2[i] = 1.0;
    }
  }

  /* Loop through data */
  startTime = -1.0e10; cntTime=0; sumTime=0.0;
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    if (doCalSelect) retCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else retCode = ObitUVRead (inUV, inUV->buffer, err);
    done =  (retCode!=OBIT_IO_OK);
    if (err->error) goto cleanup;
    
    /* initialize data pointers */
    vis  = inUV->buffer+inUV->myDesc->nrparm;
    u    = inUV->buffer+inUV->myDesc->ilocu;
    v    = inUV->buffer+inUV->myDesc->ilocv;
    time = inUV->buffer+inUV->myDesc->iloct;
    base = inUV->buffer+inUV->myDesc->ilocb;

    /* Time range for first integration */
    if (startTime<-10.) {
      startTime = *time;
      endTime   = *time + solInt;
      lastTime  = *time;
    }
    /* Get source info for lin. feeds - may not need this, needs work if so */
    if (multiSU && (lastSU!=inUV->buffer[inUV->myDesc->ilocsu])){
      suID = inUV->buffer[inUV->myDesc->ilocsu];
      Source = ObitUVGetSource (Source, inUV, suID, err);
      if (err->error) goto cleanup;
      /* Setup EVPA arrays */
      SetupEVPA(numIF, numFreq, inUV, EVPA[lookupSU[suID]], RM[lookupSU[suID]], 
		EVPAFS, EVPAFC); 
    } /* end new source */

    /* New time? */
    if (time[0] > lastTime) {
      cntTime++; sumTime+=time[0];  /* For average time */
      /* linear feeds - flip phase if linP (eff. Stokes u) fn positive */
      chi  = ObitAntennaListParAng(ANList, 1,  time[0], Source);  /* Parallactic angle */
      sincosf(chi*2.0, &chisin, &chicos);  /* Need sine/cosine of 2*chi */
      /*sincosf(chi, &chisin, &chicos);  DEBUG Need sine/cosine of 2*chi */
      lastTime = time[0];
      /* Calculate phase flips */
      GetFlip(numIF, numFreq, inUV, EVPAFS, EVPAFC, chisin, chicos, flip);
      /* DEBUG 
	 day2dhms(time[0], timeString);
	 fprintf(stderr,"time=%s chi=%f\n",timeString,chi); DEBUG */
    }
	  
    numVis = inUV->myDesc->numVisBuff;
    if (done) numVis = MAX(1,numVis);
    for (i=0; i<numVis; i++) { /* loop over buffer */
      
      /* Time to write solutions? */
      dump = *time>endTime;
      dump =  dump || done;
      if (dump) {
	/* Get IF Average frequencies if not done */
	if (avgIFFreq==NULL) 
	  avgIFFreq = GetAvgFreq (numIF, numFreq, xpol1, xpol2, inUV, err);
	if (err->error) Obit_traceback_val (err, routine, outUV->name, outSoln);
	
	/* Time info */
	if (err->prtLv>=1) {
	  if (cntTime>0) day2dhms(sumTime/cntTime, timeString);
	  else           day2dhms(time[0], timeString);
	  Obit_log_error(err, OBIT_InfoErr, "Solution at %s", timeString);
	}
	/* Data all loaded - Process,  Fit delays per IF */
	for (jif=0; jif<numIF; jif++) {   /* IF loop */
	  XYFitDelay (numFreq, BChan, dFreq, &xpol1[jif*numFreq*3], &xpol2[jif*numFreq*3], 
		      fitType, &delay[jif], &phase[jif], &snr[jif], err);
	  /* Phase to IF reference frequency */
	  /*delay[jif] = -delay[jif]; phase[jif] = -phase[jif];  DEBUG */
	  phase0 = phase[jif];
	  phase[jif] -= 2.0*G_PI*(avgIFFreq[jif]-inUV->myDesc->freqIF[jif])*delay[jif];
	  
	  /* Diagnostics */
	  if (err->prtLv>=2) {
	    if (snr[jif]<snrmin) continue;
	    for (ii=0; ii<numFreq; ii++) {
	      tr   = xpol1[(jif*numFreq+ii)*3];
	      ti   = xpol1[(jif*numFreq+ii)*3+1];
	      if ((tr==0.0) && (ti==0.0)) continue;
	      p1   = 57.295*atan2(ti, tr);
	      cp1  = p1 + fmod(360.0*(delay[jif]*dFreq)*(ii+BChan),360.0) - phase0*57.296;
	      if (xpol1[(jif*numFreq+ii)*3+2]>0.0) {
		tr /= xpol1[(jif*numFreq+ii)*3+2]; ti /= xpol1[(jif*numFreq+ii)*3+2];
		amp1 = sqrtf(tr*tr+ti*ti);
	      } else amp1 = 0.0;
	      /* channel XY residual */
	      if (cp1> 180.0) {cp1 -= 360.0;} if (cp1> 180.0) {cp1 -= 360.0;}
	      if (cp1<-180.0) {cp1 += 360.0;} if (cp1<-180.0) {cp1 += 360.0;}
	      tr   = xpol2[(jif*numFreq+ii)*3];
	      ti   = xpol2[(jif*numFreq+ii)*3+1];
	      p2   = 57.295*atan2(ti, tr);
	      if (xpol2[(jif*numFreq+ii)*3+2]>0.0) {
		tr /= xpol2[(jif*numFreq+ii)*3+2]; ti /= xpol2[(jif*numFreq+ii)*3+2];
		amp2 = sqrtf(tr*tr+ti*ti);
	      } else amp2 = 0.0;
	      /* channel YX residual */
	      cp2  = p2 + fmod(360.0*(delay[jif]*dFreq)*(ii+BChan),360.0) - phase0*57.296; 
	      if (cp2> 180.0) {cp2 -= 360.0;} if (cp2> 180.0) {cp2 -= 360.0;}
	      if (cp2<-180.0) {cp2 += 360.0;} if (cp2<-180.0) {cp2 += 360.0;}
	      Obit_log_error(err, OBIT_InfoErr, 
			     "ch %4d xc1 %8.2f xc2 %8.2f corr %8.2f xc1 %8.2f xc2 %8.2f a1 %8.5f a2 %8.5f ",
			     BChan+ii, p1, p2, 360.0*(delay[jif]*dFreq)*(ii+BChan), cp1, cp2, amp1, amp2);
	    } /* end channel loop */
	  }
	  if (err->prtLv>=1) {
	    Obit_log_error(err, OBIT_InfoErr, "IF %2d delay %8.2f nsec phase %8.2f SNR %5.1f",
			   jif+1, delay[jif]*1.0e9, phase[jif]*57.296, snr[jif]);
	    /*     		     jif+1, delay[jif]*dFreq, phase[jif]*57.296, snr[jif]);*/
	  }
	  /* Minimum SNR */
	  if (snr[jif]<snrmin) snr[jif] = 0.0;
	  
	  /* Take into account BChan
	     phase[jif] -= 2.0*G_PI * dFreq*BChan*delay[jif]; */
	  /* Set solutions into table */
	  if (cntTime>0) row->Time = sumTime/cntTime; 
	  else           row->Time = time[0];
	  row->TimeI  = solInt;
	  if (numPol>1) {
	    if (snr[jif]>snrmin) {  
	      row->Real2[jif]  =  cos(phase[jif]); row->Imag2[jif]   = -sin(phase[jif]); 
	      row->Delay2[jif] = -delay[jif];      row->Weight2[jif] =  snr[jif]; 
	    } else {
	      row->Real2[jif]  = 1.0; row->Imag2[jif]   = 0.0;
	      row->Delay2[jif] = 0.0; row->Weight2[jif] = 1.0; 
	    }
	  }
	} /* end IF loop */
	
	/* Write - Loop over antennas - write as corrections rather than solutions */
	iSNRow = -1;
	for (iAnt= 0; iAnt<numAnt; iAnt++) {
	  row->antNo  = iAnt+1; 
	  retCode = ObitTableSNWriteRow (outSoln, iSNRow, row, err);
	  if (err->error) goto cleanup;
	} /* end antenna loop */
	
	if (done) goto closeup;  /* All data read/processed? */

 	/* zero accumulators */
	for (j=0; j<numFreq*(numIF+2)*3; j++) xpol1[j] = xpol2[j] = 0.0;
 	for (j=0; j<numIF+2; j++) delay[j] = phase[j] = snr[j] = 0.0;
	/* Next timerange */
	startTime = *time; cntTime=0; sumTime=0.0;
	endTime   = *time + solInt;
      } /* End dump solution */

      /* accumulate this visibility - Get baseline length */
      bl = sqrt ((*u)*(*u) + (*v)*(*v));

      /* Weighting */
      wt = 1.0;
      if ((bl<uvrang[0]) || (bl>uvrang[1])) wt = wtuv;
      
      /* first cross pol XY - by channel then IF */
      if (fitType<=1) {  /* wanted for fitting? */
	js = 2; k = 0;
	for (jif=0; jif<numIF; jif++) {   /* IF loop */
	  for (jf=0; jf<numFreq; jf++) {  /* Frequency loop */
	    jndx = 3*(jf + jif*numFreq);
	    indx = js*incs + jif*incif + jf*incf;
	    kndx = jif*numIF + jf;
	    weight = vis[indx+2] * wt;
	    if (weight>0.0) {
	      xpol1[jndx+2] += weight;
	      if (flip[kndx]>0) { /* Flip sign of one of the real/imag parts */
		xpol1[jndx]   -= vis[indx] * weight;
		xpol1[jndx+1] += vis[indx+1]*weight;
	      } else {
		xpol1[jndx]   += vis[indx] * weight;
		xpol1[jndx+1] -= vis[indx+1]*weight;
	      }
	    }
	    k++; /* flip index */
	  } /* end Frequency loop */
	} /* end IF loop */
      } /* End if selected by fit type */
      
      /* Second cross YX pol */
      if ((fitType==0) || (fitType==2)) { /* wanted for fitting? */
	js = 3; k = 0;
	for (jif=0; jif<numIF; jif++) {   /* IF loop */
	  for (jf=0; jf<numFreq; jf++) {  /* Frequency loop */
	    jndx = 3*(jf + jif*numFreq);
	    indx = js*incs + jif*incif + jf*incf;
	    kndx = jif*numIF + jf;
	    weight = vis[indx+2] * wt;
	    if (weight>0.0) {
	      xpol2[jndx+2] += weight;
	      if (flip[kndx]>0) { /* Flip sign of one of the real/imag parts 
				     always another flip of the imaginary */
		xpol2[jndx]   -= vis[indx] * weight;
		xpol2[jndx+1] -= vis[indx+1]*weight;
	      } else {
		xpol2[jndx]   += vis[indx] * weight;
		xpol2[jndx+1] += vis[indx+1]*weight;
	      }
	    }
	    k++; /* flip index */
	  } /* end Frequency loop */
	} /* end IF loop */
      } /* End if selected by fit type */
      
      /* update data pointers */
      vis  += inUV->myDesc->lrec;
      u    += inUV->myDesc->lrec;
      v    += inUV->myDesc->lrec;
      base += inUV->myDesc->lrec;
      time += inUV->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
  /* Close */
 closeup:
  retCode = ObitUVClose (inUV, err);
  if (err->error) goto cleanup;

  /* Close output table */
  retCode = ObitTableSNClose (outSoln, err);
  if (err->error) goto cleanup;
  
 cleanup: 
  row    = ObitTableSNUnref(row);
  ANTab  = ObitTableANUnref(ANTab);
  ANList = ObitAntennaListUnref(ANList);
  Source = ObitSourceUnref(Source);

  if (antwt)  g_free(antwt);
  if (xpol1)  g_free(xpol1);
  if (xpol2)  g_free(xpol2);
  if (delay)  g_free(delay);
  if (phase)  g_free(phase);
  if (snr)    g_free(snr);
  if (flip)   g_free(flip);
  if (lookupSU)  g_free(lookupSU);
  if (EVPAFS)    g_free(EVPAFS);
  if (EVPAFC)    g_free(EVPAFC);
  if (avgIFFreq) g_free(avgIFFreq);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  return outSoln;
} /* end ObitUVXYDelayCal */

/**
 * Determine cross polarized delays and phase for given complex spectrum
 * Use direct parameter search.
 * If an amplitude of either xpol1 or xpol2 is negative, the signs of both
 * real and imaginary parts are flipped to maintain the correct orientation.
 * \param numFreq  Number of frequencies
 * \param BChan    First 0 rel channel selected
 * \param dFreq    Frequency increment in Hz
 * \param xpol1    weighted visibility spectrum, p1 x p2
 * \param xpol2    weighted visibility spectrum, p2 x p1
 * \param fitType  Fitting type, 0=joint, 1=XY only, 2=LX only
 * \param delay    [out] Delay (nsec)
 * \param phase    [out] phase (rad)
 * \param snr      [out] SNR determined from scatter of phase
 * \param err      Error/message stack, returns if error.
 */
void XYFitDelay (olong numFreq, olong BChan, ofloat dFreq,
		 ofloat *xpol1, ofloat *xpol2, olong fitType, 
		 ofloat *delay, ofloat *phase, ofloat *snr, ObitErr *err)
{
  olong i, j, cnt, cnt1, cnt2, ntry;
  ofloat td1, td2, td3, td4, d1, d2, d3, d4, test=0, best, best2, best3, best4;
  ofloat corr, cori, tr, ti, w, p, tp, dp;
  ofloat *xp1=NULL, *xp2=NULL;
  odouble sumr=0, sumi=0, sumw=0, sumr1=0, sumi1=0, sumw1=0, sumr2=0, sumi2=0, sumw2=0;
  /* error checks */
  if (err->error) return;

  /* init output */
  *delay = 0.0;
  *phase = 0.0;
  *snr   = 0.0;
  /* For now, average phase 
     for (i=0; i<numFreq; i++) {
     sumr += xpol1[i*3] + xpol2[i*3]; sumi += xpol1[i*3+1] + xpol2[i*3+1]; 
     }
     *phase = atan2(sumi,sumr);
     if (abs(sumi)>0.0) *snr = 10;
     return;*/

  /* Condition xpoln */
  for (i=0; i<numFreq; i++) {
    if (xpol1[i*3]<0.0) {
      xpol1[i*3] = -xpol1[i*3]; xpol1[i*3+1] = -xpol1[i*3+1]; 
    }
    if (xpol2[i*3]<0.0) {
      xpol2[i*3] = -xpol2[i*3]; xpol2[i*3+1] = -xpol2[i*3+1]; 
    }
  } /* end frequency loop */
 
  /* How many attempts */
  ntry =  51;
 
 /* Which input being used? */
  xp1 = xp2 = NULL;
  if (fitType==0) {xp1 = xpol1; xp2 = xpol2;}
  if (fitType==1) {xp1 = xpol1;}
  if (fitType==2) {xp2 = xpol2;}
  
  /* Iterations of direct parameter search - amplitude of weighted sum 
     sequence of ever finer searches */
  d1   = 0.0;
  best = -1.0e20;
  for (i=0; i<ntry; i++) {
    td1 = (i-ntry/2) * 0.01;  /* delay in turns per channel */
    test = XYtestDelay (numFreq, xp1, xp2, td1);
    if (test>best) {
      best = test;
      d1    = td1;
    }
  }  
  d2    = d1;
  best2 = best;
  for (i=0; i<ntry; i++) {
    td2 = d1 + (i-ntry/2) * 0.003;  /* delay in turns per channel */
    test = XYtestDelay (numFreq, xp1, xp2, td2);
    if (test>best2) {
      best2 = test;
      d2    = td2;
    }
  }
  
  d3    = d2;
  best3 = best2;
  for (i=0; i<15; i++) {
    td3 = d2 + (i-7) * 0.001;  /* delay in turns per channel */
    test = XYtestDelay (numFreq, xp1, xp2, td3);
    if (test>best3) {
      best3 = test;
      d3    = td3;
    }
  }

  d4    = d3;
  best4 = best3;
  for (i=0; i<361; i++) {
    td4 = d3 + (i-180) * 0.000003;  /* delay in turns per channel */
    test = XYtestDelay (numFreq, xp1, xp2, td4);
    if (test>best4) {
      best4 = test;
      d4    = td4;
      /*fprintf(stderr,"test %f delay %f\n",best,1.0e9*td4/dFreq);*/
    }
  }

  /* Fitted delay - unalias */
  if (d4>0.5)  d4 -= 1.0;
  if (d4<-0.5) d4 += 1.0;
  *delay = d4 / dFreq;

  /* Average phase */
  j = 0;
  sumr1 = sumi1 = sumw1 = sumr2 = sumi2 = sumw2 = 0.0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*d4*(i+BChan));
    cori = sin(2.0*G_PI*d4*(i+BChan));
    /* First poln */
    if (xp1) {
      tr = xp1[j+0]*corr - xp1[j+1]*cori;
      ti = xp1[j+0]*cori + xp1[j+1]*corr;
      w  = xp1[j+2];
      sumr1 += tr * w;
      sumi1 += ti * w;
      sumw1 += w;
    }
    /* Second poln */
    if (xp2) {
      tr = xp2[j+0]*corr - xpol2[j+1]*cori;
      ti = xp2[j+0]*cori + xpol2[j+1]*corr;
      w  = xp2[j+2];
      sumr2 += tr * w;
      sumi2 += ti * w;
      sumw2 += w;
    }
    j += 3;
  } /* end channel loop */
  /* combine */
  sumr = sumr1 + sumr2;
  sumi = sumi1 + sumi2;
  sumw = sumw1 + sumw2;
      
  if (sumw<=0.0) return;  /* Anything? */
  
  /* Average phase */
  p = atan2 (sumi, sumr);
  *phase = p;

   /* RMS difference */
  sumr = sumr1 = sumr2 = 0.0;
  cnt  = cnt1  = cnt2  = 0;
  j    = 0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*d4*(i+BChan));
    cori = sin(2.0*G_PI*d4*(i+BChan));
    /* First poln */
    if (xp1) {
      tr = xp1[j+0]*corr - xp1[j+1]*cori;
      ti = xp1[j+0]*cori + xp1[j+1]*corr;
      if (xp1[j+2]>0.0) {
	tp = atan2 (ti, tr);
	dp = (tp - p);
	if (dp> G_PI) {dp -= 2.0*G_PI;} if (dp<-G_PI) {dp += 2.0*G_PI;}
	sumr1 += dp*dp;
	cnt1++;
      }
    }
    /* Second poln  */
    if (xp2) {
      tr = xp2[j+0]*corr - xp2[j+1]*cori;
      ti = xp2[j+0]*cori + xp2[j+1]*corr;
      if (xpol2[j+2]>0.0) {
	tp = atan2 (ti, tr);
	dp = (tp - p);
	if (dp> G_PI) {dp -= 2.0*G_PI;} if (dp<-G_PI) {dp += 2.0*G_PI;}
	sumr2 += dp*dp;
	cnt2++;
      }
    }
    j += 3;
  } /* end channel loop */

  /* Combine */
  sumr = sumr1 + sumr2;
  cnt  = cnt1  + cnt2;
  /* SNR from variance of phase residuals */
  if (cnt>=2) *snr = 1.0 / sqrt(sumr/(cnt-1));
  
  return;
} /* end XYFitDelay  */

/**
 * Get amplitude^2 of weighted average with test delay
 * Use direct parameter search
 * \param numFreq  Number of frequencies
 * \param xpol1    weighted visibility spectrum, p1 x p2, if NULL, ignore
 * \param xpol2    weighted visibility spectrum, p2 x p1, if NULL, ignore
 * \param delay    Delay turns per channel
 * \return  amplitude^2 of weighted average after correcting spectra
 */
ofloat XYtestDelay (olong numFreq, ofloat *xpol1, ofloat *xpol2, ofloat delay)
{
  ofloat out = 0.0;
  olong i, j;
  ofloat corr, cori, tr, ti, w;
  odouble sumr, sumi, sumw;

  sumr = sumi = sumw = 0.0;
  j = 0;
  for (i=0; i<numFreq; i++) {
    /* Delay correction */
    corr = cos(2.0*G_PI*delay*i);
    cori = sin(2.0*G_PI*delay*i);
    /* First poln */
    if (xpol1) {
      tr = xpol1[j+0]*corr - xpol1[j+1]*cori;
      ti = xpol1[j+0]*cori + xpol1[j+1]*corr;
      w  = xpol1[j+2];
      sumr += tr * w;
      sumi += ti * w;
      sumw += w;
    }
    /* Second poln */
    if (xpol2) {
      tr = xpol2[j+0]*corr - xpol2[j+1]*cori;
      ti = xpol2[j+0]*cori + xpol2[j+1]*corr;
      w  = xpol2[j+2];
      sumr += tr * w;
      sumi += ti * w;
      sumw += w;
    }
    j += 3;
  }
  /* anything? */
  if (sumw<=0.0) return out;

  /* Average */
  sumr /= sumw;
  sumi /= sumw;
  out = sumr*sumr + sumi*sumi;
  return out;
}  /* end XYtestDelay */

/**
 * Get Average frequency per IF of frequencies with weights >0
 * in xpol1 and xpol2.
 * \param numIF    Number of IFs
 * \param numFreq  Number of frequencies per IF
 * \param xpol1    weighted visibility spectrum, XY
 * \param xpol2    weighted visibility spectrum, YX
 * \param delay    Delay turns per channel
 * \param inUV     Input UV data. 
 * \param err      Error/message stack, returns if error.
 * \return  array of weighted frequencies, g_free when done.
 */
odouble* GetAvgFreq (olong numIF, olong numFreq, ofloat *xpol1, ofloat *xpol2, 
		    ObitUV *inUV, ObitErr *err)
{
  odouble *out = NULL;
  olong i, j, off1, off2;
  odouble sumFreq, sumWt;

  if (err->error) return out;  /* Previous error? */

  /* Create outut array - some padding */
  out = g_malloc0((numIF+2)*sizeof(odouble));
  for (j=0; j<numIF; j++) { /* IF Loop */
    out[j] = inUV->myDesc->freqIF[j]; /* in case no data */
    sumFreq = sumWt = 0.0;
    off1 = j*numFreq*3;  /* offset in xpol? */
    off2 = j*numFreq;    /* offset in freqArr */
    for (i=0; i<numFreq; i++) {
      /* Anybody home? */
      if ((xpol1[off1+i*3+2]<=0.0) || (xpol1[off1+i*3+2]<=0.0)) break;
      sumFreq += xpol1[off1+i*3+2] * inUV->myDesc->freqArr[off2+i];
      sumWt   += xpol1[off1+i*3+2];
    } /* end channel loop */
    if (sumWt>0.0) out[j] = sumFreq/sumWt;
  }  /* end  IF Loop */

  return out;
} /* end GetAvgFreq */

/**
 * Get source lookup table, translates Source IDs to index in 
 * Sources ordered arrays.
 * \param SList    Source list
 * \param info     ObitInfoList with "Sources"
 * \param err      Error/message stack, returns if error.
 * \return  source lookup table
 */
olong* GetLookupSU (ObitSourceList *SList, ObitInfoList *info)
{
  olong *lookupSU = NULL;
  olong i, j, k, len;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar *Sources=NULL;

  /* Source lookup table */
  lookupSU = g_malloc0((SList->number+5)*sizeof(olong));
  if (ObitInfoListGetP(info, "Sources", &type, dim, (gpointer)&Sources)) {
    for (i=0; i<SList->number; i++) {  /* Loop over SList */
      ObitTrimTrail(SList->SUlist[i]->SourceName);
      len = strlen(SList->SUlist[i]->SourceName);
      k = 0;
      for (j=0; j<dim[1]; j++) {
	if (!strncmp(SList->SUlist[i]->SourceName,&Sources[k],len)) {
	  lookupSU[SList->SUlist[i]->SourID] = j; break;
	} /* end of if match */
	k += dim[0];
      } /* end loop over input sources */
    } /* end loop over SList */
  } else { /* end if "Sources" */
    /* Fooey no "Sources" - shouldn't get here */
    for (i=0; i<SList->number; i++)  lookupSU[0] = i;
  } /* end fake Source lookup table */
  return lookupSU;
} /* end GetLookupSU */

/**
 * Setup EVPA sin/cos arrays for a source
 * \param numIF    Number of IFs
 * \param numFreq  Number of channels per IF
 * \param inUV     Input UV data
 * \param EVPA     EVPA (deg) are reference frequency
 * \param RM       RM (rad/m^2)
 * \param EVPAFS   Sine of EVPA in each frequency channel
 * \param EVPAFC   Cosine of EVPA in each frequency channel
 */
void SetupEVPA (olong numIF, olong numFreq, ObitUV *inUV, 
		ofloat EVPA, ofloat RM, ofloat *EVPAFS, ofloat *EVPAFC)
{
  olong    nIFCh, i;
  ofloat   *args=NULL, EVPArad = EVPA*DG2RAD;
  odouble  lambda_2, lambda0_2;

  nIFCh = numIF*numFreq;
  args = g_malloc0((nIFCh+10)*sizeof(ofloat));  /* work array */
  lambda0_2 = (VELIGHT/inUV->myDesc->freq)*(VELIGHT/inUV->myDesc->freq);
  /* Calculate EVPA at each channel */
  for (i=0; i<nIFCh; i++) {
   lambda_2 = (VELIGHT/inUV->myDesc->freqArr[i])*(VELIGHT/inUV->myDesc->freqArr[i]);
   args[i] = EVPArad + RM * (lambda_2 - lambda0_2);
  }
  /* Take sine/cosine */
  ObitSinCosInit();   /* to be sure */
  ObitSinCosVec(nIFCh, args, EVPAFS, EVPAFC);

  if (args) g_free(args);
} /* end SetupEVPA */
/**
 * Setup EVPA sin/cos arrays for a source
 * \param numIF    Number of IFs
 * \param numFreq  Number of channels per IF
 * \param inUV     Input UV data
 * \param EVPAFS   Sine of EVPA in each frequency channel
 * \param EVPAFC   Cosine of EVPA in each frequency channel
 * \param chisin   Sine of twice parallactic angle
 * \param chicos   Cosine of twice parallactic angle
 *       Output
 * \param flip     1 or -1 depending on whether the phase need to be 
 *                 flipped in that channel.
 */
void GetFlip (olong numIF, olong numFreq, ObitUV *inUV, ofloat *EVPAFS, ofloat *EVPAFC, 
	      ofloat chisin, ofloat chicos, ofloat *flip)
{
  olong    nIFCh, i;
  ofloat   linP;
  gboolean flipit=FALSE;

  nIFCh = numIF*numFreq;
  for (i=0; i<nIFCh; i++) {
    /* Rotate EVPA by chi - THIS SHOULD BE CHECKED */
    linP = -EVPAFC[i] * chisin + EVPAFS[i] * chicos;
    flipit = linP>0.0; /* eff. Stokes u) fn positive */
    /*flipit = linP<0.0; /eff. Stokes u) fn negative */
    /* ObitrLDly has 
       otherwise phase[jif] = -G_PI-phase[jif];
       phase[jif] = fmod(phase[jif],2.0*G_PI);
       Not sure abut this.
     */
    /*flipit=FALSE; DEBUG*/
    if (flipit) flip[i] = -1.0;
    else        flip[i] = +1.0;
  }
  /*fprintf (stderr, "linP=%f flipit=%d ", linP,flipit);  DEBUG */
} /* end SetupEVPA */

/** 
 * Convert Time in days to a human readable form "dd/hh:mm:ss.s"
 * \param time  Time in days
 * \param timeString [out] time as string, should be >16 char
 */
void day2dhms(ofloat time, gchar *timeString)
{
  olong day, thour, tmin;
  ofloat ttim, ssec;

  /* Trap bad times */
  if ((time<-100.0) || (time>1000.0)) {
    sprintf (timeString, "Bad time");
    return;
  }

  day   = (olong)(time);
  ttim  = 24.0*(time - day);
  thour = MIN ((olong)(ttim), 23);
  ttim  = 60.0*(ttim - thour);
  tmin  = MIN ((olong)(ttim), 59);
  ssec  = 60.0*(ttim - tmin);
  /* avoid silliness */
  if (ssec>59.951) {
    tmin++;
    ssec = 0.0;
  }
  if (tmin>=60) {
    thour++;
    tmin -= 60;
  }
  if (thour>23) {
    day++;
    thour -= 24;
  }
  sprintf (timeString, "%2.2d/%2.2d:%2.2d:%4.1f", day, thour, tmin, ssec);
  /* Zero fill seconds */
  if (timeString[9]==' ') timeString[9] = '0';
} /* end day2dhms */

