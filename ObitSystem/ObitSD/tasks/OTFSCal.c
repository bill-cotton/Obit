/* $Id:  $  */
/* Obit  Radio Single dish On The Fly imaging/selfcalibration         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

#include "ObitImageMosaic.h"
#include "ObitImageUtil.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitOTF.h"
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitIOOTFFITS.h"
#include "ObitDConCleanOTF.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitDisplay.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitThread.h"
#include "ObitFITS.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitTableOTFSoln.h"
#include "ObitTableOTFCal.h"
#include "ObitOTFSoln2Cal.h"
#include "ObitDConCleanWindow.h"
#include "ObitTableCCUtil.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* OTFSCalIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void OTFSCalOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Get PSF image */
ObitImage* getPSF (ObitInfoList *myInput, ObitErr *err);
/* Define output data */
ObitOTF* setOutputData (ObitInfoList *myInput, ObitOTF* inData, 
			ObitErr *err);
/* Define output dirty image */
ObitImage* setOutputDirty (ObitInfoList *myInput, ObitOTF *outData,
			   ObitErr *err);
/* Define output clean image */
ObitImage* setOutputClean (ObitInfoList *myInput, ObitImage *dirty,
			   ObitErr *err);
/* Define output weight image */
ObitImage* setOutputWeight (ObitInfoList *myInput, ObitImage *dirty,
			    ObitErr *err);
/* Initial calibration */
void InitCal (ObitInfoList* myInput, ObitOTF* OTFData, ObitErr* err);
/* self calibration loop */
gboolean doSelfCal (ObitInfoList* myInput, ObitOTF* outData, ObitImage* dirty, 
		    ObitImage* clean, ObitImage* wt, ObitErr* err);
/* Editing */
void doEdit (ObitInfoList* myInput, ObitOTF* outData, ObitImage* clean, 
	     ObitErr* err);
/* Image/Clean */
gboolean doImage (ObitInfoList* myInput, ObitOTF* outData, ObitImage* dirty, 
		 ObitImage* clean, ObitImage* wt, gboolean noResid, ObitErr* err);
/* Write history */
void OTFSCalHistory (ObitInfoList* myInput, ObitOTF* inData, 
		      ObitData* outImage, ObitErr* err);


/* Program globals */
gchar *pgmName = "OTFSCal";       /* Program name */
gchar *infile  = "OTFSCal.in" ;   /* File with program inputs */
gchar *outfile = "OTFSCal.out";   /* File to contain program outputs */
olong  pgmNumber;        /* Program number (like POPS no.) */
olong  AIPSuser;         /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;          /* Number of AIPS directories */
gchar **AIPSdirs=NULL;   /* List of AIPS data directories */
olong  nFITS=0;          /* Number of FITS directories */
gchar **FITSdirs=NULL;   /* List of FITS data directories */
ObitInfoList *myInput  = NULL;    /* Input parameter list */
ObitInfoList *myOutput = NULL;    /* Output parameter list */
olong  highCalVer=0;              /* Highest OTFCal version in initial Cal */
gboolean doBeam = TRUE;           /* Make dirty beam? */
ObitDisplay *display=NULL;        /* Image display */
ObitDConCleanOTF *myClean=NULL;   /* CLEAN process */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task to image a uv data set                                     */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  gboolean     quit;
  ObitSystem   *mySystem= NULL;
  ObitOTF      *inData = NULL, *outData = NULL;
  ObitImage    *dirty=NULL, *clean=NULL, *wt=NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = OTFSCalIn (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input data */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Calibrate/select to output data */
  outData = setOutputData (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Initialize calibration on outData */
  InitCal (myInput, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Index */
  ObitOTFUtilIndex  (outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create output images */
  dirty = setOutputDirty  (myInput, outData, err);
  clean = setOutputClean  (myInput, dirty, err);
  wt    = setOutputWeight (myInput, dirty, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Loop imaging and calibrating */
  quit = doSelfCal (myInput, outData, dirty, clean, wt, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  if (quit) goto done;

  /* Auto Edit */
 done:
  doEdit (myInput, outData, clean, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Final Image/Clean */
  Obit_log_error(err, OBIT_InfoErr, "    ");
  Obit_log_error(err, OBIT_InfoErr, "Making final image");
  ObitErrLog(err);
  quit = doImage (myInput, outData, dirty, clean, wt, FALSE, err);
  quit = ObitDisplayShow (display, (Obit*)clean, NULL, 1, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* History */
  OTFSCalHistory (myInput, inData, (ObitData*)outData, err);
  OTFSCalHistory (myInput, inData, (ObitData*)dirty,   err);
  OTFSCalHistory (myInput, inData, (ObitData*)clean,   err);
  OTFSCalHistory (myInput, inData, (ObitData*)wt,      err);
  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  dirty     = ObitUnref(dirty);
  clean     = ObitUnref(clean);
  wt        = ObitUnref(wt);
  myClean   = ObitDConCleanOTFUnref(myClean);

  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* OTFSCalIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "OTFSCalIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs/output InfoList */
  list = defaultInputs(err);
  myOutput = defaultOutputs(err);

  /* command line arguments */
  if (argc<=1) Usage(); /* must have arguments */
  /* parse command line */
  for (ax=1; ax<argc; ax++) {

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
      
     } else if (strcmp(arg, "-inOTF") == 0) { /*inOTF */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inOTF", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

     } else if (strcmp(arg, "-dispURL") == 0) { /* Display server URL */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "dispURL", OBIT_string, dim, strTemp);
      
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
} /* end OTFSCalIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: OTFSCal -input file -output ofile [args]\n");
    fprintf(stderr, "OTFSCal Obit task to image OTF data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def OTFSCal.in\n");
    fprintf(stderr, "  -output output result file, def OTFSCal.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -inOTF input FITS OTF file\n");
    fprintf(stderr, "  -inDisk input image disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output image FITS Image file\n");  
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -dispURL display server URL \n");
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
/*     inOTF     Str [?]    input FITS uv file name [no def]              */
/*     inDisk    Int        input  FITS OTF disk no  [def 1]              */
/*     outDisk   Int        output FITS image disk no  [def 1]            */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     Targets   Str (16,1) Targets selected, blank = all                 */
/*     Scans     Int (2)    Scans selected, all                           */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     CLEANBox  Int[4,?]   Clean box, def=all                            */
/*     autoWindow Boo(1)    If true set windows automatically, def=FALSE  */
/*     Gain      Flt (1)    Clean gain, def=0.1                           */
/*     minFlux   Flt (1)    Clean minimum flux density, def=0             */
/*     Niter     Int (1)    Maximum # of CLEAN comp., def=No CLEAN        */
/*     Patch     Int (1)    Clean Min. BEAM half-width, def=100           */
/*     BeamSize  Flt (1)    Clean beam (asec)                             */
/*     deMode    Boo (1)    subtract mode when forming image, def=FALSE   */
/*     dispURL   Str(48)    Display derver URL                            */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp, iarray[4];
  ofloat ftemp, farray[3];
  gboolean btemp;
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
  strTemp = "OTFData.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inOTF", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS Image file name */
  strTemp = "OTFSCalOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS Image disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* TARGETS selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Targets", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* SCANS */
  iarray[0] = 0; iarray[1] = 0;
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Scans", OBIT_string, dim, iarray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Fraction of peak per major cycle, def=0.75 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.75; 
  ObitInfoListPut (out, "fracPeak", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Apply calibration/selection?, def=True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doCalSelect", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Subtract mode when making image?, def=False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "deMode", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  >0 => apply gain calibration, 2=> cal. wt, def=no cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doCalib", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Gain table (Cal/Soln) table to apply, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "gainUse", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "flagVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /*  Clean box, def all */
  dim[0] = 4; dim[1] = 1;
  iarray[0] = iarray[1] = iarray[2] = iarray[3] = 0;
  ObitInfoListPut (out, "CLEANBox", OBIT_oint, dim, &iarray, err);

  /* autoWindow?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "autoWindow", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Clean gain, def = 0.1 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.1; 
  ObitInfoListPut (out, "Gain", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Clean minimum flux density, def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "minFlux", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum # of CLEAN comp., def = 0 (no clean) */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "Niter", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  

  /* Clean Min. BEAM half-width, def = 200 */
  dim[0] = 1;dim[1] = 1;
  itemp = 200; 
  ObitInfoListPut (out, "Patch", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Clean beam maj, min, PA (asec, asec, deg), def = 0 (fit) */
  dim[0] = 1;dim[1] = 1;
  farray[0] = 0.0; 
  ObitInfoListPut (out, "BeamSize", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum pixels searched in inner cycle, def = 50000 */
  dim[0] = 1;dim[1] = 1;
  itemp = 50000; 
  ObitInfoListPut (out, "maxPixel", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Display URL, def = "None" */
  strTemp = "None";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "dispRL", OBIT_string, dim, strTemp, err);
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
  oint         doCalib;
  gchar        *strTemp, *tname, tmpFile[129];
  ofloat       ftemp;
  gboolean doCalSelect;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Set temporary FITS image name as outName */
  /* output FITS file name */
  if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (tmpFile, strTemp, 128);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    tname = g_strconcat ("tmp", tmpFile, NULL);
    strncpy (tmpFile, tname, 128);
    g_free(tname);
  } else { 
    strncpy (tmpFile, "tmpImage.fits", 128);
  }
  dim[0] = MIN (128, strlen(tmpFile));
  ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, tmpFile);

  /* Convert "RACenter", "DecCenter" to "RA", "Dec" */
  ftemp = 0.0; dim[0] = dim[1] = 1;
  ObitInfoListGetTest(myInput, "RACenter",  &type, dim, &ftemp);
  ObitInfoListAlwaysPut(myInput, "RA",  OBIT_float, dim, &ftemp);
  ftemp = 0.0; dim[0] = dim[1] = 1;
  ObitInfoListGetTest(myInput, "DecCenter",  &type, dim, &ftemp);
  ObitInfoListAlwaysPut(myInput, "Dec",  OBIT_float, dim, &ftemp);

  /* Data Clip level, rename ClipData to Clip */
  ftemp = 1.0e20; dim[0] = dim[1] = 1;
  ObitInfoListGetTest(myInput, "ClipData",  &type, dim, &ftemp);
  ObitInfoListAlwaysPut(myInput, "Clip",  OBIT_float, dim, &ftemp);

  /* Make sure xCells negative */
  ftemp = 0.0;
  ObitInfoListGetTest(myInput, "xCells",  &type, dim, &ftemp);
  ftemp = -fabs(ftemp);
  ObitInfoListAlwaysPut(myInput, "xCells",  OBIT_float, dim, &ftemp);

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
/*       ObitOTF with input data                                          */
/*----------------------------------------------------------------------- */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitOTF       *inData = NULL;
  ObitTableOTFTarget* targetTable=NULL;
  ObitInfoType type;
  olong         disk, nrec, nThreads, ver, qual;
  gchar        *strTemp, inFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat       RACenter, DecCenter;
  odouble      RA, Dec;
  gchar        *targets, target[24];
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Targets", "Scans", "Feeds", "timeRange", "keepCal",
    "doCalSelect", "doCalib", "gainUse", "flagVer", 
    NULL
  };
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  inData = newObitOTF("input OTF data");
  
  /* input FITS file name */
  if (ObitInfoListGetP(myInput, "inOTF", &type, dim, (gpointer)&strTemp)) {
    strncpy (inFile, strTemp, 128);
  } else { 
    strncpy (inFile, "No_Filename_Given", 128);
  }
  
  /* input FITS disk */
  ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

  /* define object */
  nrec = 5000;
  nThreads = 1;
  ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
  nrec *= nThreads;
  ObitOTFSetFITS (inData, nrec, disk, inFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Ensure inData fully instantiated and OK */
  ObitOTFFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Get input select/cal parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, inData->name, inData);
 
  /* Default position */
  RACenter = 0.0;
  ObitInfoListGetTest(myInput, "RA",  &type, dim, &RACenter);
  DecCenter = 0.0;
  ObitInfoListGetTest(myInput, "Dec",  &type, dim, &DecCenter);
  if ((RACenter==0.0) && (DecCenter==0.0)) {
    if (ObitInfoListGetP(myInput, "Targets", &type, dim, (gpointer*)&targets)) {
      strncpy (target, targets, MIN(24,dim[0])); target[dim[0]]=0;
      if (target[0]!=' ') { /* Check if given */
	ver = 1;
	targetTable = 
	  newObitTableOTFTargetValue ("TargetTable", (ObitData*)inData, &ver, OBIT_IO_ReadOnly, err);
	qual = 0;
	ObitTableOTFTargetGetByName(targetTable, target, qual, &RA, &Dec, err);
	if (err->error) Obit_traceback_val (err, routine, inData->name, inData);
	targetTable = ObitTableOTFTargetUnref(targetTable);
	RACenter  = (ofloat)RA;
	DecCenter = (ofloat)Dec;
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "RA",  OBIT_float, dim, &RACenter);
	ObitInfoListAlwaysPut(myInput, "Dec",  OBIT_float, dim, &DecCenter);
      }
    }
  }

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get PSF Image                                                         */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       PSF Image                                                        */
/*----------------------------------------------------------------------- */
ObitImage* getPSF (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage *PSF=NULL;
  olong PSFDisk=0;
  gchar PSFFile[256];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "getPSF";

  if (err->error) return PSF;

  /* Instrument PSF */  
  if (ObitInfoListGetP(myInput, "PSFFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (PSFFile, strTemp, 255);
    ObitInfoListGetTest(myInput, "PSFDisk",  &type, dim, &PSFDisk);
    PSF = newObitImage("PSF Image");
    ObitImageSetFITS (PSF, OBIT_IO_byPlane, PSFDisk, PSFFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, PSF->name, PSF);
  } else {  /* No PSF */
    Obit_log_error(err, OBIT_Error,"%s: No PSF Image File specified", routine);
  }
  return PSF;
} /* end getPSF */


/*----------------------------------------------------------------------- */
/*  Get output data                                                       */
/*  Create file anc copy/calibrate selected data from inData              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    input data to copy from                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitOTF with input data                                          */
/*----------------------------------------------------------------------- */
ObitOTF* setOutputData (ObitInfoList *myInput, ObitOTF *inData, 
			ObitErr *err)
{
  ObitOTF      *outData = NULL;
  ObitInfoType type;
  olong        disk, nrec, nThreads;
  gchar        *strTemp, *fullFITS, outFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *imgparms[] = {  /* Imaging parameters */
    "RA", "Dec", "xCells", "yCells", 
    "nx", "ny", "minWt", "deMode", "deBias", "Proj", "ConvType", "ConvParm",
    "outName", "outDisk", "doScale", "doFilter", "Clip",
    NULL
  };
  gchar *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  outData = newObitOTF("output OTF data");
  
  /* input FITS file name */
  if (ObitInfoListGetP(myInput, "outOTF", &type, dim, (gpointer)&strTemp)) {
    strncpy (outFile, strTemp, 128);
  } else { 
    strncpy (outFile, "No_Filename_Given", 128);
  }
  
  /* output FITS disk */
  ObitInfoListGet(myInput, "outOTFdsk", &type, dim, &disk, err);

  /* If the output file already exists, delete it and recreate */
  fullFITS = ObitFITSFilename (disk, outFile, err);
  if (ObitFileExist (fullFITS, err)) 
    ObitFileZapFile(fullFITS, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  g_free(fullFITS);

  /* define object */
  nrec = 5000;
  nThreads = 1;
  ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
  nrec *= nThreads;
  ObitOTFSetFITS (outData, nrec, disk, outFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

   /* Copy from inData */  
  ObitOTFCopy (inData, outData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);

  /* Delete some tables */
  ObitOTFZapTable (outData, "OTFFlag", -1, err);
  ObitOTFZapTable (outData, "History",  1, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, outData);

  /* Set imaging parameters */
  ObitInfoListCopyList (myInput, outData->info, imgparms);
    
  return outData;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/*  Create output dirty image                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   Output OTF with imaging parameters set                  */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Returns image object                                                 */
/*----------------------------------------------------------------------- */
ObitImage* setOutputDirty (ObitInfoList *myInput, ObitOTF *outData,
			   ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  gchar     outFile[129], tmpFile[129], *tname, *outName, *strTemp, *Type=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  olong      Aseq, cno;
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong      i, disk;
  gboolean  exist;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputDirty";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  
  /* Create basic output Image Object */
  outImage = newObitImage("Output Dirty Image");
  /* Fill in details */
  outImage = ObitOTFUtilCreateImage (outData, err);
 
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */ 
    /* Output AIPS name */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    if ((outName!=NULL) && (strlen(outName)>0) && (outName[0]!=' ') && (outName[1]!=' ')) {
      strncpy (Aname, outName, 13); 
    } else { /* No output name given */
      strncpy (Aname, "no name", 13); 
    }
    Aname[12] = 0;

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Dirty ", 7);
    }
    if ((strlen(Aclass)<=0) || (Aclass[0]==' ') || (Aclass[1]==' '))
      strncpy (Aclass, "Dirty ", 7);
    Aclass[6] = 0;

    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);

    /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* Get name */
    ObitInfoListGet (myInput, "outFile", &type, dim, tmpFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    
    /* add .fits */
    tname = g_strconcat (tmpFile, "Dirty.fits", NULL);
    strncpy (outFile, tname, 128);
    g_free(tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    
    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } /* End FITS output image */

  /* Instantiate */ 
   ObitImageFullInstantiate (outImage, FALSE, err);
 if (err->error) Obit_traceback_val (err, routine, outData->name, outImage);

 return outImage;
} /* end setOutputDirty  */

/*----------------------------------------------------------------------- */
/*  Create output CLEAN image                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      dirty     Output dirty image                                      */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Returns image object                                                 */
/*----------------------------------------------------------------------- */
ObitImage* setOutputClean (ObitInfoList *myInput, ObitImage *dirty,
			   ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  gchar     outFile[129], tmpFile[129], *tname, *outName, *strTemp, *Type=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  olong      Aseq, cno;
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong      i, disk;
  gboolean  exist;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputClean";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  
  /* Create basic output Image Object */
  outImage = newObitImage("Output Clean Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */ 
    /* Output AIPS name */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    if ((outName!=NULL) && (strlen(outName)>0) && (outName[0]!=' ') && (outName[1]!=' ')) {
      strncpy (Aname, outName, 13); 
    } else { /* No output name given */
      strncpy (Aname, "noname", 13); 
    }
    Aname[12] = 0;

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Clean", 7);
    }
    if ((strlen(Aclass)<=0) || (Aclass[0]==' ') || (Aclass[1]==' '))
      strncpy (Aclass, "Clean", 7);
    Aclass[6] = 0;

    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);

    /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* Get name */
    ObitInfoListGet (myInput, "outFile", &type, dim, tmpFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    
    /* add .fits */
    tname = g_strconcat (tmpFile, "Clean.fits", NULL);
    strncpy (outFile, tname, 128);
    g_free(tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    
    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } /* End FITS output image */

  /* Clone from dirty */
  ObitImageClone (dirty, outImage, err);
  if (err->error) Obit_traceback_val (err, routine, dirty->name, outImage);
 
 return outImage;
} /* end setOutputClean  */

/*----------------------------------------------------------------------- */
/*  Create output weight image                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Returns image object                                                 */
/*----------------------------------------------------------------------- */
ObitImage* setOutputWeight (ObitInfoList *myInput, ObitImage *dirty,
			    ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  gchar     outFile[129], tmpFile[129], *tname, *outName, *strTemp, *Type=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  olong     Aseq, cno;
  olong     blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong     trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong     i, disk;
  gboolean  exist;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputWeight";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  
  /* Create basic output Image Object */
  outImage = newObitImage("Output Weight Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */ 
    /* Output AIPS name */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    if ((outName!=NULL) && (strlen(outName)>0) && (outName[0]!=' ') && (outName[1]!=' ')) {
      strncpy (Aname, outName, 13); 
    } else { /* No output name given */
      strncpy (Aname, "noname", 13); 
    }
    Aname[12] = 0;

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "out2Class", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Wt  ", 7);
    }
    if ((strlen(Aclass)<=0) || (Aclass[0]==' ') || (Aclass[1]==' '))
      strncpy (Aclass, "Wt  ", 7);
    Aclass[6] = 0;

    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);

    /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* Get name */
    ObitInfoListGet (myInput, "outFile", &type, dim, tmpFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    
    /* add .fits */
    tname = g_strconcat (tmpFile, "Wt.fits", NULL);
    strncpy (outFile, tname, 128);
    g_free(tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    
    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } /* End FITS output weight image */
 
  /* Clone from dirty */
  ObitImageClone (dirty, outImage, err);
  if (err->error) Obit_traceback_val (err, routine, dirty->name, outImage);
 
 return outImage;
} /* end setOutputWeight  */

/*----------------------------------------------------------------------- */
/*  Do initial calibration                                                */
/*   1) Create initial dummy OTFCal table at 0.25 on min solution         */
/*   2) Prior calibration if requested                                    */
/*   3) Common mode (no prior) calibration if requested                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      OTFData   ObitOTF to calibrate                                    */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void InitCal (ObitInfoList* myInput, ObitOTF* OTFData, ObitErr* err)
{
  ObitTableOTFCal  *CalTab=NULL;
  ObitTableOTFSoln *SolnTab=NULL;
  ObitImage        *prior=NULL, *PSF=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       solnInt, priorInt, commonInt, Int, maxResFlx, minResFlx;
  ofloat       clipSig;
  ofloat       maxInt, maxIntFact = 10.0;  /* scale factor for maxInt from solInt */
  olong        itemp, priorDisk=0, ncoef=1, solnVer;
  gboolean     btemp, doPrior, priorMod=FALSE, doOffset=TRUE;
  gchar        *strTemp, priorFile[256], tstring[100];
  gchar        *calParms[] = {  /* Parameters to calibrate data */
    "minResFlx", "maxResFlux", "minEl", "Clip", "doScale",
    NULL };
  gchar        *routine = "InitCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitOTFIsA(OTFData));

  /* Control parameters */
  ObitInfoListGetTest(myInput,"doOffset",   &type, dim, &doOffset);
  ObitInfoListGetTest(myInput,"offsetFact", &type, dim, &maxIntFact);
  minResFlx = 0.0;
  ObitInfoListGetTest(myInput, "minResFlx",  &type, dim, &minResFlx);
  maxResFlx = 0.0;
  ObitInfoListGetTest(myInput, "maxResFlx",  &type, dim, &maxResFlx);
  solnInt = 1.0;
  ObitInfoListGetTest(myInput, "solnInt",  &type, dim, &solnInt);
  commonInt = -1.0;
  ObitInfoListGetTest(myInput, "commonInt", &type, dim, &commonInt);
  priorInt = 1.0;
  ObitInfoListGetTest(myInput, "priorInt",  &type, dim, &priorInt);
  if (ObitInfoListGetP(myInput, "priorFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (priorFile, strTemp, 255);
    /* Do prior if priorFile neither "None", nor "    " */
    doPrior = strncmp(priorFile,"None",4!=0) &&
      strncmp(priorFile,"    ",4!=0);
    ObitInfoListGetTest(myInput, "priorMod",   &type, dim, &priorMod);
    ObitInfoListGetTest(myInput, "priorDisk",  &type, dim, &priorDisk);
  } else {  /* No prior */
    doPrior = FALSE;
  }
  if (!doPrior) priorInt = 1.0e20;

  /* Instrument PSF */  
  PSF = getPSF (myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, OTFData->name);

  /* Dummy OTFCal table - need min solution interval */
  Int = MIN (solnInt, priorInt);
  if (commonInt>0.0) Int = MIN (Int, commonInt);
  /* Int *= 0.25;   Increment for Cal Table */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Make Dummy Cal table 1 with interval %6.2f sec", Int);
  ObitErrLog(err); 
  dim[0] = dim[1] = dim[2] = 1;
  Int /= 86400.0;  /* to Days */
  ObitInfoListAlwaysPut(OTFData->info, "solInt",  OBIT_float, dim, &Int);
  CalTab = ObitOTFGetDummyCal (OTFData, OTFData, 1, ncoef, err);
  if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
  highCalVer = CalTab->tabVer;
  CalTab = ObitTableOTFCalUnref(CalTab);

  /* Get input cal parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, OTFData->info, calParms);

  /* In solutions clip residuals at 10 sigma */
  clipSig = 10.0;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (OTFData->info, "clipSig", OBIT_float, dim, &clipSig);
  
  /* Self cal loop */
  /* Common mode Cal with non zero prior? */
  if (doPrior) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Common mode calibration using prior and %6.2f sec", priorInt);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Prior image %d/%s", priorDisk,priorFile);
    /* Set control parameters */
    if (doOffset) {
      maxInt = maxIntFact*priorInt/86400.0;
      ObitInfoListAlwaysPut (OTFData->info, "maxInt", OBIT_float, dim, &maxInt);
      strcpy (tstring, "Both");
       Obit_log_error(err, OBIT_InfoErr, "Detector offsets on timescale min. %6.1f sec", 
		      maxIntFact*priorInt);
    } else {
      strcpy (tstring, "Common");
    }
    dim[0] = strlen (tstring);
    ObitInfoListAlwaysPut (OTFData->info, "calType", OBIT_string, dim, tstring);
    priorInt /= 86400.0; /* to days */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (OTFData->info, "solInt", OBIT_float, dim, &priorInt);
    ObitInfoListAlwaysPut (OTFData->info, "doCC", OBIT_float, dim, &btemp);

    /* Create prior Image Object */
    prior = newObitImage("Prior Image");
    ObitImageSetFITS (prior, OBIT_IO_byPlane, priorDisk, priorFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
    SolnTab = ObitOTFUtilResidCal (OTFData, OTFData, prior, priorMod,
				   PSF, err);
    if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
    solnVer = SolnTab->tabVer;
    SolnTab = ObitTableOTFSolnUnref(SolnTab);

    /* Apply to OTFCal table */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(OTFData->info, "solnUse", OBIT_long, dim, &solnVer);
    ObitInfoListAlwaysPut(OTFData->info, "calIn",   OBIT_long, dim, &highCalVer);
    itemp = highCalVer+1;
    ObitInfoListAlwaysPut(OTFData->info, "calOut",  OBIT_long, dim, &itemp);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Apply Soln %d to Cal %d, write %d", solnVer, highCalVer,itemp);
    ObitErrLog(err); 

    CalTab = ObitOTFSoln2Cal (OTFData, OTFData, err);
    if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
    highCalVer = CalTab->tabVer;
    CalTab = ObitTableOTFCalUnref(CalTab);
  }

  /* Common mode with zero prior? */
  if (commonInt>0.0) {
    /* Set control parameters */
    if (doOffset) {
      maxInt = maxIntFact*commonInt/86400.0;
      ObitInfoListAlwaysPut (OTFData->info, "maxInt", OBIT_float, dim, &maxInt);
      strcpy (tstring, "Both");
       Obit_log_error(err, OBIT_InfoErr, "Detector offsets on timescale min. %6.1f sec", 
		      maxIntFact*commonInt);
    } else {
      strcpy (tstring, "Common");
    }
    dim[0] = strlen (tstring);
    ObitInfoListAlwaysPut (OTFData->info, "calType", OBIT_string, dim, tstring);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Common mode calibration at %6.2f sec", commonInt);

    commonInt /= 86400.0;  /* to days */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (OTFData->info, "solInt",  OBIT_float,  dim, &commonInt);

    SolnTab = ObitOTFUtilResidCal (OTFData, OTFData, NULL, FALSE,
				   PSF, err);
    if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
    solnVer = SolnTab->tabVer;
    SolnTab = ObitTableOTFSolnUnref(SolnTab);

    /* Apply to OTFCal table */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(OTFData->info, "solnUse", OBIT_long, dim, &solnVer);
    ObitInfoListAlwaysPut(OTFData->info, "calIn",   OBIT_long, dim, &highCalVer);
    itemp = highCalVer+1;
    ObitInfoListAlwaysPut(OTFData->info, "calOut",  OBIT_long, dim, &itemp);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Apply Soln %d to Cal %d, write %d", solnVer, highCalVer,itemp);
    ObitErrLog(err); 
    CalTab = ObitOTFSoln2Cal (OTFData, OTFData, err);
    if (err->error) Obit_traceback_msg (err, routine, OTFData->name);
    highCalVer = CalTab->tabVer;
    CalTab = ObitTableOTFCalUnref(CalTab);
  }
}  /* end InitCal */

/*----------------------------------------------------------------------- */
/*  self calibration loop                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitOTF to image                                        */
/*      dirty     output dirty image                                      */
/*      clean     output clean image                                      */
/*      wt        output weight image                                     */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*      TRUE if user requests quit                                        */
/*----------------------------------------------------------------------- */
gboolean doSelfCal (ObitInfoList* myInput, ObitOTF* outData, ObitImage* dirty, 
		    ObitImage* clean, ObitImage* wt, ObitErr* err)
{
  gboolean quit=FALSE;
  ObitImage *PSF=NULL;
  ObitTableOTFSoln *SolnTab=NULL;
  ObitTableOTFCal  *CalTab=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat solnInt, solnMult[100], si, maxInt, clipSig;
  ofloat maxIntFact = 10.0;  /* scale factor for maxInt from solInt */
  gboolean btemp, doOffset=TRUE;
  olong n, i, iloop, nloop, solnVer, itemp;
  gchar calType[24];
  gchar *routine="doSelfCal";

  if (err->error) return quit;  /* prior error */

  /* Control info */
  ObitInfoListGetTest(myInput,"doOffset",   &type, dim, &doOffset);
  ObitInfoListGetTest(myInput,"offsetFact", &type, dim, &maxIntFact);
  ObitInfoListGet(myInput,"solnInt",  &type, dim, &solnInt, err);
  ObitInfoListGet(myInput,"solnMult", &type, dim, solnMult, err);
  /* Number of loops? */
  n = MIN(100,dim[0]);
  nloop = 0;
  for (i=0; i<n; i++) {
    if (solnMult[i]<=0.0) break;
    nloop = i+1;
  }

  /* Get PSF */
  PSF = getPSF (myInput, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
  
  /* Be sure to apply last calibration from here on */
  dim[0] = dim[1] = dim[2] = 1;
  btemp = TRUE;
  ObitInfoListAlwaysPut (outData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = 2; 
  ObitInfoListAlwaysPut (outData->info, "doCalib", OBIT_oint, dim, &itemp);
  itemp = 0; 
  ObitInfoListAlwaysPut (outData->info, "gainUse", OBIT_oint, dim, &itemp);
  
  /* Set solution type */
  if (doOffset) strcpy (calType, "Both");
  else strcpy (calType, "Common");
  dim[0] = strlen(calType); dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outData->info, "calType", OBIT_string, dim, calType);

  /* In solutions clip residuals at 10 sigma */
  clipSig = 10.0;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outData->info, "clipSig", OBIT_float, dim, &clipSig);
  
  /* Self cal loop */
  for (iloop=0; iloop<nloop; iloop++) {

    /* Set solution interval */
    si = solnInt * solnMult[iloop];
    Obit_log_error(err, OBIT_InfoErr, "     ");
    Obit_log_error(err, OBIT_InfoErr, 
		   "*** Selfcal loop %d of %d Solution interval %6.1f sec, type %s",
		   iloop+1, nloop, si, calType);
    if (doOffset) 
      Obit_log_error(err, OBIT_InfoErr, "Detector offsets on timescale min. %6.1f sec", 
		     maxIntFact*si);
    ObitErrLog(err); 

    /* Make/CLEAN image */
    quit = doImage (myInput, outData, dirty, clean, wt, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
    if (quit) return quit;  /* Bail out? */
    
    /* Solution interval */
    si /= 86400.0; /* to days */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (outData->info, "solInt", OBIT_float, dim, &si);
    /* Interval for offset calculation */
    if (doOffset) {
      maxInt = maxIntFact*si;
      ObitInfoListAlwaysPut (outData->info, "maxInt", OBIT_float, dim, &maxInt);
   }

    /* Do solution wrt Cal version highCalVer */
    ObitInfoListAlwaysPut (outData->info, "gainUse", OBIT_oint, dim, &highCalVer);
    SolnTab = ObitOTFUtilResidCal (outData, outData, clean, TRUE,
				   PSF, err);
    if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
    solnVer = SolnTab->tabVer;
    SolnTab = ObitTableOTFSolnUnref(SolnTab);
    
    /* Delete OTFCal table to be written */
    ObitOTFZapTable (outData, "OTFCal", highCalVer+1, err);
    if (err->error) Obit_traceback_val (err, routine, outData->name, quit);

    /* Apply to OTFCal table */
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(outData->info, "solnUse", OBIT_long, dim, &solnVer);
    ObitInfoListAlwaysPut(outData->info, "calIn",   OBIT_long, dim, &highCalVer);
    itemp = highCalVer+1;
    ObitInfoListAlwaysPut(outData->info, "calOut",  OBIT_long, dim, &itemp);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Apply Soln %d to Cal %d, write %d", solnVer, highCalVer,itemp);
    ObitErrLog(err); 
    CalTab = ObitOTFSoln2Cal (outData, outData, err);
    if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
    CalTab = ObitTableOTFCalUnref(CalTab);

    /* Delete SolnTab */
    ObitOTFZapTable (outData, "OTFSoln", solnVer, err);

    /* Most recent table */
    itemp = 0; 
    ObitInfoListAlwaysPut (outData->info, "gainUse", OBIT_oint, dim, &itemp);
  } /* End self cal loop */

  /* Cleanup */
  PSF = ObitImageUnref(PSF);

  return quit;
} /* end doSelfCal */

/*----------------------------------------------------------------------- */
/*  Edit data                                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitOTF to image                                        */
/*      clean     output clean image                                      */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doEdit (ObitInfoList* myInput, ObitOTF* outData, ObitImage* clean, 
	     ObitErr* err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean doedit, btemp;
  olong itemp, FGVer;
  ofloat flagInt;
  gchar *editParm[]= {
    "maxRMS", "maxRatio",
    NULL};
  gchar *routine = "doEdit";

  if (err->error) return;  /* prior error */

  /* Anything to do? */
  ObitInfoListGetTest(myInput, "doEdit",  &type, dim, &doedit);
  if (!doedit) return;

  Obit_log_error(err, OBIT_InfoErr, "Editing data");
  ObitErrLog(err);
 
  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, outData->info, editParm);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  dim[0] = dim[1] = dim[2] = 1;
  btemp = TRUE;
  ObitInfoListAlwaysPut (outData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = 1;
  ObitInfoListAlwaysPut (outData->info, "doCalib", OBIT_long, dim, &itemp);
  itemp = 0;
  ObitInfoListAlwaysPut (outData->info, "gainUse", OBIT_long, dim, &itemp);
  FGVer = 1; dim[0] = 1; dim[1] = 1; 
  ObitInfoListGetTest(myInput, "flagver",  &type, dim, &FGVer);
  ObitInfoListAlwaysPut (outData->info, "FGVer", OBIT_long, dim, &FGVer);
  ObitInfoListGetTest(myInput, "flagInt",  &type, dim, &flagInt);
  flagInt /= 86400.0; /* to days */
  ObitInfoListAlwaysPut (outData->info, "flagInt", OBIT_float, dim, &flagInt);

  /* Flag */
  /* NOTE: THIS USES THE CLEAN IMAGE PIXELS NOT THE CC MODEL */
  ObitOTFGetSolnFlag (outData, clean, outData, FGVer, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

} /* end doEdit */

/*----------------------------------------------------------------------- */
/*  Make and Hogbom clean image                                           */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outData   ObitOTF to image                                        */
/*      dirty     output dirty image                                      */
/*      clean     output clean image                                      */
/*      wt        output weight image                                     */
/*      noResid   If TRUE don't include residuals                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   tReturn:                                                             */
/*      TRUE if user requests Quit                                        */
/*----------------------------------------------------------------------- */
gboolean doImage (ObitInfoList* myInput, ObitOTF* outData, ObitImage* dirty, 
		  ObitImage* clean, ObitImage* wt, gboolean noResid,
		  ObitErr* err)
{
  gboolean quit=FALSE;
  ObitImage *PSF=NULL;
  ObitTableCC *CCTab=NULL;
  ofloat BeamSize, Clip;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean btemp;
  olong i, ver, noParms;
  gchar *dispURL, tname[256];
  gchar *CLEANParms[] = {  /* Clean parameters */
    "CLEANBox", "autoWindow", "Gain", "minFlux", "Niter", "Patch", 
    "fracPeak", "doScale", "doRestore", "dispURL", "ClipData",
    NULL
  };
  gchar *routine = "doImage";

  if (err->error) return quit;

  /* Instrument PSF */
  if (doBeam) {
    PSF = getPSF (myInput, err);
    if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
  }

  /* Data clipping? */
  Clip = 1.0e20;
  ObitInfoListGetTest(myInput, "ClipData",  &type, dim, &Clip);
  if (fabs(Clip)<9999.0) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Imaging clipping data at +/- %12.4g Jy", Clip);
    ObitErrLog(err); 
  }

  /* Want PSF convolved with convolving kernal */
  dim[0] = dim[1] = dim[2] = 1;
  btemp = TRUE;
  ObitInfoListAlwaysPut (outData->info, "doConvBeam", OBIT_bool, dim, &btemp);

  /* Make image */
  ObitOTFUtilMakeImage (outData, dirty, doBeam, PSF, wt, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, quit);

  /* Only need to do dirty beam once 
  doBeam = FALSE; need this */

  /* Need Image display? */
  if (!display) {
    ObitInfoListGetP(myInput, "dispURL", &type, dim, (gpointer)&dispURL);
    /* dispURL not NULL terminated */
    if (dispURL) {for (i=0; i<dim[0]; i++) tname[i] = dispURL[i]; tname[i]=0;}
    if (dispURL && (strncmp(tname, "None", 4))) 
      display = ObitDisplayCreate("Display", tname, err);
    if (err->error) Obit_traceback_val (err, routine, clean->name, quit);
  }
  /* Copy dirty to clean */
  clean = ObitImageCopy (dirty, clean, err);

  /* CLEAN */
  if (myClean==NULL) {
    myClean = ObitDConCleanOTFCreate ("Clean", dirty, (ObitImage*)dirty->myBeam,
				      clean, err);
    /* Get input Clean parameters from myInput, copy to myClean */
    ObitInfoListCopyList (myInput, myClean->info, CLEANParms);
    BeamSize = 8.0; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListGetTest(myInput, "BeamSize", &type, dim, &BeamSize);
    if (BeamSize<=0.0) {  /* Default from data */
      BeamSize = outData->myDesc->beamSize*3600.0;
      ObitInfoListAlwaysPut(myInput, "BeamSize", OBIT_float, dim, &BeamSize);
    }
    BeamSize /= 3600.0;  /* to Deg */
    ObitInfoListAlwaysPut(myClean->info, "BeamSize", OBIT_float, dim, &BeamSize);
    /* Set window */
    ObitDConCleanDefWindow((ObitDConClean*)myClean, err);
    if (err->error) Obit_traceback_val (err, routine, clean->name, quit);
  }

  if (err->error) Obit_traceback_val (err, routine, outData->name, quit);

  /* Edit CLEAN box if necessary */
  quit = ObitDisplayShow (display, (Obit*)dirty, myClean->window, 1, err);
  if (err->error) Obit_traceback_val (err, routine, clean->name, quit);

  if (quit) return quit;  /* Bail out? */

  /* Want residuals? */
  dim[0] = dim[1] = dim[2] = 1;
  if (noResid) {
    btemp = TRUE;
    ObitInfoListAlwaysPut(myClean->info, "noResid", OBIT_bool, dim, &btemp);
  } else {
    btemp = FALSE;
    ObitInfoListAlwaysPut(myClean->info, "noResid", OBIT_bool, dim, &btemp);
  }
  /* Don't scale CCs */
  btemp = FALSE;
  ObitInfoListAlwaysPut(myClean->info, "doScaleCC", OBIT_bool, dim, &btemp);

  /* Do it */
  ObitDConCleanOTFDeconvolve ((ObitDCon*)myClean, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, quit);

  /* Compress CC table */
  ver     = 1;
  noParms = 0;
  CCTab   = newObitTableCCValue ("CC Table", (ObitData*)clean, &ver, 
				 OBIT_IO_ReadWrite, noParms, err);
  ObitTableCCUtilMerge (CCTab, CCTab, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, quit);
  CCTab = ObitTableCCUnref(CCTab);

  /* Cleanup */
  PSF = ObitImageUnref(PSF);

  return quit;
} /* end doImage */

/*----------------------------------------------------------------------- */
/*  Write History for OTFSCal                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to copy history from                            */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void OTFSCalHistory (ObitInfoList* myInput, ObitOTF* inData, 
		      ObitData* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "inOTF",  "inDisk", "PSFFile", "PSFDisk",
    "outDType", "outOTF", "outOTFdsk",
    "outName", "outClass", "out2Class", "outSeq", "outFile",  "outDisk", 
    "Targets", "Scans", "Feeds", "timeRange",  "keepCal", "minWt", "deMode", 
    "deBias", "doScale", "doFilter", "ClipData",
    "doCalSelect",  "doCalib",  "gainUse", "flagVer", 
    "CLEANBox",  "Gain",  "minFlux",  "Niter", "Patch",
    "BeamSize",  "fracPeak", "noResid", "doRestore", "autoWindow", 
    "priorFile", "priorDisk", "priorMod", "priorInt", "commonInt",
    "solnInt", "solnMult", "minResFlx", "maxResFlx", "doOffset", "offsetFact",
    "doEdit", "flagver", "flagInt", "maxRMS", "maxRatio",
    "nThreads",
    NULL};
  gchar *routine = "OTFSCalHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitOTFIsA(inData));
  g_assert (ObitDataIsA(outImage));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

  /* FITS - copy header */
  ObitHistoryCopyHeader (inHistory, outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end OTFSCalHistory  */

