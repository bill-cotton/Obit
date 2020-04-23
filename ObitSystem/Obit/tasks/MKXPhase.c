/* To do
 */

/* $Id$  */
/* MeerKAT X/Y phase bandpass calibration from noise diode            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2019                                               */
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
#include "ObitUVGSolve.h"
#include "ObitSkyModel.h"
#include "ObitTableSUUtil.h"
#include "ObitUVUtil.h"
#include "ObitTableSN.h"
#include "ObitTableCL.h"
#include "ObitTableCLUtil.h"
#include "ObitTableBP.h"
#include "ObitUVSoln2Cal.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* MKXPhaseIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void MKXPhaseOut (ObitInfoList* outList, ObitErr *err);
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
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Do channel solutions */
void XYBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		   ObitErr* err);
/* Write history */
void MKXPhaseHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Time averaging */
ObitUV* TimeAverage (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Create  BP table */
ObitTableBP* MakeBPTable (ObitUV* inData,  olong bchan, olong echan, 
			  olong nif, ObitErr *err);
/* Dump current phases to BP table */
void DumpBP(ObitTableBP* BPTable, ofloat *phase,  olong souId, ofloat time, 
	    olong bchan, olong echan, olong nif, olong numAnt, olong refAnt, ObitErr* err);
/* Update  BP table */
ObitTableBP* UpdateBPTable (ObitUV* inData, olong BPver, ObitUV* outData, 
			    ofloat time, olong bchan, olong echan, olong nif, 
			    ofloat *phase, ofloat RLPhase, ofloat RM, olong refAnt,
			    ObitErr *err);


/* Program globals */
gchar *pgmName = "MKXPhase";      /* Program name */
gchar *infile  = "MKXPhase.in" ;  /* File with program inputs */
gchar *outfile = "MKXPhase.out";  /* File to contain program outputs */
olong  pgmNumber;               /* Program number (like POPS no.) */
olong  AIPSuser;                /* AIPS user number number  */
olong  nAIPS=0;                 /* Number of AIPS directories */
gchar **AIPSdirs=NULL;          /* List of AIPS data directories */
olong  nFITS=0;                 /* Number of FITS directories */
gchar **FITSdirs=NULL;          /* List of FITS data directories */
ObitInfoList *myInput  = NULL;  /* Input parameter list */
ObitInfoList *myOutput = NULL;  /* Output parameter list */
olong  souNo=0;                 /* Single source number */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/* MeerKAT X-Y phase bandpass calibration from noise diode in XY autocor. */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL, *avgData=NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = MKXPhaseIn (argc, argv, err);
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

  /* Get output uvdata (for AIPS FG table) */
  outData = setOutputUV (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Save first source id in case it's the only one */
  if (inData->mySel->sources) souNo = inData->mySel->sources[0];
  
  /* Average data to SolInt */
  avgData =  TimeAverage(myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Do channel solutions, convert to BP table */
  XYBandpassCal(myInput, avgData, outData, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Write history */
  MKXPhaseHistory (myInput, outData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  avgData   = ObitUnref(avgData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* MKXPhaseIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "MKXPhaseIn";

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
} /* end MKXPhaseIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: MKXPhase -input file -output ofile [args]\n");
    fprintf(stderr, "MKXPhase Obit task determine bandpass for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def MKXPhase.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def MKXPhase.out\n");
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
  strTemp = "MKXPhase.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "MKXPhaseName";
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
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        itemp;
  gboolean     doCalSelect;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Channels */
  itemp = 1;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "BChan", OBIT_long, dim, &itemp);
  itemp = 0;
  ObitInfoListAlwaysPut (myInput, "EChan", OBIT_long, dim, &itemp);

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
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "timeRange", "BChan", "EChan", "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "refAnt",
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
/*  Create output uv data, defaults to input                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputUV (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, cno;
  oint      disk, Aseq;
  gchar     *Type, *strTemp=NULL;
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     tname[129], outFile[129];
  gchar     *routine = "setOutputUV";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */

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

    /* outClass given? */
    ObitInfoListGetP (myInput, "outClass", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "      ", 6)))
      ObitInfoListGetP (myInput, "inClass", &type, dim, (gpointer)&strTemp);
    for (i=0; i<6; i++) Aclass[i] = ' ';  Aclass[i] = 0;
    for (i=0; i<MIN(6,dim[0]); i++) Aclass[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 6;
    ObitInfoListAlwaysPut (myInput, "outClass", OBIT_string, dim, Aclass);

    /* outSeq given? */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    if (Aseq<=0) 
      ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outSeq", OBIT_oint, dim, &Aseq);

    /* outDisk given? */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) 
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* define object */
    nvis = 1;
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

    /* outDisk given? */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0) 
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    /* Save any defaulting on myInput */
    dim[0] = 1;
    ObitInfoListAlwaysPut (myInput, "outDisk", OBIT_oint, dim, &disk);

    /* define object */
    nvis = 1;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  /* Ensure outData fully instantiated and OK */
  ObitUVFullInstantiate (outUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

 return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Write History for MKXPhase                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void MKXPhaseHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "outFile",  "outDisk", "outName", "outClass", "outSeq", 
    "Sources","timeRange",  "BChan", "EChan",  "ChWid", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "Antennas",  "refAnt", "BPSoln", "nThreads",
   NULL};
  gchar *routine = "MKXPhaseHistory";

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
 
} /* end MKXPhaseHistory  */

/**
 * Apply calibration to data and averages to solint and returns averaged data.,
 * \param myInput Input parameters on InfoList    
 * \param inUV   Data to be used to determine calibration
 * \param err    ObitErr stack for reporting problems.
 * \return time averaged visibility data
 */
ObitUV* TimeAverage (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitUV       *avgData=NULL;
  ofloat       solInt;
  olong        corrType=2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *FQInclude[] = {"AIPS FQ", "AIPS AN", NULL};
  gchar *routine = "TimeAverage";

  /* error checks */
  if (err->error) return avgData;
  g_assert (ObitUVIsA(inData));
  
  /* Average data to solInt; 0=> all */
  solInt = 0.0;
  ObitInfoListGetTest (myInput, "solInt", &type, dim, &solInt); 
  if (solInt<=1.0e-5) solInt = 1000.0; 
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(inData->info, "timeAvg", OBIT_float, dim, &solInt);
  /* Only want autocorrelations */
  ObitInfoListAlwaysPut(inData->info, "corrType", OBIT_long, dim, &corrType);
  avgData = ObitUVUtilAvgT(inData, TRUE, avgData, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, avgData);
  ObitInfoListAlwaysPut(avgData->info, "corrType", OBIT_long, dim, &corrType);

  /* Be sure FQ table copied */
  ObitUVCopyTables (inData, avgData, NULL, FQInclude, err);

  return avgData;
} /* end TimeAverage */

/**
 * Determine Bandpass table, unit amplitude, phase = -XY autocorrelation phase 
 * \param myInput Input parameters on InfoList    
 * \param avgUV   Averaged data to be used to determine calibration
 * \param outUV   UV data onto which the BP table to be written
 * \param err     ObitErr stack for reporting problems.
 */
void  XYBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* outUV, 
		    ObitErr* err)
{
  ObitTableBP *BPTable=NULL, *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  olong ichan, nchan, bchan=0, echan=0, numAnt, chWid, iant;
  olong ivis, nvis, ifreq, nfreq, iif, nstok, nif, indx, nc;
  olong itemp, nchanIF, doBand, BPVer, outVer, BPSoln, refAnt;
  olong ia1, ia2, ver, suba, souId, hiLim, loLim, iVis, nVis;
  gboolean btemp, haveSome=FALSE;
  ofloat solInt=0.0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *vis, *ifvis, *fvis, *afvis, *rec, time=-1.0e20, etime;
  ofloat sumr, sumi, sumwt, wt;
  ObitUVDesc *desc=avgData->myDesc;
  gfloat *phase=NULL, fblank = ObitMagicF();
  gchar *Stokes="XY  ";
  gchar *routine = "XYBandpassCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(outUV));
  g_assert (ObitUVIsA(avgData));

  /* Any data in avgData? */
  if (avgData->myDesc->nvis<=0) {
    Obit_log_error(err, OBIT_Error, "NO averaged/calibrated data");
    return;
  }

  ObitUVGetSubA (avgData, err);
  numAnt  = avgData->myDesc->numAnt[0];/* actually highest antenna number */
  phase = g_malloc0((numAnt+10)*sizeof(ofloat));
  
  /* Range of channels -  really ALL */
  nchan = avgData->myDesc->inaxes[avgData->myDesc->jlocf];
  ObitInfoListGetTest (myInput, "BChan", &type, dim, &bchan); 
  if (bchan<=0) bchan = 1;
  ObitInfoListGetTest (myInput, "EChan", &type, dim, &echan); 
  if (echan<=0) echan = nchan;
  echan = MIN (echan, nchan);
  /* How many channels to average */
  ObitInfoListGetTest (myInput, "ChWid", &type, dim, &chWid); 
  chWid /= 2;
  chWid = MAX (0, chWid);
  
  /* Solution interval */
  solInt = 600.0; /* default 600 sec */
  ObitInfoListGetTest (myInput, "solInt", &type, dim, &solInt); 
  if (solInt<=0.0) solInt = 600.0; 
  solInt /= 86400.0; /* to days */

  /* Prior bandpass calibration */
  ObitInfoListGetTest (myInput, "doBand", &type, dim, &doBand); 
  ObitInfoListGetTest (myInput, "BPVer",  &type, dim, &BPVer); 
  ObitInfoListGetTest (myInput, "BPSoln", &type, dim, &BPSoln);
  if (BPVer<=0) BPVer = ObitTableListGetHigh (outUV->tableList, "AIPS BP");
  ObitInfoListGetTest (myInput, "refAnt", &type, dim, &refAnt); 
  
  /* Number of IFs */
  if (avgData->myDesc->jlocif>=0)  
    nif = avgData->myDesc->inaxes[avgData->myDesc->jlocif];
  else
    nif = 1;
  
  /* Number of channels, Stokes  */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  nstok = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];

  /* Accumulator array */
  nc      =  (echan-bchan+1);
  nchanIF =  nif * nc;
  phase = g_malloc0(nchanIF*(numAnt+10)*sizeof(ofloat));
		     
  /* Don't need to apply calibration */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(avgData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = -1;
  ObitInfoListAlwaysPut(avgData->info, "doCalib", OBIT_long, dim, &itemp);
  /* Only want XY */
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut(avgData->info, "Stokes", OBIT_string, dim, Stokes);
  
  /* Initialize as blanked */
  for (iant=0; iant<numAnt; iant++) {
    for (iif=0; iif<nif; iif++) {
      for (ichan=bchan-1; ichan<echan; ichan++) {
	indx = iant*nchanIF + iif*nc + ichan;
	phase[indx] = fblank;
      }
    }
  }
 
  /* Open */
  retCode = ObitUVOpen (avgData, OBIT_IO_ReadCal, err);
  if (err->error) goto cleanup;
  desc  = avgData->myDesc;

  /* Create temporary BP table */
  BPTable = MakeBPTable (avgData, bchan, echan, nif, err);
  outVer = BPSoln;
  if (outVer<=0) 
    outVer = ObitTableListGetHigh (outUV->tableList, "AIPS BP")+1;
  Obit_log_error(err, OBIT_InfoErr, 
		 "New BP table ver. %d", outVer);

  /* loop over blocks of data */
  iVis = 0; nVis = avgData->myDesc->nvis;
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    retCode = ObitUVReadSelect (avgData, NULL, err);
    if (retCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) goto cleanup;
    
    /* how much data? */
    nvis  = desc->numVisBuff;
    
    /* Loop over visibilities in buffer */
    vis = avgData->buffer+desc->nrparm;
    rec = avgData->buffer;
    /* Save first time */
    if (time<-1.0e10) {
      time = rec[desc->iloct];
      etime = time+solInt-1.156e-5; /* end time */
      souId = (olong)(rec[desc->ilocsu]+0.5);
    }
    for (ivis=0; ivis<nvis; ivis++) {
      iVis++; /* Count total */
      if (iVis>nVis) goto done;
      /* Only autoCorrelations */
      ObitUVDescGetAnts(desc, rec, &ia1, &ia2, &suba);
      if (ia1!=ia2) continue;
      /* Are we there yet? */
      time = rec[desc->iloct];
      /* Finished solution? */
      if ((time>etime) || (souId!=((olong)(rec[desc->ilocsu]+0.5)))) {
	if (haveSome)
	  DumpBP(BPTable, phase, souId, time, bchan, echan, nif, numAnt, refAnt, err);
	if (err->error) goto cleanup;
	etime = time+solInt-1.156e-5; /* new end time */
	souId = (olong)(rec[desc->ilocsu]+0.5);
	haveSome = FALSE;
      } /* end finished accumulation */
      /* loop over IFs */
      ifvis = vis;
      for (iif=0; iif<nif; iif++) {
	fvis = ifvis;
	/* Channel Loop */
	for (ichan=bchan; ichan<=echan; ichan++) {
	  indx = (ia1-1)*nchanIF + iif*nc + ichan-bchan;
	  /* loop over frequencies being averaged */
	  sumr = sumi = sumwt = 0.0;
	  loLim = -chWid; if (ichan<chWid) loLim = -ichan;
	  hiLim =  chWid; if (ichan>nfreq-chWid) hiLim = nfreq-ichan;
	  for (ifreq=loLim; ifreq<=hiLim; ifreq++) {
	    afvis = fvis + ifreq*desc->incf;
	    wt  = *(afvis+2);
	    if (wt>0.0) {
	      sumr += wt*(*(afvis)); sumi += wt*(*(afvis+1));
	      sumwt += wt;
	    }
	  } /* end loop over freq averaging */
	  if (sumwt>0.0) phase[indx] = atan2(sumi, sumr);
	  else           phase[indx] = fblank;
	  haveSome = TRUE;
	  fvis += desc->incf; /* visibility pointer */
	} /* end loop over frequencies */
	ifvis += desc->incif; /* visibility pointer */
      } /* Loop over IFs */
	/* update data pointers */
      vis += desc->lrec;
      rec += desc->lrec;
    } /* end loop over visibilities in buffer */
  } /* end loop over file */
    
  /* Close data */
 done:
  retCode = ObitUVClose (avgData, err);
  if (err->error) goto cleanup;
    
  /* Finished write remaining output */
  if (haveSome)
    DumpBP(BPTable, phase, souId, time, bchan, echan, nif, numAnt, refAnt, err);
  if (err->error) goto cleanup;
  
  /* Copy BP table to outUV */
  ver = BPSoln;
  BPOut = newObitTableBPValue ("Output BP", (ObitData*)outUV, &ver,
			       OBIT_IO_WriteOnly, nstok, nif, nchan, err);
  if (err->error) Obit_traceback_msg (err, routine,outUV->name);

  BPTable = ObitTableBPCopy (BPTable, BPOut, err);
  if (err->error) goto cleanup;
  
  /* Cleanup */
  cleanup:
  if (phase) g_free(phase);
  ObitUVZapTable (avgData, "AIPS BP", -1, err);
  BPTable = ObitTableBPUnref(BPTable );
  BPOut  = ObitTableBPUnref(BPOut);
  if (err->error) Obit_traceback_msg (err, routine, avgData->name);
} /* end XYBandpassCal  */

/** 
 * Write solution interval phases to BP table
 * \param BPTable BP table to write
 * \param phase   Array of phases in order, chan, IF, antenna
 * \param souId   Source Id.
 * \param time    Time in days to write entry
 * \param bchan   First channel calibrated
 * \param echan   Highest channel calibrated
 * \param nif     Number of IFs
 * \param numAnt  Number of antennas
 * \param refAnt  Reference antenna, if >0 then subtract from all other antennas
 * \param err     ObitErr stack for reporting problems.
 */
void DumpBP(ObitTableBP* BPTable, ofloat *phase, olong souId, ofloat time, 
	    olong bchan, olong echan, olong nif, olong numAnt, olong refAnt, 
	    ObitErr* err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  olong iant, orow, nchan, npol=2, ichan, iif;
  olong indx, jndx, rndx, nchanIF;
  ofloat fblank = ObitMagicF();
  ofloat *refPhase=NULL;
  gchar *routine = "DumpBP";

  /* error checks */
  if (err->error) return;

  nchan   = (echan-bchan+1);
  nchanIF = nif * nchan;  /* Number of channels * IFs */

  /* Save reference antenna phases */
  if (refAnt>0.0) {
    refPhase = g_malloc0(nchanIF*sizeof(ofloat));
    for (iif=0; iif<nif; iif++) {
      for (ichan=bchan; ichan<=echan; ichan++) {
	indx = (refAnt-1)*nchanIF + iif*nchan + ichan-bchan;
	rndx = iif*nchan + ichan-bchan;
	refPhase[rndx] =  phase[indx];
      }
    }
  }
  /* Open BP table */
  retCode = ObitTableBPOpen (BPTable, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTable->name);
  BPRow = newObitTableBPRow(BPTable);
  ObitTableBPSetRow (BPTable, BPRow, err);
  if (err->error) Obit_traceback_msg (err, routine, BPTable->name);

  /* Intitial values */
  BPTable->numAnt     = numAnt;
  BPRow->BW           = 0.0;
  BPRow->ChanShift[0] = 0.0;
  BPRow->ChanShift[1] = 0.0;
  BPRow->RefAnt1      = refAnt;
  BPRow->RefAnt2      = refAnt;
  BPRow->SubA         = 0;
  BPRow->FreqID       = 0;
  /* Reference? */
  if (refAnt>0) {
    for (iant=1; iant<=numAnt; iant++) {
      for (iif=0; iif<nif; iif++) {
	for (ichan=bchan; ichan<=echan; ichan++) {
	  rndx = iif*nchan + ichan-bchan;
	  indx = (iant-1)*nchanIF + iif*nchan + ichan-bchan;
	  if ((phase[indx]!=fblank) && (refPhase[rndx]!=fblank)) {
	    phase[indx] -= refPhase[rndx];
	  } else phase[indx] = fblank;
	}
      }
    }
  } /* end phase referencing */

  /* Loop over antennas  */
  for (iant=1; iant<=numAnt; iant++) {
    /* Set time, antenna etc.*/
    BPRow->Time   = time;
    BPRow->antNo  = iant;
    BPRow->SourID = souNo;

    /* loop over IFs */
    for (iif=0; iif<nif; iif++) {
      BPRow->Weight1[iif] = 1.0;
      if (npol>1) BPRow->Weight2[iif] = 1.0;
      
      /* Channel Loop */
      for (ichan=bchan; ichan<=echan; ichan++) {
	indx = (iant-1)*nchanIF + iif*nchan + ichan-bchan;
	jndx = iif*nchan + ichan -1;
	if (phase[indx]!=fblank) {
	  BPRow->Real1[jndx]   = 1.0;
	  BPRow->Imag1[jndx]   = 0.0;
	  if (npol>1) {   /* There better be or this is pretty pointless */
	    BPRow->Real2[jndx] = cos(-phase[indx]); /* As correction */
	    BPRow->Imag2[jndx] = sin(-phase[indx]);
	  }
 	} else { /* Flagged */
	  BPRow->Real1[jndx]   = fblank;
	  BPRow->Imag1[jndx]   = fblank;
	  if (npol>1) { 
	    BPRow->Real2[jndx] = fblank;
	    BPRow->Imag2[jndx] = fblank;
	  }
 	}
      } /* end channel loop */
    } /* end IF Loop */

    /* Write output table */
    orow = -1;
    retCode = ObitTableBPWriteRow (BPTable, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over antennas */

  /* Close BP table */
cleanup:
  retCode = ObitTableBPClose (BPTable, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTable->name);
  BPRow = ObitTableBPRowUnref(BPRow);
  if (refPhase) g_free(refPhase);
}  /* end DumpBP */

/** 
 * Create Bandpass table
 * \param inData  UV onto which the BP table is to be attached
 * \param bchan   First channel 
 * \param echan   Highest channel 
 * \param nif     Number of IFs
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* MakeBPTable (ObitUV* inData, olong bchan, olong echan, 
			  olong nif, ObitErr *err)
{
  ObitTableBP *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  olong nchan, npol=2, ver, numAnt;
  ObitInfoType type;
  ObitUVDesc *desc=inData->myDesc;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "MakeBPTable";

  /* error checks */
  if (err->error) return BPOut;
  g_assert (ObitUVIsA(inData));

  /* Table info - channels */
  desc = (ObitUVDesc*)inData->myIO->myDesc;
  if (desc->jlocf>=0)
    nchan = desc->inaxes[desc->jlocf];
  else
    nchan = 1;  /* Curious it */
  /* Stokes */
  if (desc->jlocs>=0)
    npol = MIN (2, desc->inaxes[desc->jlocs]);
  else
    npol = 1; 

  /* Create output BP table */
  ver = 0;
  ObitInfoListGetTest(myInput, "BPSoln",  &type, dim, &ver);
  BPOut = newObitTableBPValue ("Temp BP", (ObitData*)inData, &ver,
			       OBIT_IO_WriteOnly, npol, nif, nchan, err);
  if (err->error) goto cleanup;

  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)BPOut, err);

  /* Open BP table */
  retCode = ObitTableBPOpen (BPOut, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  BPRow = newObitTableBPRow(BPOut);
  ObitTableBPSetRow (BPOut, BPRow, err);
  if (err->error) goto cleanup;

  /* Set header values */
  numAnt  = inData->myDesc->numAnt[0];/* actually highest antenna number */
  BPOut->numAnt    = numAnt;  /* Max. antenna number */
  BPOut->numShifts = 0;
  BPOut->numChan   = nchan;
  BPOut->startChan = 1;
  BPOut->lowShift  = 1;
  BPOut->shiftInc  = 1;
  strncpy (BPOut->BPType, "          ", MAXKEYCHARTABLEBP);
  BPOut->myDesc->sort[0] = BPOut->TimeCol+1;  /* Sort order */
  BPOut->myDesc->sort[1] = BPOut->antNoCol+1;


 cleanup:

  /* Close BP table */
  retCode = ObitTableBPClose (BPOut, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = ObitTableBPRowUnref(BPRow);

  return BPOut;
} /* end MakeBPTable  */
