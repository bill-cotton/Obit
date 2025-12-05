/* $Id$  */
/* X-Y phase bandpass calibration                                     */
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
#include "ObitUVGSolve.h"
#include "ObitSkyModel.h"
#include "ObitTableSUUtil.h"
#include "ObitUVUtil.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableCL.h"
#include "ObitTableCLUtil.h"
#include "ObitTableBP.h"
#include "ObitUVSoln2Cal.h"
#include "ObitSinCos.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* XYPassIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void XYPassOut (ObitInfoList* outList, ObitErr *err);
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
/* Do channel solutions */
void XYBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		   ObitErr* err);
/* Write history */
void XYPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Initial calibration */
ObitUV* InitialCal (ObitInfoList* myInput, ObitUV* scrData, ObitErr* err);
/* Create  BP table */
ObitTableBP* MakeBPTable (ObitUV* inData, ofloat time,
			  olong bchan, olong echan, olong nif, 
			  ofloat *phase, ofloat *snr, olong refAnt, ObitErr *err);
/* Update  BP table */
ObitTableBP* UpdateBPTable (ObitUV* inData, olong BPver, ObitUV* outData, 
			    ofloat time, olong bchan, olong echan, olong nif, 
			    ofloat *phase, ofloat *snr, olong refAnt,  ObitErr *err);
/* Get source lookup table */
olong* GetLookupSU (ObitSourceList *SList, ObitInfoList *info);
/* Setup EVPA sin/cos arrays */
void SetupEVPA (olong numIF, olong numFreq, ObitUV *inUV, 
		ofloat EVPA, ofloat RM, ofloat *EVPAFS, ofloat *EVPAFC);
/* Calculate sign flips */
void GetFlip (olong numIF, olong numFreq, ObitUV *inUV, ofloat *EVPAFS, ofloat *EVPAFC, 
	      ofloat chisin, ofloat chicos, ofloat *flip);
/* Human readable time string */
void day2dhms(ofloat time, gchar *timeString);


/* Program globals */
gchar *pgmName = "XYPass";      /* Program name */
gchar *infile  = "XYPass.in" ;  /* File with program inputs */
gchar *outfile = "XYPass.out";  /* File to contain program outputs */
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
/*    X-Y phase bandpass calibration                                      */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *scrData = NULL, *avgData = NULL;;
  ObitErr      *err= NULL;
  ofloat       ftemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = XYPassIn (argc, argv, err);
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
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* Save first source id in the likely case it's the only one */
  if (inData->mySel->sources) souNo = inData->mySel->sources[0];
  
  /* Index scrData */
  dim[0] = dim[1] = 1;
  ftemp = 15.0;  /* Max scan time 15 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxScan", OBIT_float, dim, &ftemp);
  ftemp = 1.0; /* Max allowable gap 1 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxGap", OBIT_float, dim, &ftemp);
  ObitUVUtilIndex (scrData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* Initial calibration */
  avgData =  InitialCal(myInput, scrData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* Do channel solutions, convert to BP table */
  XYBandpassCal(myInput, avgData, inData, err);
  if (err->error) {ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;}

  /* Write history */
  XYPassHistory (myInput, inData, err); 
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  scrData   = ObitUnref(scrData);
  avgData   = ObitUnref(avgData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* XYPassIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "XYPassIn";

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
} /* end XYPassIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: XYPass -input file -output ofile [args]\n");
    fprintf(stderr, "XYPass Obit task determine bandpass for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def XYPass.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def XYPass.out\n");
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
  ofloat farray[3];
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
  strTemp = "XYPass.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "XYPassName";
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
  //ObitInfoType type;
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
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "FreqID", "BChan", "EChan", "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "doPol", "PDVer", "keepLin", "Mode", "ModelType", "Alpha", "refAnt",
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
/*  Write History for XYPass                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void XYPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA",
    "selBand", "selFreq", "FreqID", 
    "BChan1", "EChan1", "BChan2", "EChan2", "ChWid2", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol", "PDVer", "keepLin", "Antennas",  "refAnt",  
    "EVPA", "RM",
    "solInt1", "solInt2",  "minSNR",  "minNo", "BPSoln", "prtLv", "nThreads",
   NULL};
  gchar *routine = "XYPassHistory";

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
 
} /* end XYPassHistory  */

/**
 * Initial atmospheric calibration, calibrates data and averages to
 * solint2 and returns averaged data.
 * \param myInput Input parameters on InfoList    
 * \param scrUV   Data to be used to determine calibration
 * \param err     ObitErr stack for reporting problems.
 * \return time averaged visibility data
 */
ObitUV* InitialCal (ObitInfoList* myInput, ObitUV* scrData, ObitErr* err)
{
  ObitTableSN  *SNTable = NULL;
  ObitTableCL  *CLTable1 = NULL, *CLTable2 = NULL;
  ObitUVGSolve *solver=NULL;
  ObitUV       *avgData=NULL;
  gboolean     btemp, avgPol=TRUE;
  ofloat       ftemp, solInt;
  olong        itemp, refAnt, ver=1;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *calParms[] = {  /* Parameters to update CL with SN tables */
    "subA", "Antennas",
     NULL};
  gchar        *solverParms[] = {  /* Calibration parameters */
    "solnVer", "solType", "solMode", "avgPol", "avgIF", "doMGM", "elevMGM",
    "refAnt", "ampScalar", "minSNR",  "minNo", "prtLv",
    NULL};
  gchar *blank = "    ";
  gchar *FQInclude[] = {"AIPS FQ", "AIPS AN", NULL};
  gchar *routine = "InitialCal";

  /* error checks */
  if (err->error) return avgData;
  g_assert (ObitUVIsA(scrData));
  
  /* Initial CL table on scratch file */
  ObitInfoListGet(myInput, "solInt1", &type, dim, &ftemp, err);
  ObitInfoListAlwaysPut(scrData->info, "solInt", OBIT_float, dim, &ftemp);
  CLTable1 =  ObitTableCLGetDummy (scrData, scrData, ver, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);
  ObitInfoListGet(myInput, "refAnt", &type, dim, &refAnt, err);
 
  /* Create solver */
  solver = ObitUVGSolveCreate("Gain solver");

  /* Copy calibration control to solver */
  ObitInfoListCopyList (myInput, solver->info, solverParms);
  ObitInfoListGet(myInput, "solInt1", &type, dim, &ftemp, err);
  ObitInfoListAlwaysPut(solver->info, "solInt", OBIT_float, dim, &ftemp);
  itemp = 1;
  ObitInfoListAlwaysPut(solver->info, "solnVer", OBIT_long, dim, &itemp);
  /* Average polarizations */
  ObitInfoListAlwaysPut(solver->info, "avgPol", OBIT_bool, dim, &avgPol);
  ObitInfoListGet(myInput, "BChan1", &type, dim, &itemp, err);
  ObitInfoListAlwaysPut(scrData->info, "BChan", OBIT_long, dim, &itemp);
  ObitInfoListGet(myInput, "EChan1", &type, dim, &itemp, err);
  ObitInfoListAlwaysPut(scrData->info, "EChan", OBIT_long, dim, &itemp);
  dim[0] = strlen(blank);  /* Phase only */
  ObitInfoListAlwaysPut(scrData->info, "solMode", OBIT_string, dim, blank);

  /* Do atmospheric solution */
  SNTable = ObitUVGSolveCal (solver, scrData, scrData, scrData->mySel, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Apply calibration to CL table */
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "solnVer", OBIT_long, dim, &SNTable->tabVer);
  ObitInfoListAlwaysPut(scrData->info, "calIn",   OBIT_long, dim, &CLTable1->tabVer);
  ObitInfoListAlwaysPut(scrData->info, "refAnt",  OBIT_long, dim, &refAnt);
  itemp = CLTable1->tabVer+1;
  ObitInfoListAlwaysPut(scrData->info, "calOut",  OBIT_long, dim, &itemp);
  ObitInfoListCopyList (myInput, scrData->info, calParms);
  CLTable2 = ObitUVSoln2Cal (scrData, scrData, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Apply calibration */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = 2;
  ObitInfoListAlwaysPut(scrData->info, "doCalib", OBIT_long, dim, &itemp);
  itemp = CLTable1->tabVer+1;
  ObitInfoListAlwaysPut(scrData->info, "gainUse", OBIT_long, dim, &itemp);
  itemp = 1;   /* Copy all channels */
  ObitInfoListAlwaysPut(scrData->info, "BChan", OBIT_long, dim, &itemp);
  itemp = 0;
  ObitInfoListAlwaysPut(scrData->info, "EChan", OBIT_long, dim, &itemp);

  /* Average data to solInt2; 0=> all */
  solInt = 0.0;
  ObitInfoListGetTest (myInput, "solInt2", &type, dim, &solInt); 
  if (solInt<=1.0e-5) solInt = 1000.0; 
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "timeAvg", OBIT_float, dim, &solInt);
  avgData = ObitUVUtilAvgT(scrData, TRUE, avgData, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Be sure FQ table copied */
   ObitUVCopyTables (scrData, avgData, NULL, FQInclude, err);

  /* cleanup */
  solver    = ObitUVGSolveUnref(solver);
  SNTable   = ObitTableSNUnref(SNTable);
  CLTable1  = ObitTableSNUnref(CLTable1);
  CLTable2  = ObitTableSNUnref(CLTable2);

  return avgData;
} /* end InitialCal */

/**
 * Determine Bandpass table
 * \param myInput Input parameters on InfoList    
 * \param avgUV   Averaged data to be used to determine calibration
 * \param inUV    UV data onto which the BP table to be written
 * \param err     ObitErr stack for reporting problems.
 */
void  XYBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		    ObitErr* err)
{
  ObitTableBP *BPTable=NULL, *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableAN *ANTab=NULL;
  ObitAntennaList *ANList=NULL;
  ObitSource *Source=NULL;
  ObitSourceList *SList=NULL;
  olong nchan, bchan2=0, echan2=0, chinc2=0, nc, i, j, suID, lastSU=-10;
  olong ivis, nvis, ifreq, nfreq, iif, nif, indx, count=0;
  olong itemp, nchanIF, doBand, BPVer, outVer, npol, BPSoln, refAnt;
  olong numOrb, numPCal, numIF, iANver, lo, hi, kndx;
  olong prtLv=0, fitType, *lookupSU = NULL;
  gboolean btemp, newTime=FALSE, multiSU=FALSE;
  ofloat fblank = ObitMagicF();
  ofloat *EVPA=NULL, *EVPAFS=NULL, *EVPAFC=NULL, *RM=NULL, *flip=NULL;
  ofloat time=-1.0e20, lastTime=0.0, chi, chisin, chicos;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *vis, *ifvis, *fvis, *svis, *rec;
  ObitUVDesc *desc;
  ofloat *XYrsum=NULL, *YXrsum=NULL, *XYwt=NULL, *YXwt=NULL;
  ofloat *XYisum=NULL, *YXisum=NULL, *phase=NULL, *snr = NULL;
  ofloat *XYrsum2=NULL, *YXrsum2=NULL, *XYisum2=NULL, *YXisum2=NULL;
  ofloat *XYrzum=NULL, *YXrzum=NULL, *XYzwt=NULL, *YXzwt=NULL;
  ofloat *XYizum=NULL, *YXizum=NULL;
  ofloat xyr, xyi, yxr, yxi, xywt, yxwt, Sumr=0., Sumr2=0., Sumi=0., Sumi2=0., Sumw=0., minSNR=0.0;
  ofloat XYSumr=0., XYSumi=0., YXSumr=0., YXSumi=0., XYWt=0., YXWt=0.;
  ofloat Var_r, Var_i, Var_ph, sumTime=0.0;
  olong ant1, ant2, subA, cntTime=0;
  gchar *Stokes="    ", timeStr[24];
  gchar *routine = "XYBandpassCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(avgData));

  /* Any data in avgData? */
  if (avgData->myDesc->nvis<=0) {
    Obit_log_error(err, OBIT_Error, "NO averaged/calibrated data");
    return;
  }

  /* Make sure 4 Stokes correlations are available */
  Obit_return_if_fail(((avgData->myDesc->inaxes[avgData->myDesc->jlocs]>=4)), err,
		      "%s: MUST have cross polarized data", 
		      routine);
  
 /* Better be linear feeds? */
  Obit_return_if_fail((avgData->myDesc->crval[avgData->myDesc->jlocs]<-4.0),
		      err, "%s: MUST be linear feed data", routine);

 /* Antenna List for linear feeds */
  numOrb  = 0;  numPCal = 0; numIF = 0; iANver=1;
  ANTab = 
    newObitTableANValue (inData->name, (ObitData*)inData, &iANver, OBIT_IO_ReadOnly, 
			 numIF, numOrb, numPCal, err);
  ANList = ObitTableANGetList (ANTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Number of channels, IFs */
  desc  = avgData->myDesc;
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
 
  /* Range of channels */
  nchan = desc->inaxes[desc->jlocf];
  ObitInfoListGetTest (myInput, "BChan2", &type, dim, &bchan2); 
  if (bchan2<=0) bchan2 = 1;
  ObitInfoListGetTest (myInput, "EChan2", &type, dim, &echan2); 
  if (echan2<=0) echan2 = nchan;
  echan2 = MIN (echan2, nchan);
  /* How many channels to average */
  ObitInfoListGetTest (myInput, "ChWid2", &type, dim, &chinc2); 
  chinc2 /= 2;  /* How many above and below? */
  chinc2 = MAX (0, chinc2);
  /* Minimum SNR */
  ObitInfoListGetTest (myInput, "minSNR", &type, dim, &minSNR); 
  /* prtLv */
  ObitInfoListGetTest (myInput, "prtLv", &type, dim, &prtLv); 
 /* Stokes */
  if (avgData->myDesc->jlocs>=0)
    npol = MIN (2, avgData->myDesc->inaxes[avgData->myDesc->jlocs]);
  else
    npol = 1;

  /* EVPA, RM */
  ObitInfoListGetP(myInput, "EVPA", &type, dim, (gpointer)&EVPA);
  ObitInfoListGetP(myInput, "RM", &type, dim, (gpointer)&RM);
  fitType = 0;
  ObitInfoListGetTest(myInput, "fitType", &type, dim, &fitType);
 
  /* Create work arrays */
  EVPAFS =  g_malloc0((nfreq*(nif+2))*sizeof(ofloat));
  EVPAFC =  g_malloc0((nfreq*(nif+2))*sizeof(ofloat));
  flip   =  g_malloc0((nfreq*(nif+2))*sizeof(ofloat));
  snr    =  g_malloc0((nfreq*(nif+2))*sizeof(ofloat));
  phase  =  g_malloc0((nfreq*(nif+2))*sizeof(ofloat));
  for (j=0; j<nfreq*nif; j++) flip[j] = 1.0;  /* Default no flip */

  /* MultiSource? */
  multiSU = (avgData->myDesc->ilocsu>=0); 
  if (!multiSU) {
    suID = -1;
    Source= ObitUVGetSource (Source, inData, suID, err);
    /* Source lookup table */
    lookupSU   = g_malloc0(2*sizeof(olong));
    lookupSU[0]= 0;
    /* Setup EVPA arrays */
    SetupEVPA(nif, nfreq, inData, EVPA[0], RM[0], EVPAFS, EVPAFC);
  } else { /* end single source */
    /* Multi source, need sourceList, lookup table */
    SList = ObitUVGetSourceList (inData, err); /* Get full Source List */
    lookupSU = GetLookupSU (SList, inData->info);
  } /* end multiSource */
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Target */
  ObitInfoListGetTest (myInput, "RM",      &type, dim, RM); 
  ObitInfoListGetTest (myInput, "EVPA", &type, dim, EVPA); 
  for (i=0; i<dim[0]; i++) EVPA[i] *= DG2RAD;    /* to radians */

  /* Prior bandpass calibration */
  ObitInfoListGetTest (myInput, "doBand", &type, dim, &doBand); 
  ObitInfoListGetTest (myInput, "BPVer",  &type, dim, &BPVer); 
  ObitInfoListGetTest (myInput, "BPSoln", &type, dim, &BPSoln);
  if (BPVer<=0) BPVer = ObitTableListGetHigh (inData->tableList, "AIPS BP");
  ObitInfoListGetTest (myInput, "refAnt", &type, dim, &refAnt); 

  /* Accumulator arrays */
  nc      =  (echan2-bchan2+1);
  nchanIF =  nif * nc;
  XYrsum  = g_malloc0(nchanIF*sizeof(ofloat));
  XYwt    = g_malloc0(nchanIF*sizeof(ofloat));
  YXrsum  = g_malloc0(nchanIF*sizeof(ofloat));
  YXwt    = g_malloc0(nchanIF*sizeof(ofloat));
  XYisum  = g_malloc0(nchanIF*sizeof(ofloat));
  YXisum  = g_malloc0(nchanIF*sizeof(ofloat));
  XYrsum2 = g_malloc0(nchanIF*sizeof(ofloat));
  YXrsum2 = g_malloc0(nchanIF*sizeof(ofloat));
  XYisum2 = g_malloc0(nchanIF*sizeof(ofloat));
  YXisum2 = g_malloc0(nchanIF*sizeof(ofloat));
  /* Another set for averaging over baselines */
  XYrzum  = g_malloc0(nchanIF*sizeof(ofloat));
  XYzwt   = g_malloc0(nchanIF*sizeof(ofloat));
  YXrzum  = g_malloc0(nchanIF*sizeof(ofloat));
  YXzwt   = g_malloc0(nchanIF*sizeof(ofloat));
  XYizum  = g_malloc0(nchanIF*sizeof(ofloat));
  YXizum  = g_malloc0(nchanIF*sizeof(ofloat));

  /* Don't need to apply calibration */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(avgData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = -1;
  ObitInfoListAlwaysPut(avgData->info, "doCalib", OBIT_long, dim, &itemp);
  /* Want to convert to IQUV */
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut(avgData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Open data */
  retCode = ObitUVOpen (avgData, OBIT_IO_ReadCal, err);
 if (err->error) goto cleanup;
 
 /* loop over file */
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
   for (ivis=0; ivis<nvis; ivis++) {
     time = rec[inData->myDesc->iloct];
     ObitUVDescGetAnts(avgData->myDesc, rec, &ant1, &ant2, &subA);  /* get antennas */
     /* New source? */
     if (multiSU && (lastSU!=rec[inData->myDesc->ilocsu])){
       suID = rec[inData->myDesc->ilocsu]; lastSU = suID;
       Source = ObitUVGetSource (Source, inData, suID, err);
       if (err->error) goto cleanup;
       /* Setup EVPA arrays */
       SetupEVPA(nif, nfreq, inData, EVPA[lookupSU[suID]], RM[lookupSU[suID]], 
		 EVPAFS, EVPAFC); 
     } /* end new source */
     
     /* New time? */
     if (time > lastTime) {
       cntTime++; sumTime+=time;  /* For average time */
       /* linear feeds - flip phase if linP (eff. Stokes u) fn positive */
       chi  = ObitAntennaListParAng(ANList, 1,  time, Source);  /* Parallactic angle */
       chisin = sinf(chi*2.0); chicos = cosf(chi*2.0); 
       lastTime = time;
       /* Calculate phase flips */
       GetFlip(nif, nfreq, inData, EVPAFS, EVPAFC, chisin, chicos, flip);
       newTime = TRUE;
     }  /* end new time */

     /* if newTime accumulate average over baseine */
     if (newTime) {
       /* Accumulate average over baseline */
       /* loop over IFs */
       for (iif=0; iif<nif; iif++) {
	 /* loop over frequencies */
	 for (ifreq = 0; ifreq<nfreq; ifreq++) {
	   indx = iif*nc + ifreq-bchan2+1;
	   kndx = iif*nif + ifreq;
	   /* Normalize */
	   if (XYzwt[indx]>0.0) {XYSumr = XYrzum[indx]/XYzwt[indx]; 
	     XYSumi = XYizum[indx]/XYzwt[indx];
	     XYWt = XYzwt[indx];}
	   else {XYSumr = XYSumi = 0.0;}
	   if (YXzwt[indx]>0.0) {YXSumr = YXrzum[indx]/YXzwt[indx]; 
	     YXSumi = YXizum[indx]/YXzwt[indx];
	     YXWt = YXzwt[indx];}
	   else {YXSumr = YXSumi = 0.0;}
	   // DEBUG
	   //if ((prtLv>=2) && (indx==91) && (XYWt>0.0)) {
	     //day2dhms(time, timeStr); /* Human readable time */
	     //XYSumr = XYrzum[indx]/XYzwt[indx]; XYSumi = XYizum[indx]/XYzwt[indx]; 
	     //YXSumr = YXrzum[indx]/YXzwt[indx]; YXSumi = YXizum[indx]/YXzwt[indx]; 
	     //fprintf(stderr,"%s ch %d %d %d phases %f %f  %f %f  %f %f  %f %f \n", 
	//	     timeStr,indx+1,iif,ifreq,
	//	     57.296*atan2f(XYSumi,XYSumr), 57.296*atan2f(YXSumi,YXSumr),
	//	     XYSumr,XYSumi,YXSumr,YXSumi, XYzwt[indx], XYzwt[indx]);
	   //} // end DEBUG
	   /* Reset accumulators */
	   XYrzum[indx] = XYizum[indx] = XYzwt[indx] = 0.0;
	   YXrzum[indx] = YXizum[indx] = YXzwt[indx] = 0.0;
	   /* Accumulate baseline averages, adjust both XY and YX to the same */
	   if (flip[kndx]>0) {
	     /* XY */
	     if (XYWt>0.0) {
	       XYrsum[indx] += XYSumr * XYWt; XYrsum2[indx] += XYSumr * XYSumr * XYWt;
	       XYisum[indx] += XYSumi * XYWt; XYisum2[indx] += XYSumi * XYSumi * XYWt;
	       XYwt[indx]   += XYWt;
	     }
	     /* YX */
	     if (YXWt>0.0) {
	       YXrsum[indx] += YXSumr * YXWt; YXrsum2[indx] += YXSumr * YXSumr * YXWt;
	       YXisum[indx] -= YXSumi * YXWt; YXisum2[indx] += YXSumi * YXSumi * YXWt;
	       YXwt[indx]   += YXWt;
	     }
	   } else { /* negate the sign of the other */
	     /* XY */
	     if (XYWt>0.0) {
	       XYrsum[indx] += XYSumr * XYWt; XYrsum2[indx] += XYSumr * XYSumr * XYWt;
	       XYisum[indx] -= XYSumi * XYWt; XYisum2[indx] += XYSumi * XYSumi * XYWt;
	       XYwt[indx]   += XYWt;
	     }
	     /* YX */
	     if (YXWt>0.0) {
	       YXrsum[indx] += YXSumr * YXWt; YXrsum2[indx] += YXSumr * YXSumr * YXWt;
	       YXisum[indx] += YXSumi * YXWt; YXisum2[indx] += YXSumi * YXSumi * YXWt;
	       YXwt[indx]   += YXWt;
	     }
	   } /* end flip */
	 } /* end frequency loop */
       } /* end IF Loop */
       newTime = FALSE; /* No longer TRUE */
     } /* end newTime accumupation */
     
     /* Average over baseline */
     /* loop over IFs */
     ifvis = vis;
     for (iif=0; iif<nif; iif++) {
       /* loop over frequencies */
       fvis = ifvis;
       for (ifreq = 0; ifreq<nfreq; ifreq++) {
	 indx = iif*nc + ifreq-bchan2+1;
	 /* Loop over blocks of channels being averaged */
	 lo = MAX(0, ifreq-chinc2); hi = MIN(nfreq, ifreq+chinc2);
	 for (j=lo; j<=hi;j++) {
	   /*  Stokes correlations */
	   svis = ifvis + j*desc->incf;
	   /* XX,YY - not needed */
	   /* extract data */
	   svis += 2*desc->incs; /* visibility pointer */
	   xyr = *svis; xyi = *(svis+1); xywt = *(svis+2);
	   svis += desc->incs; 
	   yxr =  *svis; yxi = *(svis+1); yxwt = *(svis+2);
	   // DEBUG  IF 1, ch 92, baseline 40-50
	   // if ((prtLv>=2) && (ant1==40) && (ant2==50)&& (indx==91)) {
	   //   day2dhms(time, timeStr); /* Human readable time */
	   //   XYSumr = xyr; XYSumi = xyi; 
	   //   YXSumr = yxr; YXSumi = yxi; 
	   //   fprintf(stderr,"%s Vis ch %d %d %d phases %f %f  %f %f  %f %f  \n", 
	// 	   timeStr,indx+1,iif,ifreq,
	// 	   57.296*atan2f(XYSumi,XYSumr), 57.296*atan2f(YXSumi,YXSumr),
	// 	   XYSumr,XYSumi,YXSumr,YXSumi);
	 // } // end DEBUG
	   /* XY */
	   if (xywt>0.0) {
	     XYrzum[indx] += xyr * xywt;  XYizum[indx] += xyi * xywt; 
	     XYzwt[indx]   += xywt;
	   }
	   /* YX */
	   if (yxwt>0.0) {
	     YXrzum[indx] += yxr * yxwt;  YXizum[indx] += yxi * yxwt; 
	     YXzwt[indx]   += yxwt;
	   }
	 } /* end inner loop */
	 fvis += desc->incf; /* visibility pointer */
       } /* end loop over frequencies */
       ifvis += desc->incif; /* visibility pointer */
     } /* Loop over IFs */    
     /* update data pointers */
     vis += desc->lrec; rec += desc->lrec;
   } /* end loop over buffer */
 } /* end loop over data */
 /* Accumulate final average over baseline if any */
 /* loop over IFs */
 for (iif=0; iif<nif; iif++) {
   /* loop over frequencies */
   for (ifreq = 0; ifreq<nfreq; ifreq++) {
     indx = iif*nc + ifreq-bchan2+1;
     kndx = iif*nif + ifreq;
     /* Normalize */
     if (XYzwt[indx]>0.0) {XYSumr = XYrzum[indx]/XYzwt[indx]; 
       XYSumi = XYizum[indx]/XYzwt[indx];
       XYWt = XYzwt[indx];}
     else {XYSumr = XYSumi = 0.0;}
     if (YXzwt[indx]>0.0) {YXSumr = YXrzum[indx]/YXzwt[indx]; 
       YXSumi = YXizum[indx]/YXzwt[indx];
       YXWt = YXzwt[indx];}
     else {YXSumr = YXSumi = 0.0;}
     // DEBUG
     //if ((prtLv>=2) && (indx==91) && (XYWt>0.0)) {
       //day2dhms(time, timeStr); /* Human readable time */
       //XYSumr = XYrzum[indx]/XYzwt[indx]; XYSumi = XYizum[indx]/XYzwt[indx]; 
       //YXSumr = YXrzum[indx]/YXzwt[indx]; YXSumi = YXizum[indx]/YXzwt[indx]; 
       //fprintf(stderr,"%s ch %d %d %d phases %f %f  %f %f  %f %f  %f %f \n", 
	//       timeStr,indx+1,iif,ifreq,
	    //   57.296*atan2f(XYSumi,XYSumr), 57.296*atan2f(YXSumi,YXSumr),
	    //   XYSumr,XYSumi,YXSumr,YXSumi, XYzwt[indx], XYzwt[indx]);
     //} // end DEBUG
     /* Accumulate baseline averages, adjust both XY and YX to the same */
     if (flip[kndx]>0) {
       /* XY */
       if (XYWt>0.0) {
	 XYrsum[indx] += XYSumr * XYWt; XYrsum2[indx] += XYSumr * XYSumr * XYWt;
	 XYisum[indx] += XYSumi * XYWt; XYisum2[indx] += XYSumi * XYSumi * XYWt;
	 XYwt[indx]   += XYWt;
       }
       /* YX */
       if (YXWt>0.0) {
	 YXrsum[indx] += YXSumr * YXWt; YXrsum2[indx] += YXSumr * YXSumr * YXWt;
	 YXisum[indx] -= YXSumi * YXWt; YXisum2[indx] += YXSumi * YXSumi * YXWt;
	 YXwt[indx]   += YXWt;
       }
     } else { /* negate the sign of the other */
       /* XY */
       if (XYWt>0.0) {
	 XYrsum[indx] += XYSumr * XYWt; XYrsum2[indx] += XYSumr * XYSumr * XYWt;
	 XYisum[indx] -= XYSumi * XYWt; XYisum2[indx] += XYSumi * XYSumi * XYWt;
	 XYwt[indx]   += XYWt;
       }
       /* YX */
       if (YXWt>0.0) {
	 YXrsum[indx] += YXSumr * YXWt; YXrsum2[indx] += YXSumr * YXSumr * YXWt;
	 YXisum[indx] += YXSumi * YXWt; YXisum2[indx] += YXSumi * YXSumi * YXWt;
	 YXwt[indx]   += YXWt;
       }
     } /* end flip */
   } /* end frequency loop */
 } /* end IF Loop */

 /* Close data */
 retCode = ObitUVClose (avgData, err);
 if (err->error) goto cleanup;

 /* Average time */
 sumTime /= cntTime;

 /* Convert to angle, SNR */
 for (i=0; i<nchanIF; i++) {
   // Check signs
   if (fitType==0)      {Sumr=XYrsum[i]+YXrsum[i];    Sumi=XYisum[i]+YXisum[i];
                         Sumr2=XYrsum2[i]+XYrsum2[i]; Sumi2=XYisum2[i]+YXisum2[i];
                         Sumw = XYwt[i] + YXwt[i];}
   else if (fitType==1) {Sumr=XYrsum[i];     Sumi=XYisum[i];
                         Sumr2=XYrsum2[i];   Sumi2=XYisum2[i];;
                         Sumw = XYwt[i];}
   else if (fitType==2) {Sumr=-YXrsum[i];    Sumi=-YXisum[i];
                         Sumr2=YXrsum2[i];   Sumi2=YXisum2[i];;
                         Sumw = YXwt[i];}
   /* Normalize, get phase, SNR */
   if (Sumw>0.0) {
     Sumr /= Sumw; Sumi /= Sumw; Sumr2 /= Sumw; Sumi2 /= Sumw;
     if (Sumr<0.0) {Sumr = -Sumr; Sumi = -Sumi;}  /* want solution closer to 0 */
     phase[i] = atan2f(Sumi,Sumr);
     /* Statistics: Signal to noise assuming a phase RMS of 1 radian is SNR = 1
	Var_ri = (Sum2-Sum1*Sum1)
	Var_phase=(1/(Sumr+Sumi*Sumi/Sumr))*Var_i + (Sumi/(1+Sumi*Sumi))*Var_r
	i.e. snr = 1/sqrt(Var_phase)      */
     Var_r = Sumr2 - Sumr*Sumr; Var_i = Sumi2 - Sumi*Sumi;
     Var_ph =(1/(Sumr+Sumi*Sumi/Sumr))*Var_i - (Sumi/(1+Sumi*Sumi))*Var_r;
     if (Var_ph>0.0) snr[i] = 1./ sqrtf(Var_ph);
     else            snr[i] = 0.0;
   } else {phase[i] = fblank; snr[i] = 0.0;}
   /* Check minSNR */
   if (snr[i]<minSNR) {phase[i] = fblank; snr[i] = 0.0;}
   if (snr[i]>0.0) count++;
 } /* End channel loop */
 
 /* tell about it */
 if (count>0) {
   Obit_log_error(err, OBIT_InfoErr, "Calibrated %d channels", count);
 } else {
   Obit_log_error(err, OBIT_Error, "NO channels calibrated");
   goto cleanup;
 }

  /* Diagnostic listing */
  if (prtLv>=2) {
    day2dhms(sumTime, timeStr); /* Human readable time */
    Obit_log_error(err, OBIT_InfoErr, 
		   "Average time %s phase, SNR ", timeStr);
    for (i=0; i<nchanIF; i+=4) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "%d %6.3f %5.0f, %6.3f %5.0f, %6.3f %5.0f, %6.3f %5.0f",
		     i,57.296*phase[i],snr[i],57.296* phase[i+1],snr[i+1],
		     57.296*phase[i+2],snr[i+2], 57.296*phase[i+3],snr[i+3]);
   } /* end diagnostic listing */
  }
  
  /* Create new or update prior depending on doBand */
  if (doBand<=0) {
   /* Create/populate BP table */
   BPTable = MakeBPTable (avgData, sumTime, bchan2, echan2, nif, 
			  phase,  snr, refAnt, err);
   outVer = BPSoln;
   if (outVer<=0) 
     outVer = ObitTableListGetHigh (inData->tableList, "AIPS BP")+1;
   Obit_log_error(err, OBIT_InfoErr, 
		  "New BP table ver. %d", outVer);
  } else {
   BPTable = UpdateBPTable (inData, BPVer, avgData, sumTime, bchan2, echan2, nif, 
			    phase,  snr, refAnt, err);
    outVer = BPSoln;
    if (outVer<=0) 
      outVer = ObitTableListGetHigh (inData->tableList, "AIPS BP")+1;
    Obit_log_error(err, OBIT_InfoErr, 
		   "Update BP table ver. %d to ver. %d", BPVer, outVer);
  }

  /* Create output BP table */
  BPOut = newObitTableBPValue ("Temp BP", (ObitData*)inData, &outVer,
			       OBIT_IO_WriteOnly, npol, nif, nchan, err);
  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)BPOut, err);
  if (err->error) goto cleanup;

  /* Copy BP table to inData */
  BPTable = ObitTableBPCopy (BPTable, BPOut, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  ObitErrLog(err); 
  if (XYrsum)    {g_free(XYrsum);}  if (YXrsum)    g_free(YXrsum);
  if (XYrsum2)   {g_free(XYrsum2);} if (YXrsum2)   g_free(YXrsum2);
  if (XYisum)    {g_free(XYisum);}  if (YXisum)    g_free(YXisum);
  if (XYisum2)   {g_free(XYisum2);} if (YXisum2)   g_free(YXisum2);
  if (XYrzum)    {g_free(XYrzum);}  if (YXrzum)    g_free(YXrzum);
  if (XYizum)    {g_free(XYizum);}  if (YXizum)    g_free(YXizum);
  if (XYwt)      {g_free(XYwt);}    if (YXwt)      g_free(YXwt);
  if (XYzwt)     {g_free(XYzwt);}   if (YXzwt)     g_free(YXzwt);
  if (flip)      {g_free(flip);}
  if (phase)     {g_free(phase);}   if (snr)       g_free(snr);
  if (lookupSU)  {g_free(lookupSU);}
  if (EVPAFS)    {g_free(EVPAFS);}  if (EVPAFC)    g_free(EVPAFC);
  ObitUVZapTable (avgData, "AIPS BP", -1, err);
  BPTable = ObitTableBPUnref(BPTable );
  ObitUVZapTable (avgData, "AIPS SN", -1, err);
  if (err->error) Obit_traceback_msg (err, routine, avgData->name);
  } /* end XYBandpassCal  */

/** 
 * Create Bandpass table, creates entries at time and for antennas numbers
 * up to the max in the data.
 * \param inData  UV onto which the BP table is to be attached
 * \param time    Time in days to write entry
 * \param bchan   First channel calibrated
 * \param echan   Highest channel calibrated
 * \param nif     Number of IFs
 * \param phase   Channel phase in radians
 * \param snr     Channel snr
 * \param refAnt  Reference antenna
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* MakeBPTable (ObitUV* inData, ofloat time,
			  olong bchan, olong echan, olong nif, 
			  ofloat *phase,  ofloat *snr, olong refAnt, ObitErr *err)
{
  ObitTableBP *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  ObitUVDesc *desc;
  olong i, irow, orow, nchan, npol, ver, numAnt, ichan, iif, nc;
  olong indx, jndx, suba=1;
  ofloat fblank = ObitMagicF();
  ObitInfoType type;
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
  if (err->error) Obit_traceback_val (err, routine,inData->name, BPOut);

  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)BPOut, err);

  /* Open BP table */
  retCode = ObitTableBPOpen (BPOut, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = newObitTableBPRow(BPOut);
  ObitTableBPSetRow (BPOut, BPRow, err);
  if (err->error) Obit_traceback_val (err, routine,inData->name, BPOut);

  /* Set header values */
  numAnt  = inData->myDesc->numAnt[suba-1];/* actually highest antenna number */
  BPOut->numAnt    = numAnt;               /* Max. antenna number */
  BPOut->numShifts = 0;
  BPOut->numChan   = nchan;
  BPOut->startChan = 1;
  BPOut->lowShift  = 1;
  BPOut->shiftInc  = 1;
  strncpy (BPOut->BPType, "          ", MAXKEYCHARTABLEBP);
  BPOut->myDesc->sort[0] = BPOut->TimeCol+1;  /* Sort order */
  BPOut->myDesc->sort[1] = BPOut->antNoCol+1;

  /* Initialize BP Row */
  desc = (ObitUVDesc*)inData->myIO->myDesc;
  BPRow->BW           = desc->cdelt[desc->jlocf];
  BPRow->ChanShift[0] = 0.0;
  BPRow->ChanShift[1] = 0.0;
  BPRow->RefAnt1      = refAnt;
  BPRow->RefAnt2      = refAnt;
  BPRow->SubA         = 0;
  BPRow->FreqID       = 0;
  BPRow->TimeI        = 24.0;
  BPRow->SourID       = souNo;
  for (i=0; i<nif; i++) BPRow->ChanShift[i] = 0.0;
  for (i=0; i<nif; i++) BPRow->Weight1[i]   = 1.0;
  if (npol>1) for (i=0; i<nif; i++) BPRow->Weight2[i] = 1.0;
  for (i=0; i<nchan*nif; i++) { 
    BPRow->Real1[i]   = fblank;
    BPRow->Imag1[i]   = fblank;
    if (npol>1) {
      BPRow->Real2[i]   = fblank;
      BPRow->Imag2[i]   = fblank;
    }
  }

  /* Fill actual values */
  nc      =  (echan-bchan+1);
  for (iif=0; iif<nif; iif++) {
    for (ichan=bchan; ichan<=echan; ichan++) {
      indx = iif*nc + ichan-bchan;
      jndx = iif*nchan + ichan -1;
      if (phase[indx]!=fblank) {
	BPRow->Real1[jndx]   = 1.0;
	BPRow->Imag1[jndx]   = 0.0;
	if (npol>1) {   /* There better be or this is pretty pointless */
	  BPRow->Real2[jndx]   = cos(phase[indx]);
	  BPRow->Imag2[jndx]   = sin(phase[indx]);
	  //SHIT BPRow->Weight2[jndx] = snr[indx];
	}
      } else {
	BPRow->Real1[jndx]   = fblank;
	BPRow->Imag1[jndx]   = fblank;
	if (npol>1) { 
	  BPRow->Real2[jndx]   = fblank;
	  BPRow->Imag2[jndx]   = fblank;
	   //SHIT BPRow->Weight2[jndx] = 0.0;
	}
      }
    } /* end loop over channel */
  } /* end loop over IF */
  
  /* Loop over antennas - entries all the same */
  for (irow=1; irow<=numAnt; irow++) {

    /* Set time, antenna etc.*/
    BPRow->Time   = time;
    BPRow->antNo  = irow;

    /* Write output table */
    orow = -1;
    retCode = ObitTableBPWriteRow (BPOut, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over antennas */

 cleanup:
  /* Close BP table */
  retCode = ObitTableBPClose (BPOut, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = ObitTableBPRowUnref(BPRow);

  return BPOut;
} /* end MakeBPTable  */
/** 
 * Create Bandpass table, creates entries at time and for antennas numbers
 * up to the max in the data.
 * \param inData  UV onto which the input BP table is attached
 * \param outData UV onto which the BP table is to be attached
 * \param time    Time in days to write entry
 * \param bchan   First channel calibrated
 * \param echan   Highest channel calibrated
 * \param nif     Number of IFs
 * \param phase   Channel phase in radians
 * \param snr     Channel snr
 * \param refAnt  Reference antenna
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* UpdateBPTable (ObitUV* inData, olong BPver, ObitUV* outData, 
			    ofloat time, olong bchan, olong echan, olong nif, 
			    ofloat *phase,  ofloat *snr, olong refAnt, ObitErr *err)
{
  ObitTableBP *BPOut=NULL, *BPIn=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  ObitUVDesc *desc;
  olong irow, orow, nchan, npol, ver, ichan, iif, nc;
  olong indx, jndx;
  ofloat xr, xi, yr, yi, fblank = ObitMagicF();
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "MakeBPTable";

  /* error checks */
  if (err->error) return BPOut;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Table info - channels */
  desc = (ObitUVDesc*)outData->myIO->myDesc;
  if (desc->jlocf>=0)
    nchan = desc->inaxes[desc->jlocf];
  else
    nchan = 1;  /* Curious it */
  /* Stokes */
  if (desc->jlocs>=0)
    npol = MIN (2, desc->inaxes[desc->jlocs]);
  else
    npol = 1; 

  /* Create input BP table */
  ver = BPver;
  ObitInfoListGetTest(myInput, "BPVer",  &type, dim, &ver);
  BPIn = newObitTableBPValue ("Temp BP", (ObitData*)inData, &ver,
			      OBIT_IO_ReadOnly, npol, nif, nchan, err);
  /* Does it exist? */
  if (BPIn==NULL) 
    Obit_log_error(err, OBIT_Error, "Failure finding input BP ver %d", BPver);
  if (err->error) Obit_traceback_val (err, routine, inData->name, BPOut);
  
  /* Create output BP table */
  ver = 0;
  ObitInfoListGetTest(myInput, "BPSoln",  &type, dim, &ver);
  BPOut = newObitTableBPValue ("Temp BP", (ObitData*)outData, &ver,
			       OBIT_IO_WriteOnly, npol, nif, nchan, err);
  if (err->error) Obit_traceback_val (err, routine,outData->name, BPOut);
  
  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)BPOut, err);
  
  /* Open BP tables */
  retCode = ObitTableBPOpen (BPIn, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPIn->name, BPOut);
  retCode = ObitTableBPOpen (BPOut, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = newObitTableBPRow(BPOut);
  ObitTableBPSetRow (BPOut, BPRow, err);
  if (err->error) Obit_traceback_val (err, routine,outData->name, BPOut);

  /* Set header values */
  BPOut->numAnt    = BPIn->numAnt;  /* Max. antenna number */
  BPOut->numShifts = BPIn->numShifts;
  BPOut->numChan   = BPIn->numChan;
  BPOut->startChan = BPIn->startChan;
  BPOut->lowShift  = BPIn->lowShift;
  BPOut->shiftInc  = BPIn->shiftInc;
  strncpy (BPOut->BPType, "          ", MAXKEYCHARTABLEBP);
  BPOut->myDesc->sort[0] = BPIn->myDesc->sort[0];  /* Sort order */
  BPOut->myDesc->sort[1] = BPIn->myDesc->sort[1];

  /* reset now invalid pointer to IO descriptor */
  desc = (ObitUVDesc*)outData->myIO->myDesc;

  /* Loop over input table updating */
  for (irow=1; irow<=BPIn->myDesc->nrow; irow++) {
    retCode = ObitTableBPReadRow (BPIn, irow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Fill actual values */
    BPRow->RefAnt1 = refAnt;
    BPRow->RefAnt2 = refAnt;
    nc      =  (echan-bchan+1);
    for (iif=0; iif<nif; iif++) {
      for (ichan=bchan; ichan<=echan; ichan++) {
	indx = iif*nc + ichan-bchan;
	jndx = iif*nchan + ichan -1;
	if ((phase[indx]!=fblank) && (BPRow->Real2[jndx]!=fblank) && 
	    (BPRow->Imag2[jndx]!=fblank))  {
	  if (npol>1) {   /* There better be or this is pretty pointless */
	    yr = BPRow->Real2[jndx];
	    yi = BPRow->Imag2[jndx];
	    xr = cos(phase[indx]);
	    xi = sin(phase[indx]);
	    BPRow->Real2[jndx]   = xr*yr - xi*yi;
	    BPRow->Imag2[jndx]   = xr*yi + xi*yr;
	    BPRow->Weight2[jndx] = snr[indx];
	  }
	} else {  /* No value - blank 2nd poln */
	  if (npol>1) { 
	    BPRow->Real2[jndx]   = fblank;
	    BPRow->Imag2[jndx]   = fblank;
	    BPRow->Weight2[jndx] = 0.0;
	  }
	}
      } /* end loop over channel */
    } /* end loop over IF */
    
    /* Write output table */
    orow = -1;
    retCode = ObitTableBPWriteRow (BPOut, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over antennas */

 cleanup:

  /* Close BP table */
  retCode = ObitTableBPClose (BPOut, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  retCode = ObitTableBPClose (BPIn, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPIn->name, BPOut);
  BPRow = ObitTableBPRowUnref(BPRow);
  BPIn  = ObitTableBPUnref(BPIn);

  return BPOut;
} /* end UpdateBPTable  */

/**
 * Get source lookup table, translates Source IDs to index in 
 * Source ordered arrays.
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
 * Setup EVPA sin/cos arrays for a source, how to flip signs of data
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
