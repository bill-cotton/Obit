/* $Id$  */
/* Instrumental polarization calibration                              */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012,2014                                          */
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
#include "ObitPolnCalFit.h"
#include "ObitComplex.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* PCalIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void PCalOut (ObitInfoList* outList, ObitErr *err);
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
/* Get input sky model */
ObitSkyModel* getInputSkyModel (ObitInfoList *myInput, ObitErr *err);
/* Write history */
void PCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Average data */
ObitUV* AverData (ObitInfoList* myInput, ObitUV* scrData, 
		  ObitPolnCalFit *PolnFitter, ObitErr* err);


/* Program globals */
gchar *pgmName = "PCal";      /* Program name */
gchar *infile  = "PCal.in" ;  /* File with program inputs */
gchar *outfile = "PCal.out";  /* File to contain program outputs */
olong  pgmNumber;               /* Program number (like POPS no.) */
olong  AIPSuser;                /* AIPS user number number  */
olong  nAIPS=0;                 /* Number of AIPS directories */
gchar **AIPSdirs=NULL;          /* List of AIPS data directories */
olong  nFITS=0;                 /* Number of FITS directories */
gchar **FITSdirs=NULL;          /* List of FITS data directories */
ObitInfoList *myInput  = NULL;  /* Input parameter list */
ObitInfoList *myOutput = NULL;  /* Output parameter list */
olong  souNo=0;                 /* Single source number */
olong  nCal=1;                  /* Number of calibrators selected */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*    Instrumental polarization calibration                               */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *scrData = NULL, *avgData = NULL;;
  ObitSkyModel *skyModel=NULL;
  ObitPolnCalFit *PolnFitter=NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = PCalIn (argc, argv, err);
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

  /* Get input sky model */
  skyModel = getInputSkyModel (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Copy/select/calibrate to scratch file */
  scrData = newObitUVScratch (inData,err);
  scrData = ObitUVCopy (inData, scrData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  /* Save number of selected calibrators */
  nCal = MAX(1, inData->mySel->numberSourcesList);
  
  /* Save first source id in case it's the only one */
  if (inData->mySel->sources) souNo = inData->mySel->sources[0];
  
  /* Divide if model given */
  if (skyModel) ObitSkyModelDivUV (skyModel, scrData, scrData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create Fitter */
  PolnFitter = ObitPolnCalFitCreate ("PolnFitter");

  /* Average data */
  avgData = AverData (myInput, scrData, PolnFitter, err);

  /* Do channel solutions */
  ObitPolnCalFitFit (PolnFitter, avgData, inData, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Write history */
  PCalHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* Shutdown Obit */
 exit: 
  /* Remove KeepSou keyword on imput descriptor */
  ObitUVOpen(inData, OBIT_IO_ReadWrite, err);
  ObitInfoListRemove (inData->myDesc->info, "KeepSou");
  ObitInfoListRemove (((ObitUVDesc*)inData->myIO->myDesc)->info, "KeepSou");
  inData->myStatus = OBIT_Modified;
  ObitUVClose(inData,  err);
  /* cleanup */
  PolnFitter = ObitPolnCalFitUnref(PolnFitter); 
  myInput    = ObitInfoListUnref(myInput); 
  inData     = ObitUnref(inData);
  scrData    = ObitUnref(scrData);
  avgData    = ObitUnref(avgData);
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* PCalIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "PCalIn";

  /* error checks */
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
    if (strcmp(arg, "-output") == 0){ /* output results */
      outfile = argv[++ax];
      ObitReturnDumpRetCode (-999, outfile, myOutput, err);
      if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

    } else if (strcmp(arg, "-input") == 0){ /* input parameters */
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

  return list;
} /* end PCalIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: PCal -input file -output ofile [args]\n");
    fprintf(stderr, "PCal Instrumental polarization calibration for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def PCal.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def PCal.out\n");
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
  gchar ctemp[5];
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
  strTemp = "PCal.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "PCalName";
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
  
  /* No average IFs in initial calibration */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "avgIF", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* No average Poln in  initial calibration*/
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "avgPol", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* A&P solution in  initial calibration*/
  dim[0] = 4; dim[1] = 1;
  ctemp[0]='A'; ctemp[1]='&'; ctemp[2]='P'; ctemp[3]=' '; ctemp[4]=0; 
  ObitInfoListPut (out, "solType", OBIT_bool, dim, ctemp, err);
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *strTemp;
  ObitSkyModelMode modelMode;
  ObitSkyModelType modelType;
  ofloat       modelFlux;
  gboolean     doCalSelect;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Convert test Cmethod to enum  Mode */
  ObitInfoListGetP (myInput, "Cmethod", &type, dim, (gpointer)&strTemp);
  if (!strncmp (strTemp, "GRID", 4)) modelMode = OBIT_SkyModel_Grid;
  else if (!strncmp (strTemp, "DFT", 3)) modelMode = OBIT_SkyModel_DFT;
  else modelMode = OBIT_SkyModel_Fastest;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "Mode", OBIT_long, dim, &modelMode);

  /* Convert test Cmodel to enum  ModelType */
  ObitInfoListGetP (myInput, "Cmodel", &type, dim, (gpointer)&strTemp);
  modelFlux = 0.0;
  ObitInfoListGetTest (myInput, "modelFlux", &type, dim, &modelFlux); 
  if (!strncmp (strTemp, "COMP", 4)) modelType = OBIT_SkyModel_Comps;
  else if (!strncmp (strTemp, "IMAG", 3)) modelType = OBIT_SkyModel_Image;
  else modelType = OBIT_SkyModel_Comps;
  /* Is a model given in the parameters? */
  if (modelFlux!=0.0)  modelType = OBIT_SkyModel_Point;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "ModelType", OBIT_long, dim, &modelType);

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
  gboolean     doCalSelect, KeepSou;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "UVRange", "FreqID", 
    "BChan", "EChan", "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "doPol", "Mode", "ModelType", "Alpha", "refAnt",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Make sure multisource nature kept even if only selecting one source */
  KeepSou = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (myInput, "inKeepSou", OBIT_bool, dim, &KeepSou);

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
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
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Make sure multisource nature kept even if only selecting one source */
  KeepSou = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (inData->myDesc->info, "KeepSou", OBIT_bool, dim, &KeepSou);
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get input sky model                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*      Sky Model to be used, NULL if none given                          */
/*----------------------------------------------------------------------- */
ObitSkyModel* getInputSkyModel (ObitInfoList *myInput, ObitErr *err)
{
  ObitSkyModel *skyModel=NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitImage    **image=NULL;
  ObitInfoType type;
  gboolean     do3D=TRUE;
  olong         Aseq, disk, cno,i=0, nmaps;
  gchar        *Type, *Type2, *strTemp, inFile[129], inRoot[129];
  gchar        Aname[13], Aclass[7], Aroot[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat       smodel[20], modptflx,  modptxof, modptyof, modptypm[4];
  gchar        name[101];
  gchar        *dataParms[] = {  /* Control parameters */
    "CCVer",  "BComp",  "EComp",  "Flux", "PBCor", "antSize", "Factor", 
    "minFlux", "Mode", "ModelType", "REPLACE", "Stokes", 
    "MODPTFLX", "MODPTXOF", "MODPTYOF", "MODPTYPM", 
    NULL};
  gchar *routine = "getInputSkyModel";

  /* error checks */
  if (err->error) return skyModel;
  g_assert (ObitInfoListIsA(myInput));

  /* How many fields? */
  nmaps = 1;
  ObitInfoListGetTest(myInput, "nfield", &type, dim, &nmaps);

  /* Image model of model in parameters? */
  smodel[0] = 0.0;
  ObitInfoListGetTest (myInput, "modelFlux", &type, dim, &smodel[0]);
  if ((smodel!=NULL) && (smodel[0]!=0.0)) {
    /* Model passed - get rest */
    ObitInfoListGetTest (myInput, "modelPos",  &type, dim, &smodel[1]);
    ObitInfoListGetTest (myInput, "modelParm", &type, dim, &smodel[3]);
    modptflx = smodel[0];
    modptxof = smodel[1] / 3600.0;
    modptyof = smodel[2] / 3600.0;
    modptypm[0] = smodel[3];
    modptypm[1] = smodel[4];
    modptypm[2] = smodel[5];
    modptypm[3] = smodel[6];
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "MODPTFLX", OBIT_float, dim, &modptflx);
    ObitInfoListAlwaysPut (myInput, "MODPTXOF", OBIT_float, dim, &modptxof);
    ObitInfoListAlwaysPut (myInput, "MODPTYOF", OBIT_float, dim, &modptyof);
    dim[0] = 4;
    ObitInfoListAlwaysPut (myInput, "MODPTYPM", OBIT_float, dim, modptypm);

    /* Create Sky Model */
    skyModel = newObitSkyModel ("Sky Model");

    Obit_log_error(err, OBIT_InfoErr, "Using input model parameters");
  } else if (nmaps>0) {  /* Image given */
    
    Obit_log_error(err, OBIT_InfoErr, "Using image model");

    /* Allocate Image array */
    image = g_malloc0(nmaps*sizeof(ObitImage));
    
    /* Create image mosaic */
    mosaic = newObitImageMosaic ("Mosaic", nmaps);
    
    /* File type - could be either AIPS or FITS use DataType2 (default DataType) */
    ObitInfoListGetP (myInput, "DataType",  &type, dim, (gpointer)&Type);
    ObitInfoListGetP (myInput, "DataType2", &type, dim, (gpointer)&Type2);
    if (!strncmp (Type2, "    ", 4)) Type2 = Type;
    if (!strncmp (Type2, "AIPS", 4)) { /* AIPS input */
      /* input AIPS disk */
      ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);
      /* input AIPS name */
      if (ObitInfoListGetP(myInput, "in2Name", &type, dim, (gpointer)&strTemp)) {
	strncpy (Aname, strTemp, 13);
      } else { /* Didn't find */
	strncpy (Aname, "No Name ", 13);
      } 
      Aname[12] = 0;
      /* input AIPS class */
      if  (ObitInfoListGetP(myInput, "in2Class", &type, dim, (gpointer)&strTemp)) {
	strncpy (Aroot, strTemp, 7);
      } else { /* Didn't find */
	strncpy (Aroot, "NoClas", 7);
      }

      /* input AIPS sequence */
      ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);
      
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	/* If only one field use class given */
	if ((nmaps==1) && (i==0)) {
	} else { /* derive class from field number */
	  Aroot[2] = 0;
	  g_snprintf (Aclass, 7, "%s%4.4d",Aroot,i+1);
	} /* end one or many fields */
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
      }
      
      /* Loop over fields */
      for (i=0; i<nmaps; i++) {
	g_snprintf (name, 100, "Input image %d",i+1);
	/* If only one field use class given */
	if ((nmaps==1) && (i==0)) {
	  g_snprintf (Aclass, 7, "%s",Aroot);
	} else { /* derive class from field number */
	  Aroot[2] = 0;
	  g_snprintf (Aclass, 7, "%s%4.4d",Aroot,i+1);
	} /* end one or many fields */
	
	  /* Find catalog number */
	cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
	if (cno<0) Obit_log_error(err, OBIT_Error, "Failure looking up %s", name);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	/* define object */
	image[i] = newObitImage(name);
	ObitImageSetAIPS(image[i], OBIT_IO_byPlane, disk, cno, AIPSuser,  blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, i, image[i], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
      } /* end loop over fields */
      
      /* get do3D from first image */
      do3D = image[0]->myDesc->do3D;
      
    } else if (!strncmp (Type2, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      if (ObitInfoListGetP(myInput, "in2File", &type, dim, (gpointer)&strTemp)) {
	strncpy (inRoot, strTemp, 128);
      } else { 
	strncpy (inRoot, "No_Filename_Given", 128);
      }
      ObitTrimTrail(inRoot);  /* remove trailing blanks */
   
      /* input FITS disk */
      ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);

      if (nmaps==1) {
 	
	/* Set file name */
	g_snprintf (inFile, 128, "%s",inRoot);
 
	/* define object */
	g_snprintf (name, 100, "Input image");
	image[0] = newObitImage(name);
	ObitImageSetFITS(image[0], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	
	  /* Attach Image */
	  ObitImageMosaicSetImage (mosaic, 0, image[0], err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
     } else { /* Multiple fields */
	
	/* Loop over fields */
	for (i=0; i<nmaps; i++) {
	  /* Set file name */
	  g_snprintf (inFile, 128, "%s%d",inRoot,i);
	  
	  /* define object */
	  g_snprintf (name, 100, "Input image %d",i+1);
	  image[i] = newObitImage(name);
	  ObitImageSetFITS(image[i], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	  
	  /* Attach Image */
	  ObitImageMosaicSetImage (mosaic, i, image[i], err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
	} /* end loop over fields */
      }

      /* get do3D from first image */
      do3D = image[0]->myDesc->do3D;
      
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		     pgmName, Type2);
      return skyModel;
    }

    /* Create Sky Model */
    skyModel = ObitSkyModelCreate ("Sky Model", mosaic);

    /* deallocate images */
    for (i=0; i<nmaps; i++) image[i] = ObitImageUnref(image[i]);
    g_free(image);  /* Deallocate array */

   /* End image or components model */
  } else {  /* Neither model given just return */
    Obit_log_error(err, OBIT_InfoErr, 
		   "Using source fluxes from SU table for point model");
    return skyModel;
  }
 
  /* Get input parameters from myInput, copy to skyModel */
  ObitInfoListCopyList (myInput, skyModel->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, skyModel->name, skyModel);
  
  /* Save do3D */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut (skyModel->info, "do3D", OBIT_bool, dim, &do3D);
  
  return skyModel;
} /* end getInputSkyModel */
  
/*----------------------------------------------------------------------- */
/*  Write History for PCal                                                */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void PCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "UVRange", "subA",
    "selBand", "selFreq", "FreqID",  "BChan", "EChan", "BIF", "EIF",  
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol", "Antennas", 
    "DataType2", "in2File", "in2Disk", "in2Name", "in2Class", "in2Seq", 
    "nfield", "CCVer", "BComp", "EComp", "Cmethod", "Cmodel", "Flux",
    "modelFlux", "modelPos", "modelParm", "solInt", 
    "solnType", "Sources", "Qual", "souCode", "doFitI", "doFitPol", "doFitV", 
    "doFitGain", "doFitRL", "RLPhase", "RM", "PPol", "ChWid", "ChInc", 
    "CPSoln", "PDSoln", "BPSoln", "doBlank", 
    "nThreads",
   NULL};
  gchar *routine = "PCalHistory";

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
 
} /* end PCalHistory  */

/**
 * Averages data to solint and returns averaged data.
 * \param myInput      Input parameters on InfoList    
 * \param scrUV        Data to be used to determine calibration
 * \param PolnFitter   Fitter to be used, control parameters copied to info
 * \param err          ObitErr stack for reporting problems.
 * \return time averaged visibility data
 */
ObitUV* AverData (ObitInfoList* myInput, ObitUV* scrData, 
		  ObitPolnCalFit *PolnFitter, ObitErr* err)
{
  ObitUV       *avgData=NULL;
  ofloat       solInt;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *PolCalParms[] = {     /* Polarization calibration parameters */
    "solnType", "Sources", "Qual", "souCode", "doFitI", "doFitPol", "doFitV", 
    "doFitGn", "doFitRL", "RLPhase", "RM", "PPol", "ChWid", "ChInc", 
    "CPSoln", "PDSoln", "BPSoln",  "doBlank", "doFitGain", 
    "doBand", "BPVer", "refAnt", "prtLv", "BIF", "BChan",
    NULL};
  gchar *FQInclude[] = {"AIPS FQ", "AIPS AN", NULL};
  gchar *routine = "AverData";

  /* error checks */
  if (err->error) return avgData;
  g_assert (ObitUVIsA(scrData));
  
  /* Average data to solInt; 0=> all */
  solInt = 0.0;
  ObitInfoListGetTest (myInput, "solInt", &type, dim, &solInt); 
  if (solInt<=1.0e-5) solInt = 1000.0; 
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "timeAvg", OBIT_float, dim, &solInt);
  avgData = ObitUVUtilAvgT(scrData, TRUE, avgData, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Be sure FQ table copied */
   ObitUVCopyTables (scrData, avgData, NULL, FQInclude, err);

   /* Copy poln cal calibration control to PolnFitter */
   ObitInfoListCopyList (myInput, PolnFitter->info, PolCalParms);

   /* Set number of calibrators */
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  /* nCal from global */
  ObitInfoListAlwaysPut(PolnFitter->info, "nCal", OBIT_long, dim, &nCal);

  return avgData;
} /* end AverData */
