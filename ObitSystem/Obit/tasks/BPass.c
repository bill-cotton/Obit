/* $Id$  */
/* Obit Radio interferometry calibration software                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009-2011                                          */
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
#include "ObitSkyModelMF.h"
#include "ObitTableSUUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitUVUtil.h"
#include "ObitTableSN.h"
#include "ObitTableCL.h"
#include "ObitTableCLUtil.h"
#include "ObitTableBP.h"
#include "ObitUVSoln2Cal.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* BPassIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void BPassOut (ObitInfoList* outList, ObitErr *err);
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
/* Do channel solutions */
void  BandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		  ObitErr* err);
/* Divide by SU table flux density */
void DivideSource (ObitUV *inUV, ObitUV *scrUV, ObitErr *err);
/* Divide buffer table flux density */
static void DivideBuffer (ObitSourceList *sList, olong sourId, ObitUV *uvdata);
/* Write history */
void BPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Initial calibration */
ObitUV* InitialCal (ObitInfoList* myInput, ObitUV* scrData, ObitErr* err);
/* Create initial BP table */
ObitTableBP* DummyBPTable (ObitUV* inData, ObitTableSN *SNTab, ObitErr *err);
/* Copy SN info to BP table */
void SN2BPTable (ObitTableSN *SNTab, ObitTableBP *BPTab, olong chan,
		 ObitErr *err);
/* Bandpass from AutoCorrelations */
void AutoCorrBP (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		 ObitErr* err);

/* Program globals */
gchar *pgmName = "BPass";       /* Program name */
gchar *infile  = "BPass.in" ;   /* File with program inputs */
gchar *outfile = "BPass.out";   /* File to contain program outputs */
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
/*   Obit task to determine gains for a uv data set                       */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *scrData = NULL, *avgData = NULL;;
  ObitSkyModel *skyModel=NULL;
  ObitErr      *err= NULL;
  ofloat       ftemp;
  gboolean     doAuto = FALSE;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *SCParms[] = {  /* Selfcal Parameters  */
    "refAnt",
     NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = BPassIn (argc, argv, err);
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

   /* Use autocorrelations rather than cross correlations? */
  ObitInfoListGetTest (myInput, "doAuto", &type, dim, &doAuto); 
  if (doAuto) {
    AutoCorrBP (myInput, inData, inData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  } else { /* cross correlations */
    
    /* Get input sky model */
    skyModel = getInputSkyModel (myInput, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    
    /* Copy/select/calibrate to scratch file */
    scrData = newObitUVScratch (inData,err);
    scrData = ObitUVCopy (inData, scrData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    
    /* Save first source id in case it's the only one */
    if (inData->mySel->sources) souNo = inData->mySel->sources[0];
    
    /* Divide if model given */
    if (skyModel) ObitSkyModelDivUV (skyModel, scrData, scrData, err);
    /* Otherwise divide by source flux density if given */
    else DivideSource (inData, scrData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    /* Self cal parameters to scrData */
    ObitInfoListCopyList (myInput, scrData->info, SCParms);
  
    /* Index scrData */
    dim[0] = dim[1] = 1;
    ftemp = 15.0;  /* Max scan time 15 min. */
    ObitInfoListAlwaysPut(scrData->info, "maxScan", OBIT_float, dim, &ftemp);
    ftemp = 1.0; /* Max allowable gap 1 min. */
    ObitInfoListAlwaysPut(scrData->info, "maxGap", OBIT_float, dim, &ftemp);
    ObitUVUtilIndex (scrData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    
    /* Initial calibration */
    avgData =  InitialCal(myInput, scrData, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    
    /* Do channel solutions, convert to BP table */
    BandpassCal(myInput, avgData, inData, err);
    if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;
  } /* end calibration */

  /* Write history */
  BPassHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
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

ObitInfoList* BPassIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "BPassIn";

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
} /* end BPassIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: BPass -input file -output ofile [args]\n");
    fprintf(stderr, "BPass Obit task determine bandpass for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def BPass.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def BPass.out\n");
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
  strTemp = "BPass.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "BPassName";
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

  /* Apply PBCor?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "PBCor", OBIT_bool, dim, &btemp, err);
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
  olong        nvis, nThreads, doCalib;
  gboolean     doCalSelect;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "FreqID", "BChan", "EChan",   "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "doPol", "Mode", "ModelType", "Alpha", "refAnt",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1000;
  nThreads = 1;
  ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
  nvis *= nThreads;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  type, dim,  &nvis);
    
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

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

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
  ObitCCCompType CCType;
  ObitInfoType  type;
  gboolean      do3D=TRUE;
  olong         Aseq, disk, cno,i=0, nmaps, ver;
  gchar        *Type, *Type2, *strTemp, inFile[129], inRoot[129];
  gchar         Aname[13], Aclass[7], Aroot[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat        smodel[20], modptflx,  modptxof, modptyof, modptypm[4];
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

    /* Create Sky Model for appropriate type */
    ver = 0;
    ObitInfoListGetTest(myInput, "CCVer", &type, dim, &ver);
    CCType   = ObitTableCCUtilGetType ((ObitData*)mosaic->images[0], ver, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", skyModel);
    if ((CCType==OBIT_CC_PointModTSpec)|| (CCType==OBIT_CC_GaussModTSpec) ||
	(CCType==OBIT_CC_CGaussModTSpec) || (CCType==OBIT_CC_USphereModTSpec)) {
      skyModel = (ObitSkyModel*)ObitSkyModelMFCreate ("Sky Model", mosaic);
      Obit_log_error(err, OBIT_InfoErr, "Using tabulated spectrum sky model");
    } else
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
  
/**
 * Divide model into uvdata in UVin and rewrite
 * \param inUV    Input UV data, possibly with SU table
 * \param scrUV   Uv data object to be corrected.
 * \param err     ObitErr stack for reporting problems.
 */
void DivideSource (ObitUV *inUV, ObitUV *scrUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableSU *SUTable = NULL;
  ObitSourceList *sList=NULL;
  olong firstVis, sid, sourId, iver, i, j;
  gchar *routine = "DivideSource";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(scrUV));

  /* See if there is a source table */
  iver = 1;
  SUTable = newObitTableSUValue (inUV->name, (ObitData*)inUV, &iver, OBIT_IO_ReadOnly, 0, err);
  if (err->error) return;
  if (SUTable==NULL) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: NO source table found on %s", routine, inUV->name);
    return;
  }

  sList = ObitTableSUGetList (SUTable, err);
  if (err->error) goto cleanup;

  /* index in sList if only one source selected */
  if (inUV->mySel->sources) sid    = inUV->mySel->sources[0];
  else sid = 0;
  sourId = sid;

  /* Convert flux densities to reciprocal, 0=>1 */
  for (i=0; i<sList->number; i++) {
    if (sList->SUlist[i]->SourID == sid) sourId = i;
    for (j=0; j<sList->SUlist[i]->numIF; j++) {
      if (sList->SUlist[i]->IFlux[j] == 0.0) sList->SUlist[i]->IFlux[j] = 1.0;
      sList->SUlist[i]->IFlux[j] = 1.0 / sList->SUlist[i]->IFlux[j];
    }
  }

  retCode = ObitUVOpen (scrUV, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* loop correcting data */
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitUVRead (scrUV, NULL, err);
    if (retCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) goto cleanup;
    firstVis = scrUV->myDesc->firstVis;
    
    /* Divide buffer load */
    DivideBuffer (sList, sourId, scrUV);

    /* rewrite buffer */
    retCode = ObitUVWrite (scrUV, NULL, err);
    if (err->error) goto cleanup;
    scrUV->myDesc->firstVis = firstVis;  /* reset first vis in buffer */
    ((ObitUVDesc*)scrUV->myIO->myDesc)->firstVis = firstVis;
    
  } /* end loop dividing data */

  /* Cleanup */
 cleanup:
  /* Close data */
  retCode = ObitUVClose (scrUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  if (SUTable) ObitTableSUUnref(SUTable);
  if (sList) ObitSourceListUnref(sList);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

} /* end DivideSource  */

/**
 * Divide buffer load of data by reciprocal of flux density in SourceList
 * \param sList   Source list with 1/flux densities
 * \param sourId  Source ID (0-rel) in sList is uvdata a siingle source file
 * \param uvdata  Object with uv data in buffer, prepared for correcrion.
 */
static void DivideBuffer (ObitSourceList *sList, olong sourId, ObitUV *uvdata)
{
  olong ivis, nvis, ifreq, nfreq, iif, nif, nstok, istok, sid, i;
  olong sour, lastsour;
  ofloat *vis, *ifvis, *fvis, *svis, *rec, factor;
  ObitUVDesc *desc;

  /* how much data? */
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  nstok = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  sid = sourId;
  lastsour = -10;
 
  /* Loop over visibilities in buffer */
  vis = uvdata->buffer+desc->nrparm;
  rec = uvdata->buffer;
  for (ivis=0; ivis<nvis; ivis++) {

    /* Which source? */
    if (desc->ilocsu>=0) {
      sour = (olong)(*(rec+desc->ilocsu) + 0.5);
      if (lastsour!=sour) {
	lastsour = sour;
	/* find in sList */
	for (i=0; i<sList->number; i++) {
	  if (sList->SUlist[i]->SourID == sour) {sid=i; break;}
	}
      }
    }
    
    /* loop over IFs */
    ifvis = vis;
    for (iif=0; iif<nif; iif++) {
      factor = sList->SUlist[sid]->IFlux[iif];
      
      /* loop over frequencies */
      fvis = ifvis;
      for (ifreq = 0; ifreq<nfreq; ifreq++) {

	/* Loop over Stokes */
	svis = fvis;
	for (istok=0; istok<nstok; istok++) {
	  
	  /* divide */
	  *(svis)   *= factor;
	  *(svis+1) *= factor;
	  
	  svis += desc->incs; /* visibility pointer */
	} /* end loop over Stokes */
	
	  fvis += desc->incf; /* visibility pointer */
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */

    /* update data pointers */
    vis += desc->lrec;
    rec += desc->lrec;
  } /* end loop over visibilities */
} /* end DivideBuffer */

/*----------------------------------------------------------------------- */
/*  Write History for BPass                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void BPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA",
    "selBand", "selFreq", "FreqID", 
    "BChan1", "EChan1", "BChan2", "EChan2", "ChWid2", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol", "Antennas",  "refAnt",  "doAuto",
    "DataType2", "in2File", "in2Disk", "in2Name", "in2Class", "in2Seq", 
    "nfield", "CCVer", "BComp", "EComp", "Cmethod", "Cmodel", "Flux",
    "modelFlux", "modelPos", "modelParm", "Alpha",
    "solInt1", "solInt2", "solType", "solMode", "avgPol", "avgIF", "doMGM", "minSNR",
    "minNo", "prtLv",
    "nThreads",
   NULL};
  gchar *routine = "BPassHistory";

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
 
} /* end BPassHistory  */

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
  gboolean     btemp;
  ofloat       ftemp, solInt;
  olong        itemp, ver=1;
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
  gchar *FQInclude[] = {"AIPS FQ", NULL};
  gchar *routine = "InitialCal";

  /* error checks */
  if (err->error) return avgData;
  g_assert (ObitUVIsA(scrData));
  
  /* Initial CL table on scratch file */
  ObitInfoListGet(myInput, "solInt1", &type, dim, &ftemp, err);
  ObitInfoListAlwaysPut(scrData->info, "solInt", type, dim, &ftemp);
  CLTable1 =  ObitTableCLGetDummy (scrData, scrData, ver, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Create solver */
  solver = ObitUVGSolveCreate("Gain solver");

  /* Copy calibration control to solver */
  ObitInfoListCopyList (myInput, solver->info, solverParms);
  ObitInfoListGet(myInput, "solInt1", &type, dim, &ftemp, err);
  ObitInfoListAlwaysPut(solver->info, "solInt", OBIT_long, dim, &ftemp);
  itemp = 1;
  ObitInfoListAlwaysPut(solver->info, "solnVer", OBIT_long, dim, &itemp);
  ObitInfoListGet(myInput, "BChan1", &type, dim, &itemp, err);
  ObitInfoListAlwaysPut(scrData->info, "BChan", OBIT_long, dim, &itemp);
  ObitInfoListGet(myInput, "EChan1", &type, dim, &itemp, err);
  ObitInfoListAlwaysPut(scrData->info, "EChan", OBIT_long, dim, &itemp);
  dim[0] = strlen(blank);  /* Phase only */
  ObitInfoListAlwaysPut(solver->info, "solMode", OBIT_string, dim, blank);


  /* Do atmospheric solution */
  SNTable = ObitUVGSolveCal (solver, scrData, scrData, scrData->mySel, err);
  if (err->error) Obit_traceback_val (err, routine, scrData->name, avgData);

  /* Apply calibration to CL table */
  dim[0] = dim[1], dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(scrData->info, "solnVer", OBIT_long, dim, &SNTable->tabVer);
  ObitInfoListAlwaysPut(scrData->info, "calIn",   OBIT_long, dim, &CLTable1->tabVer);
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
  ftemp = 0.0;
  ObitInfoListGetTest(myInput, "Alpha", &type, dim, &ftemp);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut (scrData->info, "Alpha", OBIT_float, dim, &ftemp);

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
void  BandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		  ObitErr* err)
{
  ObitTableSN **SNTables = NULL;
  ObitTableBP *BPTable = NULL;
  ObitUVGSolve *solver=NULL;
  olong ichan, nchan, bchan2=0, echan2=0, chinc2=0;
  olong nif, maxch, maxno, itemp, highVer, ver;
  gboolean btemp;
  ofloat ftemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar        *solverParms[] = {  /* Calibration parameters */
    "solnVer", "solType", "solMode", "avgPol", "avgIF", "doMGM", "elevMGM",
    "refAnt", "ampScalar", "minSNR",  "minNo", "prtLv",
    NULL};
  gchar *copyBPTable[] = {"AIPS BP", NULL};
  gchar *routine = "BandpassCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(avgData));

  /* Range of channels */
  nchan = avgData->myDesc->inaxes[avgData->myDesc->jlocf];
  ObitInfoListGetTest (myInput, "BChan2", &type, dim, &bchan2); 
  if (bchan2<=0) bchan2 = 1;
  ObitInfoListGetTest (myInput, "EChan2", &type, dim, &echan2); 
  if (echan2<=0) echan2 = nchan;
  echan2 = MIN (echan2, nchan);
  /* How many channels to use at a time */
  ObitInfoListGetTest (myInput, "ChWid2", &type, dim, &chinc2); 
  chinc2 /= 2;
  chinc2 = MAX (0, chinc2);

  /* Number of IFs */
  if (avgData->myDesc->jlocif>=0)  
    nif = avgData->myDesc->inaxes[avgData->myDesc->jlocif];
  else
    nif = 1;

  /* Array of SN tables */
  SNTables = g_malloc0(nchan*sizeof(*SNTables));

  /* Create solver */
  solver = ObitUVGSolveCreate("Gain solver");
  /* Save first source id in case it's the only one */
  if (inData->mySel->sources) solver->curSource = inData->mySel->sources[0];

  /* Copy calibration parameters */
  ObitInfoListCopyList (myInput, solver->info, solverParms);
  ObitInfoListGet(myInput, "solInt2", &type, dim, &ftemp, err);
  ObitInfoListAlwaysPut(solver->info, "solInt", type, dim, &ftemp);

  /* Don't need to apply calibration */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(avgData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = -1;
  ObitInfoListAlwaysPut(avgData->info, "doCalib", OBIT_long, dim, &itemp);
  /* Don't care how many work */
  ftemp = 0.0;
  ObitInfoListAlwaysPut(solver->info, "minOK", OBIT_float, dim, &ftemp);

  /* Loop over channels getting channel calibration */
  for (ichan=bchan2; ichan<=echan2; ichan++) {
    /* Select channel */
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    itemp = MAX (1, ichan-chinc2);
    ObitInfoListAlwaysPut(avgData->info, "BChan", OBIT_long, dim, &itemp);
    itemp = MIN (nchan, ichan+chinc2);
    ObitInfoListAlwaysPut(avgData->info, "EChan", OBIT_long, dim, &itemp);
    
    /* Solve one channel and all IFs */
    SNTables[ichan-1] = ObitUVGSolveCal (solver, avgData, avgData, avgData->mySel, err);
    if (err->error) Obit_traceback_msg (err, routine, avgData->name);
  }  /* end channel loop */

  /* Find channel with maximum number of solutions - 
     use it to create initial BP table */
  maxch = bchan2; 
  maxno = 0;
  for (ichan=bchan2; ichan<=echan2; ichan++) {
    if (SNTables[ichan-1] && (SNTables[ichan-1]->myDesc->nrow>maxno)) {
      maxch = ichan; 
      maxno = SNTables[ichan-1]->myDesc->nrow;
    }
  }
  /* Create BP table */
  BPTable = DummyBPTable (avgData, SNTables[maxch-1], err);

  /* Loop over SN tables copying to BP */
  for (ichan=bchan2; ichan<=echan2; ichan++) {
    SN2BPTable (SNTables[ichan-1], BPTable, ichan-1, err);
    if (err->error) Obit_traceback_msg (err, routine, avgData->name);
  } /* end  end loop copying */  

  /* Copy BP table to inData */
  ObitUVCopyTables (avgData, inData, NULL, copyBPTable, err);
  if (err->error) Obit_traceback_msg (err, routine, avgData->name);

  /* Tell which BP version */
  ver = 0;
  ObitInfoListGetTest(myInput, "BPSoln",  &type, dim, &ver);
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS BP");
  if (ver<=0) ver = highVer;
  Obit_log_error(err, OBIT_InfoErr, "Writing crosscorrelation BP Table %d", ver);

  /* Cleanup */
  solver = ObitUVGSolveUnref(solver);
  ObitUVZapTable (avgData, "AIPS BP", -1, err);
  BPTable = ObitTableBPUnref(BPTable );
  ObitUVZapTable (avgData, "AIPS SN", -1, err);
  if (err->error) Obit_traceback_msg (err, routine, avgData->name);
  for (ichan=0; ichan<nchan; ichan++) 
    SNTables[ichan] = ObitTableSNUnref(SNTables[ichan]);

} /* end BandpassCal  */

/**
 * Create dummy Bandpass table, creates entries at times and for antennas
 * of the entries in SNTmpl, bandpass calibration blanked.
 * \param inData  UV onto which the BP table is to be attached
 * \param SNTmpl  Template SN table
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* DummyBPTable (ObitUV* inData, ObitTableSN *SNTmpl, ObitErr *err)
{
  ObitTableBP *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableSNRow *SNRow=NULL;
  ObitTableBPRow *BPRow=NULL;
  ObitUVDesc *desc;
  olong i, irow, orow, nchan, nif, npol, ver;
  ofloat fblank = ObitMagicF();
  ObitInfoType type;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "DummyBPTable";

  /* error checks */
  if (err->error) return BPOut;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitTableSNIsA(SNTmpl));

  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine,inData->name, BPOut);

  /* Table info - channels */
  desc = (ObitUVDesc*)inData->myIO->myDesc;
  if (desc->jlocf>=0)
    nchan = desc->inaxes[desc->jlocf];
  else
    nchan = 1;  /* Curious */
  /* IFs */
  if (desc->jlocif>=0)
    nif = desc->inaxes[desc->jlocif];
  else
    nif = 1;
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

  /* Open SN table */
  retCode = ObitTableOpen ((ObitTable*)SNTmpl, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, SNTmpl->name, BPOut);
  SNRow = newObitTableSNRow(SNTmpl);

  /* Open BP table */
  retCode = ObitTableOpen ((ObitTable*)BPOut, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = newObitTableBPRow(BPOut);
  ObitTableBPSetRow (BPOut, BPRow, err);
  if (err->error) Obit_traceback_val (err, routine,inData->name, BPOut);

  /* Set header values */
  BPOut->numAnt    = SNTmpl->numAnt;  /* Max. antenna number */
  BPOut->numShifts = 0;
  BPOut->numChan   = nchan;
  BPOut->startChan = 1;
  BPOut->lowShift  = 1;
  BPOut->shiftInc  = 1;
  strncpy (BPOut->BPType, "          ", MAXKEYCHARTABLEBP);
  if ((SNTmpl->myDesc->sort[0]-1)==SNTmpl->TimeCol)
    BPOut->myDesc->sort[0] = BPOut->TimeCol+1;  /* Sort order */
  if ((SNTmpl->myDesc->sort[1]-1)==SNTmpl->antNoCol)
    BPOut->myDesc->sort[1] = BPOut->antNoCol+1;

  /* Initialize BP Row */
  desc = (ObitUVDesc*)inData->myIO->myDesc;
  BPRow->BW           = desc->cdelt[desc->jlocf];
  BPRow->ChanShift[0] = 0.0;
  BPRow->ChanShift[1] = 0.0;
  BPRow->RefAnt1      = 0;
  BPRow->RefAnt2      = 0;
  for (i=0; i<nif; i++) BPRow->ChanShift[i] = 0.0;
  for (i=0; i<nif; i++) BPRow->Weight1[i]   = 0.0;
  if (npol>1) for (i=0; i<nif; i++) BPRow->Weight2[i] = 0.0;
  for (i=0; i<nchan*nif; i++) { 
    BPRow->Real1[i]   = fblank;
    BPRow->Imag1[i]   = fblank;
    if (npol>1) {
      BPRow->Real2[i]   = fblank;
      BPRow->Imag2[i]   = fblank;
    }
  }

  /* Loop over SN table */
  for (irow=1; irow<=SNTmpl->myDesc->nrow; irow++) {
    retCode = ObitTableSNReadRow (SNTmpl, irow, SNRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Set time, antenna etc.*/
    BPRow->Time   = SNRow->Time;
    BPRow->TimeI  = fabs(SNRow->TimeI);
    BPRow->SourID = SNRow->SourID;
    BPRow->antNo  = SNRow->antNo;
    BPRow->SubA   = SNRow->SubA;
    BPRow->FreqID = SNRow->FreqID;
    /* single source in calibration? */
    if (BPRow->SourID==0) BPRow->SourID = souNo;

    /* Write output table */
    orow = -1;
    retCode = ObitTableBPWriteRow (BPOut, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over SN table */

  /* Close SN table */
 cleanup:
  retCode = ObitTableClose ((ObitTable*)SNTmpl, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, SNTmpl->name, BPOut);
  SNRow = ObitTableSNRowUnref(SNRow);

  /* Close BP table */
  retCode = ObitTableClose ((ObitTable*)BPOut, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, BPOut->name, BPOut);
  BPRow = ObitTableBPRowUnref(BPRow);

  return BPOut;
} /* end DummyBPTable  */

/**
 * Copy the contents of an SN table to a BP table for a given channel
 * The solutions in the SN table are kept as corrections and in the 
 * BP table as solutions so the gains are inverted.
 * \param SNTab  InputSN table
 * \param BPTab  BP table to be updated
 * \param chan   0-rel channel number 
 * \param err    ObitErr stack for reporting problems.
 */
void SN2BPTable (ObitTableSN *SNTab, ObitTableBP *BPTab, olong chan,
		 ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableSNRow *SNRow=NULL;
  ObitTableBPRow *BPRow=NULL;
  olong irow, orow, nchan, nif, npol, count, i, indx;
  ofloat amp, fblank = ObitMagicF();
  gboolean found = FALSE;
  gchar *routine = "SN2BPTable";

  /* error checks */
  if (err->error) return;
  g_assert (ObitTableBPIsA(BPTab));
  if (SNTab==NULL) return;   /* Any data? */

  /* Table info - channels */
  nchan = BPTab->numChan;
  nif   = BPTab->numIF;
  npol  = BPTab->numPol;

  /* Open SN table */
  retCode = ObitTableSNOpen (SNTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, SNTab->name);
  SNRow = newObitTableSNRow(SNTab);

  /* Open BP table */
  retCode = ObitTableBPOpen (BPTab, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTab->name);
  BPRow = newObitTableBPRow(BPTab);
  ObitTableBPSetRow (BPTab, BPRow, err);
  if (err->error) Obit_traceback_msg (err, routine, BPTab->name);

  /* Make sure valid number on number of antennas */
  if ((BPTab->numAnt<=0) && (SNTab->numAnt>0))
    BPTab->numAnt    = SNTab->numAnt;  /* Max. antenna number */

  /* Loop over SN table */
  orow = 0;
  for (irow=1; irow<=SNTab->myDesc->nrow; irow++) {
    retCode = ObitTableSNReadRow (SNTab, irow, SNRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    SNRow->TimeI  = fabs(SNRow->TimeI); /* Grumble, grumble */

    /* Set time, antenna etc
       BPRow->Time.  = SNRow->Time;
       BPRow->TimeI. = SNRow->TimeI;
       BPRow->SourID.= SNRow->SourID;
       BPRow->antNo  = SNRow->antNo;
       BPRow->SubA.  = SNRow->;
       BPRow->FreqID.= SNRow->FreqID;.*/

    /* Find corresponding BP row */
    orow++;
    count = 0;
    retCode = ObitTableBPReadRow (BPTab, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
    /*found = (BPRow->Time==SNRow->Time) && (BPRow->antNo==SNRow->antNo);*/
    found = (fabs(BPRow->Time-SNRow->Time)<=MAX(BPRow->TimeI,SNRow->TimeI)) 
      && (BPRow->antNo==SNRow->antNo);
    while (!found) {
      /* Forward?? */
      if ((BPRow->Time<SNRow->Time) || 
	  ((fabs(BPRow->Time-SNRow->Time)<=MAX(BPRow->TimeI,SNRow->TimeI)) 
	   && (BPRow->antNo<SNRow->antNo))) {
	orow++;
	retCode = ObitTableBPReadRow (BPTab, orow, BPRow, err);
	if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
	found = (fabs(BPRow->Time-SNRow->Time)<=MAX(BPRow->TimeI,SNRow->TimeI)) 
	  && (BPRow->antNo==SNRow->antNo);
	/* Reverse? */
      } else if ((BPRow->Time>SNRow->Time) || 
		 ((fabs(BPRow->Time-SNRow->Time)<=MAX(BPRow->TimeI,SNRow->TimeI)) 
		 && (BPRow->antNo>SNRow->antNo))) {
	orow--;
	retCode = ObitTableBPReadRow (BPTab, orow, BPRow, err);
	if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
	found = (fabs(BPRow->Time-SNRow->Time)<MAX(BPRow->TimeI,SNRow->TimeI)) 
	  && (BPRow->antNo==SNRow->antNo);
      }
      found = (fabs(BPRow->Time-SNRow->Time)<=MAX(BPRow->TimeI,SNRow->TimeI)) 
	&& (BPRow->antNo==SNRow->antNo);
      /* Not forever */
      count++;
      if (count>(2*SNTab->myDesc->nrow)) goto cleanup;
    }

    /* Copy info to BP Row */
    BPRow->RefAnt1   = SNRow->RefAnt1[0];
    if (npol>1) BPRow->RefAnt2   = SNRow->RefAnt2[0];

    for (i=0; i<nif; i++) {
      /* Invert if valid */
      if ((SNRow->Real1[i]!=fblank) && (SNRow->Imag1[i]!=fblank)) {
	amp = SNRow->Real1[i]*SNRow->Real1[i] + SNRow->Imag1[i]*SNRow->Imag1[i];
	if (amp>1.0e-20) {
	  SNRow->Real1[i] /= amp;
	  SNRow->Imag1[i] /= amp;
	}
      } /* end invert */
      indx = i*nchan + chan;
      BPRow->Real1[indx]   = SNRow->Real1[i];
      BPRow->Imag1[indx]   = SNRow->Imag1[i];
      if (SNRow->Weight1[i]!=fblank) BPRow->Weight1[i]    = MAX (SNRow->Weight1[i], BPRow->Weight1[i]);
      if (npol>1) {
	/* Invert if valid */
	if ((SNRow->Real2[i]!=fblank) && (SNRow->Imag2[i]!=fblank)) {
	  amp = SNRow->Real2[i]*SNRow->Real2[i] + SNRow->Imag2[i]*SNRow->Imag2[i];
	  if (amp>1.0e-20) {
	    SNRow->Real2[i] /= amp;
	    SNRow->Imag2[i] /= amp;
	  }
	} /* end invert */
	BPRow->Real2[indx]   = SNRow->Real2[i];
	BPRow->Imag2[indx]   = SNRow->Imag2[i];
	if (SNRow->Weight2[i]!=fblank) BPRow->Weight2[i]    = MAX (SNRow->Weight2[i], BPRow->Weight2[i]);
      }
    }

    /* reWrite output table */
    retCode = ObitTableBPWriteRow (BPTab, orow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  } /* end loop over SN table */

  /* Close SN table */
 cleanup:
  retCode = ObitTableSNClose (SNTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, SNTab->name);
  SNRow = ObitTableSNRowUnref(SNRow);

  /* Close BP table */
  retCode = ObitTableBPClose (BPTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTab->name);
  BPRow = ObitTableBPRowUnref(BPRow);

} /* end SN2BPTable  */

/*----------------------------------------------------------------------- */
/*  Derive bandpass from autocorrelations                                 */
/*  Each spectrum is normalized to an average value of 1.0                */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to derive Bandpass for                           */
/*      outData   ObitUV to Write BP table on                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void AutoCorrBP (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		 ObitErr* err)
{
  ObitIOCode retCode;
  ObitTableBP* BPTable;
  ObitTableBPRow* BPRow=NULL;
  ObitUVDesc *inDesc;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, j, lrec, ver, suba;
  olong lastSubA=0, lastSourceID=0, lastFQID=0, curSourceID=0;
  olong  numAnt, numChan, numIF, numStoke, indx, jndx;
  olong corrType=2, refAnt = 1, iant, iif, istok, ichan, incs, incf, incif;
  olong ant1, ant2, ivis, orow, count = 0, cnt1, cnt2;
  ofloat cbase, curTime, norm1=1.0, norm2=1.0;
  ofloat **BPSum=NULL, **BPWt=NULL, solInt, *inBuffer;
  ofloat fblank = ObitMagicF();
  odouble startTime=0.0, endTime=0.0, lastTime = -1.0e20; 
  gboolean OK;
  gchar *routine = "AutoCorrBP";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));

  /* Only want autocorrelations */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inData->info, "corrType", OBIT_int, dim, &corrType);
  
  /* Solution interval */
  solInt = 0.0;
  ObitInfoListGetTest(myInput, "solInt2", &type, dim, &solInt);
  solInt /= 1440.0;     /* to days */
  if (solInt<=1.0e-5) solInt = 1000.0; /* 0 => whole data set */

  /* Reference antenna */
  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);
  ver = 0;
  ObitInfoListGetTest(myInput, "BPSoln", &type, dim, &ver);

  /* Open */
  retCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Get descriptor */
  inDesc  = inData->myDesc;

  /* Save first source id in case it's the only one */
  if (inData->mySel->sources) curSourceID = inData->mySel->sources[0];

  /* Number of channels, antennas */
  numChan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>0) numIF   = inDesc->inaxes[inDesc->jlocif];
  else numIF = 1;
  numStoke = inDesc->inaxes[inDesc->jlocs];
  numStoke = MIN (2, numStoke);
  suba = 1;
  numAnt  = inData->myDesc->numAnt[suba-1];/* actually highest antenna number */
  incs  = inDesc->incs;
  incf  = inDesc->incf;
  incif = inDesc->incif;
  lrec  = inDesc->lrec;

  /* Create output BP table */
  BPTable = newObitTableBPValue ("Output BP", (ObitData*)outData, &ver,
			       OBIT_IO_WriteOnly, numStoke, numIF, numChan, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Tell which BP version */ 
  Obit_log_error(err, OBIT_InfoErr, "Writing autocorrelation BP Table %d", ver);

  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)BPTable, err);
  if (err->error) goto cleanup;

  /* Open BP table */
  retCode = ObitTableBPOpen (BPTable, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTable->name);
  BPRow = newObitTableBPRow(BPTable);
  ObitTableBPSetRow (BPTable, BPRow, err);
  if (err->error) Obit_traceback_msg (err, routine, BPTable->name);

  /* Set header values */
  BPTable->numAnt    = numAnt;  /* Max. antenna number */
  BPTable->numShifts = 0;
  BPTable->numChan   = numChan;
  BPTable->startChan = 1;
  BPTable->lowShift  = 1;
  BPTable->shiftInc  = 1;
  strncpy (BPTable->BPType, "          ", MAXKEYCHARTABLEBP);
  BPTable->myDesc->sort[0] = BPTable->TimeCol+1;  /* Sort order */
  BPTable->myDesc->sort[1] = BPTable->antNoCol+1;

  /* Initialize BP Row */
  BPRow->BW        = inDesc->cdelt[inDesc->jlocf];
  BPRow->ChanShift[0] = 0.0;
  BPRow->ChanShift[1] = 0.0;
  BPRow->RefAnt1   = refAnt;
  BPRow->RefAnt2   = refAnt;
  for (i=0; i<numIF; i++) BPRow->ChanShift[i] = 0.0;
  for (i=0; i<numIF; i++) BPRow->Weight1[i] = 0.0;
  if (numStoke>1) for (i=0; i<numIF; i++) BPRow->Weight2[i] = 0.0;
  for (i=0; i<numChan*numIF; i++) { 
    BPRow->Real1[i]   = fblank;
    BPRow->Imag1[i]   = fblank;
    if (numStoke>1) {
      BPRow->Real2[i]   = fblank;
      BPRow->Imag2[i]   = fblank;
    }
  }

  /* Create accumulators, is order: channel, Stokes, IF */
  BPSum = g_malloc0(numAnt*sizeof(ofloat*));
  BPWt = g_malloc0(numAnt*sizeof(olong*));
  for (i=0; i<numAnt; i++) {
    BPSum[i] = g_malloc0(numIF*numStoke*numChan*sizeof(ofloat));
    BPWt[i]  = g_malloc0(numIF*numStoke*numChan*sizeof(olong));
  }

  /* Loop over file */
  startTime = -9999.9;
  while (retCode==OBIT_IO_OK) {
      retCode = ObitUVReadSelect (inData, inData->buffer, err);
      if (retCode > OBIT_IO_EOF) goto cleanup;

      inBuffer = inData->buffer;  /* Buffer */
      /* Hack to write last accumulation */
      if (retCode==OBIT_IO_EOF) inDesc->numVisBuff = 1;

      /* loop over visibilities */
      for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
	curTime = inBuffer[inDesc->iloct]; /* Time */
	if (inDesc->ilocsu>=0) curSourceID = (olong)(inBuffer[inDesc->ilocsu]+0.5);
	cbase = inBuffer[inData->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	if (ant1!=ant2) {inBuffer += lrec; continue;}
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
	if (startTime < -1000.0) {  /* Set time window etc. if needed */
	  startTime = curTime;
	  lastTime  = curTime;
	  endTime   = startTime + solInt;
	  lastSourceID = curSourceID;
	}
	
	/* Still in current interval/source? */
	if ((curTime>=endTime) || (curSourceID != lastSourceID) ||
	    (inDesc->firstVis>inDesc->nvis) || (retCode==OBIT_IO_EOF)) {
	  /* No - write BP table - setup row */
	  BPRow->SourID = lastSourceID;
	  BPRow->SubA   = lastSubA;
	  BPRow->FreqID = lastFQID;
	  BPRow->Time   = 0.5*(lastTime+startTime);
	  BPRow->TimeI  = lastTime-startTime;
	  for (iant=0; iant<numAnt; iant++) {
	    BPRow->antNo   = iant+1;
	    OK = FALSE;
	    for (iif=0; iif<numIF; iif++) {
	      BPRow->Weight1[iif] = 0.0;
	      /* Normalize */
	      norm1 = 0.0; cnt1 = 0;
	      for (ichan=0; ichan<numChan; ichan++) { 
		jndx = iif*numChan*numStoke + ichan;
		if (BPWt[iant][jndx]>0.0) {
		  norm1 += (BPSum[iant][jndx]/BPWt[iant][jndx]);
		  cnt1++;
		}
	      }
	      if (cnt1>=1) norm1 = cnt1/norm1;
	      else norm1 = 1.0;
	      if (numStoke>1) {
		BPRow->Weight2[iif] = 0.0;
		norm2 = 0.0; cnt2 = 0;
		for (ichan=0; ichan<numChan; ichan++) { 
		  jndx = iif*numChan*numStoke + numChan + ichan;
		  if (BPWt[iant][jndx]>0.0) {
		    norm2 += (BPSum[iant][jndx]/BPWt[iant][jndx]);
		    cnt2++;
		  }
		}
		if (cnt2>=1) norm2 = cnt2/norm2;
		else norm2 = 1.0;
	      } /* end setup for 2nd Poln */
	      for (ichan=0; ichan<numChan; ichan++) { 
		indx = iif*numChan + ichan;
		jndx = iif*numChan*numStoke + ichan;
		if (BPWt[iant][jndx]>0.0) {
		  BPRow->Real1[indx]  = sqrt(norm1*BPSum[iant][jndx]/BPWt[iant][jndx]);
		  BPRow->Imag1[indx]  = 0.0;
		  BPRow->Weight1[iif] = MAX (norm1*BPRow->Weight1[indx], BPWt[iant][jndx]);
		  OK = TRUE;
		} else {
		  BPRow->Real1[indx]   = fblank;
		  BPRow->Imag1[indx]   = fblank;
		}
	      } /* end loop over channel */
	      if (numStoke>1) {
		for (ichan=0; ichan<numChan; ichan++) { 
		  indx = iif*numChan + ichan;
		  jndx = iif*numChan*numStoke + numChan + ichan;
		  if (BPWt[iant][jndx]>0.0) {
		    BPRow->Real2[indx]  = sqrt(norm2*BPSum[iant][jndx]/BPWt[iant][jndx]);
		    BPRow->Imag2[indx]  = 0.0;
		    BPRow->Weight2[iif] = MAX (norm2*BPRow->Weight2[indx], BPWt[iant][jndx]);
		    OK = TRUE;
		  } else {
		    BPRow->Real2[indx]   = fblank;
		    BPRow->Imag2[indx]   = fblank;
		  }
		} /* end loop over channel */
	      } /* end if two stokes */
	      if (BPRow->Weight1[iif]>0.0) 
		BPRow->Weight1[iif] = MAX (1.0, BPRow->Weight1[iif]);
	      if ((numStoke>1) && (BPRow->Weight2[iif]>0.0) )
		BPRow->Weight2[iif] = MAX (1.0, BPRow->Weight2[iif]);
	    } /* end loop over IF */
	    
	    /* Write table */
	    if (OK) {
	      orow = -1;
	      count++;  /* How much good data? */
	      ObitTableBPWriteRow (BPTable, orow, BPRow, err);
	      if (err->error) goto cleanup;
	    }
	  } /* end loop over antenna */

	  /* reset accumulators */
	  for (i=0; i<numAnt; i++) {
	    for (j=0; j<numChan*numIF*numStoke; j++) {BPSum[i][j]=0.0; BPWt[i][j]=0.0;}
	  }
	  startTime = curTime;
	  lastTime  = curTime;
	  endTime   = startTime + solInt;
	  lastSourceID = curSourceID;
	} /* end if end of accumulation */

        /* Done? */
	if (retCode==OBIT_IO_EOF) break;
	  
	/* Accumulate spectra */
	/* loop over IF */
	for (iif=0; iif<numIF; iif++) {
	  /* Loop over polarization */
	  for (istok=0; istok<numStoke; istok++) {
	    /* Loop over frequency channel */
	    for (ichan=0; ichan<numChan; ichan++) { 
	      	indx = inDesc->nrparm + istok*incs + iif*incif + ichan*incf;
		jndx = iif*numChan*numStoke + istok*numChan + ichan;
		if (inBuffer[indx+2]>0.0) {
		  BPSum[ant1-1][jndx] += inBuffer[indx];
		  BPWt[ant1-1][jndx]  += inBuffer[indx+2];
		}
	    } /* end loop over Channel */
	  } /* end loop over Stokes */
	} /* end loop over IF */
	/* Save descriptive info*/
	if (inDesc->ilocfq>=0) lastFQID = (olong)(inBuffer[inDesc->ilocfq]+0.5);
	else lastFQID = 0;
	lastTime = curTime;
	lastSourceID = curSourceID;
	inBuffer += lrec;
      } /* end loop over buffer */
      /* Done? */
      if (retCode==OBIT_IO_EOF) break;
  } /* end loop over file */

  /* Close UV data */
  retCode = ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Close BP table */
  retCode = ObitTableBPClose (BPTable, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, BPTable->name);

  /* Cleanup */
 cleanup:
  BPTable = ObitTableBPUnref(BPTable );
  BPRow   = ObitTableBPRowUnref(BPRow);
  /* Delete accumulators */
  if (BPSum) {
    for (i=0; i<numAnt; i++) {
      if (BPSum[i]) g_free(BPSum[i]);
    }
    g_free(BPSum);
  }
  if (BPWt) {
    for (i=0; i<numAnt; i++) {
      if (BPWt[i]) g_free(BPWt[i]);
    }
    g_free(BPWt);
  }

  /* something actually done? */
   Obit_return_if_fail((count>0), err,
		      "%s: No valid autocorrelation data in %s", 
		      routine, inData->name);

}  /* end AutoCorrBP */
