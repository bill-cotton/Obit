/* $Id: MazrCal.c 225 2010-09-24 00:53:00Z bill.cotton $  */
/* VLBI Maser calibration                                 */
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
#include "ObitUVGSolve.h"
#include "ObitUVGSolveWB.h"
#include "ObitSkyModel.h"
#include "ObitSkyModelMF.h"
#include "ObitTableSUUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitUVUtil.h"
#include "ObitSkyGeom.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* MazrCalIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void MazrCalOut (ObitInfoList* outList, ObitErr *err);
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
/* Divide by FT of channel models */
void DivideModel (ObitInfoList *myInput, ObitSkyModel *skyModel, 
		  ObitUV *scrUV, ObitErr *err);
/* Write history */
void MazrCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Program globals */
gchar *pgmName = "MazrCal";       /* Program name */
gchar *infile  = "MazrCal.in" ;   /* File with program inputs */
gchar *outfile = "MazrCal.out";   /* File to contain program outputs */
olong  pgmNumber;                 /* Program number (like POPS no.) */
olong  AIPSuser;                  /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;                   /* Number of AIPS directories */
gchar **AIPSdirs=NULL;            /* List of AIPS data directories */
olong  nFITS=0;                   /* Number of FITS directories */
gchar **FITSdirs=NULL;            /* List of FITS data directories */
ObitInfoList *myInput  = NULL;    /* Input parameter list */
ObitInfoList *myOutput = NULL;    /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   VLBI Maser calibration                                               */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  const ObitUVGSolveClassInfo *solnClass;
  ObitUV       *inData = NULL, *scrData = NULL;
  ObitSkyModel *skyModel=NULL;
  ObitErr      *err= NULL;
  ObitUVGSolve *solver=NULL;
  ObitTableSN  *SNTable = NULL;
  ofloat       ftemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *solverParms[] = {  /* Calibration parameters */
    "solInt", "solnVer", "solType", "solMode", "avgPol", "avgIF", "doMGM", "elevMGM",
    "refAnt", "ampScalar", "minSNR",  "minNo", "prtLv",
    NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = MazrCalIn (argc, argv, err);
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
  if (err->prtLv>=3) {
    Obit_log_error(err, OBIT_InfoErr, "Calibrate/edit/select data");
    ObitErrLog(err);
  }
  scrData = newObitUVScratch (inData,err);
  scrData = ObitUVCopy (inData, scrData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Divide by FT of model, channel at a time */
  DivideModel (myInput, skyModel, scrData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Index scrData */
  dim[0] = dim[1] = 1;
  ftemp = 15.0;  /* Max scan time 15 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxScan", OBIT_float, dim, &ftemp);
  ftemp = 1.0; /* Max allowable gap 1 min. */
  ObitInfoListAlwaysPut(scrData->info, "maxGap", OBIT_float, dim, &ftemp);
  ObitUVUtilIndex (scrData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create Calibration solution object */
  solver = ObitUVGSolveCreate("Gain solver");

  /* Save first source id in case it's the only one */
  if (inData->mySel->sources) solver->curSource = inData->mySel->sources[0];

  /* Copy calibration control to solver */
  ObitInfoListCopyList (myInput, solver->info, solverParms);

  /* Do solution */
  solnClass = (ObitUVGSolveClassInfo*)solver->ClassInfo;
  SNTable   = solnClass->ObitUVGSolveCal (solver, scrData, inData, inData->mySel, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Write history */
  MazrCalHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  solver    = ObitUVGSolveUnref(solver);
  SNTable   = ObitTableSNUnref(SNTable);
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  scrData   = ObitUnref(scrData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* MazrCalIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "MazrCalIn";

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
} /* end MazrCalIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: MazrCal -input file -output ofile [args]\n");
    fprintf(stderr, "MazrCal VLBI Maser calibration \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def MazrCal.in\n");
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
  strTemp = "MazrCal.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "MazrCalName";
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

  /* default Spectral index */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Alpha", OBIT_float, dim, &ftemp, err);
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
  gchar        *strTemp;
  ObitSkyModelMode modelMode;
  ObitSkyModelType modelType;
  ofloat       modelFlux, refAnt, *refAnts;
  gboolean     doCalSelect;
  oint         doCalib;
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
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Copy first refAnts value to refAnt */
   if (ObitInfoListGetP(myInput, "refAnts",  &type, dim, (gpointer)&refAnts)) {
     refAnt = refAnts[0];
   } else {
     refAnt = 0;
     /* For old inputs */
     ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);
   }
   dim[0] = dim[1] = dim[2] = 1;
   ObitInfoListAlwaysPut (myInput, "refAnt", OBIT_long, dim, &refAnt);

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
    "Sources", "souCode", "Qual", "subA", "selBand", "selFreq", "FreqID",
    "Stokes", "timeRange", "BChan", "EChan",   "BIF", "EIF", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol",
    "Antennas", "Mode", "ModelType", "Alpha",
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
    nThreads = 1;
    ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
    nvis *= nThreads;
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
  ObitInfoType type;
  gboolean      do3D=TRUE;
  olong         Aseq, disk, cno,i=0, ver, nmaps;
  gchar        *Type, *Type2, *strTemp, inFile[129], inRoot[129];
  gchar        Aname[13], Aclass[7], Aroot[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar        name[101];
  gchar        *dataParms[] = {  /* Control parameters */
    "BComp",  "EComp",  "Flux", "PBCor", "antSize", "Factor", 
    "minFlux", "Mode", "ModelType", "Stokes", "prtLv", "noNeg",
    NULL};
  gchar *routine = "getInputSkyModel";

  /* error checks */
  if (err->error) return skyModel;
  g_assert (ObitInfoListIsA(myInput));

  /* How many fields? */
  nmaps = 1;
  ObitInfoListGetTest(myInput, "nfield", &type, dim, &nmaps);
  nmaps = MAX (1,nmaps);   /* Must have at least one */

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
      
      /* Ensure image fully instantiated and OK */
      ObitImageFullInstantiate (image[i], TRUE, err);
      
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
      
      /* Ensure image fully instantiated and OK */
      ObitImageFullInstantiate (image[0], TRUE, err);
      
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
	
	/* Ensure image fully instantiated and OK */
	ObitImageFullInstantiate (image[i], TRUE, err);
	
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
 * \param myInput   Input parameters on InfoList          
 * \param skyModel  Basic sky model.
 * \param scrUV     Uv data object to be corrected.
 * \param err       ObitErr stack for reporting problems.
 */
void DivideModel (ObitInfoList *myInput, ObitSkyModel *skyModel, 
		  ObitUV *scrUV, ObitErr *err)
{
  olong ichan, nchan, hiver, itemp, prtLv=0, CCVer1=1;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat xShift, yShift, shift;
  gchar *routine = "DivideModel";

  /* error checks */
  if (err->error) return;
  g_assert (ObitSkyModelIsA(skyModel));
  g_assert (ObitUVIsA(scrUV));

  /* Make sure there are enough CC tables for all channels */
  ObitInfoListGetTest (myInput, "prtLv", &type, dim, &prtLv); 
  ObitInfoListGetTest (myInput, "CCVer1", &type, dim, &CCVer1); 
  hiver = ObitTableListGetHigh(skyModel->mosaic->images[0]->tableList, "AIPS CC");
  nchan = scrUV->myDesc->inaxes[scrUV->myDesc->jlocf];
  Obit_return_if_fail ((hiver>=(CCVer1+nchan-1)), err, 
			"%s not enough CC Tables %d < %d", 
			routine, hiver, CCVer1+nchan-1);

  /* Make sure the shift is less than 0.1 degree */
  ObitSkyGeomShiftXY (scrUV->myDesc->crval[scrUV->myDesc->jlocr], 
		      scrUV->myDesc->crval[scrUV->myDesc->jlocd], 
		      scrUV->myDesc->crota[scrUV->myDesc->jlocd],
		      skyModel->mosaic->images[0]->myDesc->crval[0],
		      skyModel->mosaic->images[0]->myDesc->crval[1],
		      &xShift, &yShift);
  shift = sqrt (xShift*xShift + yShift*yShift);
  Obit_return_if_fail ((shift<0.1), err, 
			"%s Excessive shift %f > %f", 
			routine, shift, 0.1);

  /* scrUV better be asingle source file */
  Obit_return_if_fail ((scrUV->myDesc->ilocsu<=0), err, 
		       "You must select a single source");
  
  /* Loop over channels dividing */
  for (ichan=0; ichan<nchan; ichan++) {
    /* Set info on skymodel */
    dim[0] = dim[1] = dim[2] = 1;
    itemp = CCVer1+ichan;
    ObitInfoListAlwaysPut (skyModel->info, "CCVer", OBIT_long, dim, &itemp);
    if (prtLv>=3) {
      Obit_log_error(err, OBIT_InfoErr, "Dividing channel %d", itemp);
      ObitErrLog(err);
    }
    itemp = ichan+1;
    ObitInfoListAlwaysPut (skyModel->info, "BChan", OBIT_long, dim, &itemp);
    ObitInfoListAlwaysPut (skyModel->info, "EChan", OBIT_long, dim, &itemp);
    ObitSkyModelDivUV (skyModel, scrUV, scrUV, err);
    if (err->error) Obit_traceback_msg (err, routine, scrUV->name);
  } /* end loop over channels*/

} /* end DivideModel  */

/*----------------------------------------------------------------------- */
/*  Write History for MazrCal                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void MazrCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA",
    "selBand", "selFreq", "FreqID", "BChan", "EChan", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ", "BPVer", "flagVer", 
    "doPol", "Antennas",  "refAnt", 
    "DataType2", "in2File", "in2Disk", "in2Name", "in2Class", "in2Seq", 
    "nfield", "CCVer1", "BComp", "EComp", "Cmethod", "Cmodel", "Flux", "Alpha",
    "solInt", "solType", "solMode", "avgPol", "avgIF", "noNeg", "doMGM", "minSNR",
    "minNo", "nThreads",
   NULL};
  gchar *routine = "MazrCalHistory";

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
 
} /* end MazrCalHistory  */

