/* $Id$  */
/* R-L phase bandpass calibration                                     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2014                                          */
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
ObitInfoList* RLPassIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void RLPassOut (ObitInfoList* outList, ObitErr *err);
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
void RLBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		   ObitErr* err);
/* Divide by SU table flux density */
void DivideSource (ObitUV *inUV, ObitUV *scrUV, ObitErr *err);
/* Divide buffer table flux density */
static void DivideBuffer (ObitSourceList *sList, olong sourId, ObitUV *uvdata);
/* Write history */
void RLPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);
/* Initial calibration */
ObitUV* InitialCal (ObitInfoList* myInput, ObitUV* scrData, ObitErr* err);
/* Create  BP table */
ObitTableBP* MakeBPTable (ObitUV* inData, ofloat time,
			  olong bchan, olong echan, olong nif, 
			  ofloat *phase, ofloat RLPhase, ofloat RM, olong refAnt,
			  ObitErr *err);
/* Update  BP table */
ObitTableBP* UpdateBPTable (ObitUV* inData, olong BPver, ObitUV* outData, 
			    ofloat time, olong bchan, olong echan, olong nif, 
			    ofloat *phase, ofloat RLPhase, ofloat RM, olong refAnt,
			    ObitErr *err);


/* Program globals */
gchar *pgmName = "RLPass";      /* Program name */
gchar *infile  = "RLPass.in" ;  /* File with program inputs */
gchar *outfile = "RLPass.out";  /* File to contain program outputs */
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
/*    R-L phase bandpass calibration                                      */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL, *scrData = NULL, *avgData = NULL;;
  ObitSkyModel *skyModel=NULL;
  ObitErr      *err= NULL;
  ofloat       ftemp;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *SCParms[] = {  /* Selfcal Parameters  */
    "refAnt",
     NULL};

  /* Startup - parse command line */
  err = newObitErr();
  myInput = RLPassIn (argc, argv, err);
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
  RLBandpassCal(myInput, avgData, inData, err);
  if (err->error) ierr = 1;   ObitErrLog(err);  if (ierr!=0) goto exit;

  /* Write history */
  RLPassHistory (myInput, inData, err); 
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

ObitInfoList* RLPassIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "RLPassIn";

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
} /* end RLPassIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: RLPass -input file -output ofile [args]\n");
    fprintf(stderr, "RLPass Obit task determine bandpass for UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def RLPass.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def RLPass.out\n");
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
  strTemp = "RLPass.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "RLPassName";
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
  gboolean     doCalSelect;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "FreqID", "BChan", "EChan", "BIF", "EIF", 
    "subA", "Antennas", "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", 
    "flagVer", "doPol", "PDVer", "Mode", "ModelType", "Alpha", "refAnt",
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
  if (err->error) Obit_traceback_msg (err, routine, scrUV->name);
  
  if (SUTable) ObitTableSUUnref(SUTable);
  if (sList) ObitSourceListUnref(sList);
  if (err->error) Obit_traceback_msg (err, routine, scrUV->name);

} /* end DivideSource  */

/**
 * Divide buffer load of data by reciprocal of flux density in SourceList
 * \param sList   Source list with 1/flux densities
 * \param sourId  Source ID (0-rel) in sList is uvdata a single source file
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
/*  Write History for RLPass                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void RLPassHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "Sources", "Qual", "souCode", "timeRange",  "subA",
    "selBand", "selFreq", "FreqID", 
    "BChan1", "EChan1", "BChan2", "EChan2", "ChWid2", 
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol", "PDVer", "Antennas",  "refAnt",  
    "DataType2", "in2File", "in2Disk", "in2Name", "in2Class", "in2Seq", 
    "nfield", "CCVer", "BComp", "EComp", "Cmethod", "Cmodel", "Flux",
    "modelFlux", "modelPos", "modelParm", "RLPhase", "RM",
    "solInt1", "solInt2",  "minSNR",  "minNo", "BPSoln", "prtLv", "nThreads",
   NULL};
  gchar *routine = "RLPassHistory";

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
 
} /* end RLPassHistory  */

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
void  RLBandpassCal(ObitInfoList* myInput, ObitUV* avgData, ObitUV* inData, 
		    ObitErr* err)
{
  ObitTableBP *BPTable=NULL, *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  olong ichan, nchan, bchan2=0, echan2=0, chinc2=0, nc, i;
  olong ivis, nvis, ifreq, nfreq, iif, nstok, nif, indx, count=0;
  olong itemp, nchanIF, doBand, BPVer, outVer, npol, BPSoln, refAnt;
  gboolean btemp;
  ofloat RLPhase=0.0, RM=0.0, fblank = ObitMagicF();
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *vis, *ifvis, *fvis, *svis, *rec, time=-1.0e20;
  ObitUVDesc *desc;
  gfloat *Qsum=NULL, *Usum=NULL, *Qwt=NULL, *Uwt=NULL;
  gchar *Stokes="IQUV";
  gchar *routine = "RLBandpassCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(avgData));


  /* Any data in avgData? */
  if (avgData->myDesc->nvis<=0) {
    Obit_log_error(err, OBIT_Error, "NO averaged/calibrated data");
    return;
  }

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
  /* Stokes */
  if (avgData->myDesc->jlocs>=0)
    npol = MIN (2, avgData->myDesc->inaxes[avgData->myDesc->jlocs]);
  else
    npol = 1; 

  /* Target */
  ObitInfoListGetTest (myInput, "RM",      &type, dim, &RM); 
  ObitInfoListGetTest (myInput, "RLPhase", &type, dim, &RLPhase); 
  RLPhase *= DG2RAD;    /* to radians */

  /* Prior bandpass calibration */
  ObitInfoListGetTest (myInput, "doBand", &type, dim, &doBand); 
  ObitInfoListGetTest (myInput, "BPVer",  &type, dim, &BPVer); 
  ObitInfoListGetTest (myInput, "BPSoln", &type, dim, &BPSoln);
  if (BPVer<=0) BPVer = ObitTableListGetHigh (inData->tableList, "AIPS BP");
  ObitInfoListGetTest (myInput, "refAnt", &type, dim, &refAnt); 

  /* Number of IFs */
  if (avgData->myDesc->jlocif>=0)  
    nif = avgData->myDesc->inaxes[avgData->myDesc->jlocif];
  else
    nif = 1;

  /* Accumulator arrays */
  nc      =  (echan2-bchan2+1);
  nchanIF =  nif * nc;
  Qsum = g_malloc0(nchanIF*sizeof(ofloat));
  Qwt  = g_malloc0(nchanIF*sizeof(ofloat));
  Usum = g_malloc0(nchanIF*sizeof(ofloat));
  Uwt  = g_malloc0(nchanIF*sizeof(ofloat));

  /* Don't need to apply calibration */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(avgData->info, "doCalSelect", OBIT_bool, dim, &btemp);
  itemp = -1;
  ObitInfoListAlwaysPut(avgData->info, "doCalib", OBIT_long, dim, &itemp);
  /* Want to convert to IQUV */
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut(avgData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Loop over channels getting channel calibration */
  for (ichan=bchan2; ichan<=echan2; ichan++) {
    /* Select channel(s) */
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    itemp = MAX (1, ichan-chinc2);
    ObitInfoListAlwaysPut(avgData->info, "BChan", OBIT_long, dim, &itemp);
    itemp = MIN (nchan, ichan+chinc2);
    ObitInfoListAlwaysPut(avgData->info, "EChan", OBIT_long, dim, &itemp);
    
    /* Sum one block of channels and all IFs */
    retCode = ObitUVOpen (avgData, OBIT_IO_ReadCal, err);
    if (err->error) goto cleanup;
    desc  = avgData->myDesc;

    nfreq = desc->inaxes[desc->jlocf];
    nif = 1;
    if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
    nstok = 1;
    if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];

    /* loop over blocks of data */
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
      if (time<-1.0e10) time = rec[desc->iloct];
      for (ivis=0; ivis<nvis; ivis++) {
	
	/* loop over IFs */
	ifvis = vis;
	for (iif=0; iif<nif; iif++) {
	  indx = iif*nc + ichan-bchan2;
	  
	  /* loop over frequencies being averaged */
	  fvis = ifvis;
	  for (ifreq = 0; ifreq<nfreq; ifreq++) {
	    
	    /*  Stokes */
	    svis = fvis;
	    /* I - not needed */
	    svis += desc->incs; /* visibility pointer */
	    /* Q */
	    if (*(svis+2)>0.0) {
	      Qsum[indx] += *svis * *(svis+2);
	      Qwt[indx]  += *(svis+2);
	    }
	    svis += desc->incs; 
	    /* U */
	    if (*(svis+2)>0.0) {
	      Usum[indx] += *svis * *(svis+2);
	      Uwt[indx]  += *(svis+2);
	    }

	    fvis += desc->incf; /* visibility pointer */
	  } /* end loop over frequencies */
	  ifvis += desc->incif; /* visibility pointer */
	} /* Loop over IFs */
	
	/* update data pointers */
	vis += desc->lrec;
	rec += desc->lrec;
      } /* end loop over buffer */
  
    } /* end loop over data */

    /* Close data */
    retCode = ObitUVClose (avgData, err);
    if (err->error) goto cleanup;
    }  /* end channel loop */

  /* Normalize - convert to angle */
  for (i=0; i<nchanIF; i++) {
    if (Qwt[i]>0.0) Qsum[i] /= Qwt[i];
    if (Uwt[i]>0.0) Usum[i] /= Uwt[i];
    if ((Qwt[i]>0.0) && (Uwt[i]>0.0)) {
      Qsum[i] = atan2(Usum[i], Qsum[i]);
      count++;
    }
    else Qsum[i] = fblank;
  }
		      
  /* tell about it */
  if (count>0) {
    Obit_log_error(err, OBIT_InfoErr, "Calibrated %d channels", count);
  } else {
    Obit_log_error(err, OBIT_Error, "NO channels calibrated");
    goto cleanup;
 }

  /* Create new or update prior depending on doBand */
  if (doBand<=0) {
    /* Create/populate BP table */
    BPTable = MakeBPTable (avgData, time, bchan2, echan2, nif, 
			   Qsum,  RLPhase, RM, refAnt, err);
    outVer = BPSoln;
    if (outVer<=0) 
      outVer = ObitTableListGetHigh (inData->tableList, "AIPS BP")+1;
    Obit_log_error(err, OBIT_InfoErr, 
		   "New BP table ver. %d", outVer);
  } else {
    BPTable = UpdateBPTable (inData, BPVer, avgData, time, bchan2, echan2, nif, 
			     Qsum,  RLPhase, RM, refAnt, err);
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
  if (Qsum) g_free(Qsum);
  if (Usum) g_free(Usum);
  if (Qwt)  g_free(Qwt);
  if (Uwt)  g_free(Uwt);
  ObitUVZapTable (avgData, "AIPS BP", -1, err);
  BPTable = ObitTableBPUnref(BPTable );
  ObitUVZapTable (avgData, "AIPS SN", -1, err);
  if (err->error) Obit_traceback_msg (err, routine, avgData->name);
  } /* end RLBandpassCal  */

/** 
 * Create Bandpass table, creates entries at time and for antennas numbers
 * up to the max in the data.
 * \param inData  UV onto which the BP table is to be attached
 * \param time    Time in days to write entry
 * \param bchan   First channel calibrated
 * \param echan   Highest channel calibrated
 * \param nif     Number of IFs
 * \param phase   Channel phase in radians
 * \param RLPhase Target R-L phase (rad) at reference freq
 * \param RM      Rotation measure (rad/M^2)
 * \param refAnt  Reference antenna
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* MakeBPTable (ObitUV* inData, ofloat time,
			  olong bchan, olong echan, olong nif, 
			  ofloat *phase, ofloat RLPhase, ofloat RM,
			  olong refAnt, ObitErr *err)
{
  ObitTableBP *BPOut=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  ObitUVDesc *desc;
  olong i, irow, orow, nchan, npol, ver, numAnt, ichan, iif, nc;
  olong indx, jndx, suba=1;
  odouble freq, reffreq, lambda, reflambda;
  ofloat pdif, fblank = ObitMagicF();
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
  BPOut->numAnt    = numAnt;  /* Max. antenna number */
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

  /* Reference wavelength */
  reffreq   = desc->freq;
  reflambda = VELIGHT/reffreq;

  /* Fill actual values */
  nc      =  (echan-bchan+1);
  for (iif=0; iif<nif; iif++) {
    for (ichan=bchan; ichan<=echan; ichan++) {
      /* channel wavelength */
      freq   = reffreq + desc->freqIF[iif] + 
	(ichan-desc->crpix[desc->jlocf]+1) * desc->cdelt[desc->jlocf];
      lambda = VELIGHT/freq;
      /* Target phase difference */
      pdif = RLPhase + 2.0 * (lambda*lambda - reflambda*reflambda) * RM;
      indx = iif*nc + ichan-bchan;
      jndx = iif*nchan + ichan -1;
      if (phase[indx]!=fblank) {
	BPRow->Real1[jndx]   = 1.0;
	BPRow->Imag1[jndx]   = 0.0;
	if (npol>1) {   /* There better be or this is pretty pointless */
	  BPRow->Real2[jndx] = cos(pdif-phase[indx]);
	  BPRow->Imag2[jndx] = sin(pdif-phase[indx]);
	}
      } else {
	BPRow->Real1[jndx]   = fblank;
	BPRow->Imag1[jndx]   = fblank;
	if (npol>1) { 
	  BPRow->Real2[jndx] = fblank;
	  BPRow->Imag2[jndx] = fblank;
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
 * \param RLPhase Target R-L phase (rad) at reference freq
 * \param RM      Rotation measure (rad/M^2)
 * \param refAnt  Reference antenna
 * \param err     ObitErr stack for reporting problems.
 */
ObitTableBP* UpdateBPTable (ObitUV* inData, olong BPver, ObitUV* outData, 
			    ofloat time, olong bchan, olong echan, olong nif, 
			    ofloat *phase, ofloat RLPhase, ofloat RM,
			    olong refAnt, ObitErr *err)
{
  ObitTableBP *BPOut=NULL, *BPIn=NULL;
  ObitIOCode retCode = OBIT_IO_OK;
  ObitTableBPRow *BPRow=NULL;
  ObitUVDesc *desc;
  olong irow, orow, nchan, npol, ver, ichan, iif, nc;
  olong indx, jndx;
  odouble freq, reffreq, lambda, reflambda;
  ofloat pdif, xr, xi, yr, yi, fblank = ObitMagicF();
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

  /* Loop over input table updating */
  for (irow=1; irow<=BPIn->myDesc->nrow; irow++) {
    retCode = ObitTableBPReadRow (BPIn, irow, BPRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Reference wavelength */
    reffreq   = desc->freq;
    reflambda = VELIGHT/reffreq;
    
    /* Fill actual values */
    BPRow->RefAnt1 = refAnt;
    BPRow->RefAnt2 = refAnt;
    nc      =  (echan-bchan+1);
    for (iif=0; iif<nif; iif++) {
      for (ichan=bchan; ichan<=echan; ichan++) {
	/* channel wavelength */
	freq   = reffreq + desc->freqIF[iif] + 
	  (ichan-desc->crpix[desc->jlocf]+1) * desc->cdelt[desc->jlocf];
	lambda = VELIGHT/freq;
	/* Target phase difference */
	pdif = RLPhase + 2.0 * (lambda*lambda - reflambda*reflambda) * RM;
	indx = iif*nc + ichan-bchan;
	jndx = iif*nchan + ichan -1;
	if (phase[indx]!=fblank) {
	  BPRow->Real1[jndx]   = 1.0;
	  BPRow->Imag1[jndx]   = 0.0;
	  if (npol>1) {   /* There better be or this is pretty pointless */
	    yr = BPRow->Real2[jndx];
	    yi = BPRow->Imag2[jndx];
	    xr = cos(pdif-phase[indx]);
	    xi = sin(pdif-phase[indx]);
	    BPRow->Real2[jndx] = xr*yr - xi*yi;
	    BPRow->Imag2[jndx] = xr*yi + xi*yr;
	  }
	} else {  /* No value - blank */
	  BPRow->Real1[jndx]   = fblank;
	  BPRow->Imag1[jndx]   = fblank;
	  if (npol>1) { 
	    BPRow->Real2[jndx] = fblank;
	    BPRow->Imag2[jndx] = fblank;
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

