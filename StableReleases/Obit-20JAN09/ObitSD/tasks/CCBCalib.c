/* $Id$  */
/* Obit Radio Single Disk OTF calibration software                    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitOTF.h"
#include "ObitOTFUtil.h"
#include "ObitOTFCal.h"
#include "ObitOTFCalUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitOTFGetAtmCor.h"
#include "ObitOTFSoln2Cal.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* CCBCalibIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void CCBCalibOut (ObitInfoList* outList, ObitErr *err);
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
/* Write history */
void CCBCalibHistory (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);
/* Clear any previous calibration */
void CCBClearCal (ObitOTF* inData, ObitErr* err);
/* Initial Cal table */
void CCBInitCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);
/* Nodding gain cal */
void CCBNodGainCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);
/* Atmospheric/gain/pointing calibration  */
olong CCBAtmCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);
/* Apply Soln to Cal */
void CCBUpdateCal (ObitInfoList* myInput, ObitOTF* inData, olong solnVer, 
		   olong CalIn, olong CalOut, ObitErr* err);
/* Baseline calibration  */
olong CCBBaselineCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);

/* Program globals */
gchar *pgmName = "CCBCalib";       /* Program name */
gchar *infile  = "CCBCalib.in" ;   /* File with program inputs */
gchar *outfile = "CCBCalib.out";   /* File to contain program outputs */
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
/*   Obit Radio Single Disk CCB calibration software                       */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  olong         snver;
  ObitSystem   *mySystem= NULL;
  ObitOTF       *inData = NULL;
  ObitErr      *err= NULL;

  /* Startup - parse command line */
  err = newObitErr();
  myInput = CCBCalibIn (argc, argv, err);
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

  /* Get input OTF data */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Clear previous calibration tables */
  CCBClearCal (inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create initial Cal table */
  CCBInitCal (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Nodding gain cal */
  CCBNodGainCal (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Atmospheric/gain/pointing calibration */
  snver = CCBAtmCal (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Apply Soln 1 to Cal 1 => Cal 2 */
  CCBUpdateCal (myInput, inData, snver, 1, 2, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Baseline fitting */
  snver = CCBBaselineCal(myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Apply Soln 2 to Cal 2 => Cal 3 */
  CCBUpdateCal (myInput, inData, snver, 2, 3, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  CCBCalibHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput); 
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* CCBCalibIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "CCBCalibIn";

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
} /* end CCBCalibIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: CCBCalib -input file -output ofile [args]\n");
    fprintf(stderr, "CCBCalib Obit task determine gains for OTF data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def CCBCalib.in\n");
    fprintf(stderr, "  -output uv data onto which to attach FG table, def CCBCalib.out\n");
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
  strTemp = "CCBCalib.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "CCBCalibName";
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
  btemp = FALSE;
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
    gboolean     doCalSelect;
    oint         doCCBCalib;*/
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));


} /* end digestInputs */

/*----------------------------------------------------------------------- */
/*  Get input data                                                        */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitOTF with input data                                           */
/*----------------------------------------------------------------------- */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitOTF       *inData = NULL;
  ObitInfoType type;
  olong         disk, nvis=1000;
  gchar        *strTemp, inFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Targets", "Scans", "timeRange", "flagVer",
    NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  inData = newObitOTF("input OTF data");
  

  /* Only FITS supported */
  if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (inFile, strTemp, 128);
  } else { 
    strncpy (inFile, "No_Filename_Given", 128);
  }
  
  /* input FITS disk */
  ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
  
  /* define object */
  ObitOTFSetFITS (inData, nvis, disk, inFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated and OK */
  ObitOTFFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Check that this is CCB data */
  Obit_retval_if_fail((inData->myDesc->OTFType==OBIT_GBTOTF_CCB), err, inData,
		      "%s: This is NOT CCB data", routine);
  
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Write History for CCBCalib                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to write history to                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CCBCalibHistory (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", 
    "Targets", "Scans", "timeRange", "calInter", 
    "NodScan", "NodFlux", "NodIndex", "tRx", "calJy",
    "tau0", "RAOff", "DecOff", "minEl", "order", "solInt", 
    NULL};
  gchar *routine = "CCBCalibHistory";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitOTFIsA(inData));
  
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
  
} /* end CCBCalibHistory  */

/*----------------------------------------------------------------------- */
/*  Delete and Cal or Soln Tables                                         */
/*   Input:                                                               */
/*      inData    ObitOTF to clear                                        */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CCBClearCal (ObitOTF* inData, ObitErr* err)
{
  olong ver;
  gchar *CalType="OTFCal", *SolnType="OTFSoln";
  gchar *routine = "CCBClearCal";
  
  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Removing previous calibration tables");
  ObitErrLog(err);  

  /* Delete Cal tables */
  ver = -1;
  ObitOTFZapTable (inData, CalType, ver, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Delete Soln tables */
  ver = -1;
  ObitOTFZapTable (inData, SolnType, ver, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
} /* end  CCBClearCal */

/*----------------------------------------------------------------------- */
/*  Create Cal Table 1   with interval calInter                           */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to add Cal table to                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CCBInitCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  olong ver;
  oint ncoef;
  ofloat calInt;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableOTFCal *CalTable=NULL;
  gchar *routine = "CCBInitCal";
  
  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Create initial Cal table");
  ObitErrLog(err);  

  /* Convert calInter to Days */
  calInt = 1.0;  /*Default 1 min. */
  ObitInfoListGetTest(myInput, "calInter", &type, dim, &calInt);
  if (calInt<=0) calInt = 1.0;  
  calInt /= 1440.0;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inData->info, "solInt", OBIT_float, dim, &calInt);
  
  /* Create Table */
  ver = 1;
  ncoef = 1;
  CalTable = ObitOTFGetDummyCal (inData, inData, ver, ncoef, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Don't need table */
  CalTable = ObitTableOTFCalUnref(CalTable);
} /* end  CCBInitCal */

/*----------------------------------------------------------------------- */
/*  Get calibration info from Nodding scan on calibrator                  */
/*  leaves results on inData info                                         */
/*  "tRx"      OBIT_float (*,1,1) Receiver temperature per detector in    */
/*                                units of the cal                        */
/*  "calJy"    OBIT_float (*,1,1) Noise cal value in Jy, per detector     */
/*  If no calibrator nodding scan is specified, these values are taken    */
/*  from myInput, if a nodding scan is used, these values are copied      */
/*  to myInput to be saved in the history                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to add Cal table to                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CCBNodGainCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  olong calScan, scans[2], ndetect;
  ofloat calFlux, calIndex, xtemp[32];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar   *dataParms[] = {  /* Calibration parameters from myInput */
    "calJy", "tRx", 
    NULL};
  gchar *routine = "CCBNodGainCal";
  
  /* Get calibration info using nodding */
  ObitInfoListGet(myInput, "NodScan", &type, dim, &calScan, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  if (calScan<=0) {
    ndetect = MAX (1, inData->geom->numberDetect);
    if (ObitInfoListGetTest(myInput, "calJy", &type, dim, xtemp)) {
      dim[0] = MIN (dim[0],ndetect);  /* No more than actual number of detectors */
      ObitInfoListAlwaysPut (inData->info, "calJy", type, dim, xtemp);}

    if (ObitInfoListGetTest(myInput, "tRx", &type, dim, xtemp)) {
      dim[0] = MIN (dim[0],ndetect);
      ObitInfoListAlwaysPut (inData->info, "tRx", type, dim, xtemp);}
    
    Obit_log_error(err, OBIT_InfoErr, 
		   "Using input parameters in lieu of a nodding scan");
    ObitErrLog(err);  
    return;
  }

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Fitting Nodding Scan");
  ObitErrLog(err);  

  /* Calibrator flux density/spectral index */
  calFlux = 1.0;  /* Default 1 Jy */
  ObitInfoListGetTest(myInput, "NodFlux", &type, dim, &calFlux);
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inData->info, "calFlux", OBIT_float, dim, &calFlux);
  calIndex = 1.0;  /* Default -0.7 */
  ObitInfoListGetTest(myInput, "NodIndex", &type, dim, &calIndex);
  ObitInfoListAlwaysPut(inData->info, "calIndex", OBIT_float, dim, &calIndex);
  
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inData->info, "Scan", OBIT_long, dim, &calScan);
  
  /* Do fitting */
  ObitOTFCalUtilFitNod (inData, -1, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Reset scan selection */
  ObitInfoListGet(myInput, "Scans", &type, dim, scans, err);
  if (scans[0]<=0) scans[0] = -1000000;
  if (scans[1]<=0) scans[1] =  1000000;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inData->info, "Scans", OBIT_long, dim, scans);
  /* Save calibration parameters on myInput */
  ObitInfoListCopyList (inData->info, myInput, dataParms);
} /* end  CCBNodGainCal */

/*----------------------------------------------------------------------- */
/*  Apply Soln table to a Cal Table writing a new one                     */
/*   Input:                                                               */
/*      inData    ObitOTF to update Cal table                             */
/*      solnVer   Input Soln table version                                */
/*      CalIn     Input Cal table version                                 */
/*      CalOut    Output Cal table version                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void CCBUpdateCal (ObitInfoList* myInput, ObitOTF* inData, olong solnVer, 
		   olong CalIn, olong CalOut, ObitErr* err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitTableOTFCal *CalTable=NULL;
  gchar *routine = "CCBUpdateCal";

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Apply Soln %d to Cal %d, write Cal %d",
		 solnVer, CalIn, CalOut);
  ObitErrLog(err);  

  /* Set tables */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inData->info, "solnUse", OBIT_long, dim, &solnVer);
  ObitInfoListAlwaysPut(inData->info, "calIn", OBIT_long, dim, &CalIn);
  ObitInfoListAlwaysPut(inData->info, "calOut", OBIT_long, dim, &CalOut);

  CalTable = ObitOTFSoln2Cal (inData, inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Don't need table */
  CalTable = ObitTableOTFCalUnref(CalTable);
 
} /* end  CCBUpdateCal */

/*----------------------------------------------------------------------- */
/*  Atmospheric/gain/pointing calibration                                 */
/*  Writes Soln Table 1                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to add Cal table to                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
olong CCBAtmCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  olong snver=0;
  olong i, ndetect;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean keepCal;
  ObitTableOTFSoln *solnTable;
  ofloat farr[1000], solint, tau0;
  gchar   *dataParms[] = {  /* Control Parameters */
    "RAoff", "Decoff", "minEl", 
    NULL};
  gchar *routine = "CCBAtmCal";

  /* Error check */  
  if (err->error) return snver;

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Atmospheric/gain/pointing calibration");
  ObitErrLog(err);  

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);

  /* Some need changing */
  /* Solution interval */
  solint = 1.0;  /* Default 1 min. */
  ObitInfoListGetTest(myInput, "solInt", &type, dim, &solint);
  if (solint<=0) solint = 1.0;  
  solint /= 1440.0;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inData->info, "solInt",  OBIT_float, dim, &solint);

  /* Tau0 */
  ObitInfoListGetTest(myInput, "tau0", &type, dim, &tau0);
  dim[0] = 1;
  ObitInfoListAlwaysPut(inData->info, "Tau0",  OBIT_float, dim, &tau0);

  /* Arbitrary value for air temperature/cal */
  ndetect = MAX (1, inData->geom->numberDetect);
  for (i=0; i<ndetect; i++) farr[i] = 1.0;
  dim[0] = ndetect;
  ObitInfoListAlwaysPut(inData->info, "aTemp",   OBIT_float, dim, farr);

  /* calType "GainOffset" */
  dim[0] = 10;
  ObitInfoListAlwaysPut(inData->info, "calType", OBIT_string, dim, "GainOffset");

  /* Need cal-on data */
  keepCal = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListPut (inData->info, "keepCal", OBIT_bool, dim, &keepCal, err);

  /* Determine atmospheric calibration */
  solnTable = ObitOTFGetAtmCor (inData, inData, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, snver);

  /* Get version number */
  snver = solnTable->tabVer;

  /* Don't need table */
  solnTable = ObitTableOTFSolnUnref(solnTable);
  return snver;
} /* end  CCBAtmCal */

/*----------------------------------------------------------------------- */
/*  Baseline calibration                                                  */
/*  Writes Soln Table 2                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to add Cal table to                             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return version of Soln table                                         */
/*----------------------------------------------------------------------- */
olong CCBBaselineCal (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  olong snver=0;
  olong gainuse, itemp;
  ofloat solint;
  gboolean doCalSelect, keepCal;
  ObitTableOTFSoln *solnTable;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar   *dataParms[] = {  /* Control Parameters */
    "ORDER", "MINEL", 
    NULL};
  gchar *routine = "CCBNodGainCal";
  
  /* Error check */  
  if (err->error) return snver;

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Baseline calibration");
  ObitErrLog(err);  

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, inData->name, snver);

  /* Solution interval - convert to min. */
  solint = 1.0;  /* Default 1 min. */
  ObitInfoListGetTest(myInput, "solInt", &type, dim, &solint);
  if (solint<=0) solint = 1.0;  
  solint /= 1440.0;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inData->info, "solInt",  OBIT_float, dim, &solint);

  /* Apply previous calibration  */
  dim[0] = 1; dim[1] = 1;
  gainuse = 0;
  ObitInfoListPut (inData->info, "gainUse", OBIT_long, dim, &gainuse, err);
  itemp = 1;
  ObitInfoListAlwaysPut (inData->info, "doCalib", OBIT_long, dim, &itemp);
  doCalSelect = TRUE;
  ObitInfoListAlwaysPut (inData->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Drop cal-on data */
  keepCal = FALSE;
  ObitInfoListAlwaysPut (inData->info, " keepCal", OBIT_bool, dim, &keepCal);
  
  /* Fit polynomial to baseline */
  solnTable = ObitOTFGetSolnPolyBL (inData, inData, err);

  /* Get version number */
  snver = solnTable->tabVer;

  /* Don't need table */
  solnTable = ObitTableOTFSolnUnref(solnTable);
  return snver;
} /* end CCBBaselineCal */

