/* $Id$  */
/* Obit task to Flag selected UV-data       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
#include "ObitUVUtil.h"
#include "ObitUVEdit.h"
#include "ObitTableFG.h"
#include "ObitHistory.h"
#include "ObitData.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVFlagIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void UVFlagOut (ObitInfoList* outList, ObitErr *err);
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
/* Flag all selected data */
void  FlagData(ObitUV* inData, ObitErr* err);
/* Write history */
void UVFlagHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Program globals */
gchar *pgmName = "UVFlag";       /* Program name */
gchar *infile  = "UVFlag.in" ;   /* File with program inputs */
gchar *outfile = "UVFlag.out";   /* File to contain program outputs */
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
/*   Obit task to automatically edit a uv data set                        */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitTableFG  *FlagTab=NULL;
  olong        flagTab=0;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitErr      *err= NULL;
  gchar        opCode[5];

   /* Startup - parse command line */
  err = newObitErr();
  myInput = UVFlagIn (argc, argv, err);
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

  /* Make sure flagTab exists, assign if passed 0 */
  ObitInfoListGetTest(myInput, "flagTab",  &type, dim, &flagTab);
  FlagTab = newObitTableFGValue("tmpFG", (ObitData*)inData, &flagTab, OBIT_IO_ReadWrite, 
				err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  /* Open/close table */
  ObitTableFGOpen (FlagTab, OBIT_IO_ReadWrite, err);
  ObitTableFGClose (FlagTab, err);
  FlagTab = ObitTableFGUnref(FlagTab);
  /* Save for history */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(myInput, "flagTab",  OBIT_long, dim, &flagTab);

  /* opCode */
  strncpy (opCode, "FLAG", 4);
  ObitInfoListGetTest(myInput, "opCode",  &type, dim, opCode); opCode[4] = 0;

  /* Do editing */
  if (!strncmp(opCode, "FLAG", 4)) {                                       /* Straight flagging */
    FlagData(inData, err);
  } else if (!strncmp(opCode, "ELEV", 4)) {                                /* Elevation flagging */
    ObitUVEditElev (inData, inData, err);
  } else if (!strncmp(opCode, "SHAD", 4) || !strncmp(opCode, "CROS", 4) || /* Shadow/cross talk */
	    !strncmp(opCode, "SHCR", 4)) {
    ObitUVEditShadCross (inData, inData, err);  
  } else {                                                                 /* Unknown */
    Obit_log_error(err, OBIT_Error, "%s: Unknown opCode %s", pgmName, opCode);
  }
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Write history */
  UVFlagHistory (myInput, inData, err); 
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  FlagTab   = ObitUnref(FlagTab);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);
  myOutput  = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem  = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* UVFlagIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "UVFlagIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

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
} /* end UVFlagIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVFlag -input file -output ofile [args]\n");
    fprintf(stderr, "UVFlag Obit task flag Selected UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVFlag.in\n");
    fprintf(stderr, "  -output parameter file, def UVFlag.out\n");
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
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
  gboolean btemp;
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultInputs";*/

  /* error checks */
  if (err->error) return out;

  /* add parser items */
  /* Program number */
  dim[0] = 1; dim[1] = 1;
  itemp = 1;
  ObitInfoListAlwaysPut (out, "pgmNumber", OBIT_oint, dim, &itemp);

  /* Default FITS directories - same directory */
  dim[0] = 1; dim[1] = 1;
  itemp = 0; /* number of FITS directories */
  ObitInfoListAlwaysPut (out, "nFITS", OBIT_oint, dim, &itemp);

  /* AIPS user number */
  dim[0] = 1; dim[1] = 1;
  itemp = 2;
  ObitInfoListAlwaysPut (out, "AIPSuser", OBIT_oint, dim, &itemp);

  /* Default AIPS directories */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; /* number of AIPS directories */
  ObitInfoListAlwaysPut (out, "nAIPS", OBIT_oint, dim, &itemp);

  /* Default type "FITS" */
  strTemp = "FITS";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "DataType", OBIT_string, dim, strTemp);

  /* input FITS file name */
  strTemp = "UVFlag.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inFile", OBIT_string, dim, strTemp);

  /* input AIPS file name */
  strTemp = "UVFlagName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inName", OBIT_string, dim, strTemp);

  /* input AIPS file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "inClass", OBIT_string, dim, strTemp);

  /* AIPS sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListAlwaysPut (out, "inSeq", OBIT_oint, dim, &itemp);

  /* AIPS or FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListAlwaysPut (out, "inDisk", OBIT_oint, dim, &itemp);

  /* Default opCode = "FLAG" */
  strTemp = "FLAG";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "opCode", OBIT_string, dim, strTemp);

  /* Stokes parameter to edit */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (out, "Stokes", OBIT_string, dim, strTemp);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListAlwaysPut (out, "timeRange", OBIT_float, dim, farray);

  /*  Apply calibration/selection?, def=True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListAlwaysPut (out, "doCalSelect", OBIT_bool, dim, &btemp);

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
  ofloat ftemp;
  gchar *strTemp, Stokes[5], opCode[5];
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Initialize Threading */
  ObitThreadInit (myInput);

  /* Change Stokes to PFlag */
  strncpy (Stokes, "    ", 4);
  ObitInfoListGetTest(myInput, "Stokes",  &type, dim, Stokes); Stokes[4] = 0;
  ObitInfoListAlwaysPut (myInput, "PFlag", OBIT_string, dim, Stokes);
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListAlwaysPut (myInput, "Stokes", OBIT_string, dim, strTemp);

  /* Check Shadow/Crosstalk limits */
  strncpy (opCode, "FLAG", 4);
  ObitInfoListGetTest(myInput, "opCode",  &type, dim, opCode); opCode[4] = 0;
  ftemp = 0.0;
  dim[0] = dim[1] = dim[2] = 1;
  if (!strncmp(opCode, "SHAD", 4)) {
    ObitInfoListAlwaysPut (myInput, "minCross", OBIT_float, dim, &ftemp);
  }
  if (!strncmp(opCode, "CROS", 4)) {
    ObitInfoListAlwaysPut (myInput, "minShad", OBIT_float, dim, &ftemp);
  }

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
  olong        nvis=1000, flagver;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual","PFlag", "timeRange", 
    "BChan", "EChan", "BIF", "EIF", 
    "subA", "Antennas", "Baseline", "flagTab", "FreqID", 
    "opCode", "minElev", "minShad", "minCross", "Reason",
    "doCalSelect", NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Set buffer size */
  nvis = 1000;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO",  OBIT_long, dim,  &nvis);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Ensure inData fully instantiated and OK and selector set 
     Don't want flagging */
  flagver = -1;
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(inData->info, "flagVer", OBIT_int, dim, &flagver);
  
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Flag all selected data                                                */
/*   Input:                                                               */
/*      inData    ObitUV to flag                                          */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void  FlagData(ObitUV* inData, ObitErr* err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        flagTab=0;
  ObitTableFG  *FlagTab=NULL;
  ObitTableFGRow *row=NULL;
  olong *Antennas=NULL, nant, *Baseline=NULL, nbl, nsou, *souId;
  olong i, iFGRow, isou, iant, ibl, count;
  gchar Stokes[5], Reason[25];
  gchar *routine = "FlagData";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  
  /* Get flag table */
  ObitInfoListGetTest(inData->info, "flagTab",  &type, dim, &flagTab);
  FlagTab = newObitTableFGValue("outFG", (ObitData*)inData, &flagTab, OBIT_IO_ReadWrite, 
				err);
  ObitTableFGOpen (FlagTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Create Row */
  row = newObitTableFGRow (FlagTab);
  
  /* Attach row to output buffer */
  ObitTableFGSetRow (FlagTab, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Get selected data */
  row->TimeRange[0] = 0.0; row->TimeRange[1] = 1.0e20;
  ObitInfoListGetTest(inData->info, "timeRange",  &type, dim, row->TimeRange); 
  ObitInfoListGetTest(inData->info, "subA",  &type, dim, &row->SubA);
  ObitInfoListGetTest(inData->info, "FreqID",  &type, dim, &row->freqID);
  row->ifs[0] = 1; row->ifs[1] = 0; 
  ObitInfoListGetTest(inData->info, "BIF",  &type, dim, &row->ifs[0]);
  ObitInfoListGetTest(inData->info, "EIF",  &type, dim, &row->ifs[1]);
  row->ifs[0] =  MAX (1, row->ifs[0]);
  row->chans[0] = 1; row->chans[1] = 0; 
  ObitInfoListGetTest(inData->info, "BChan",  &type, dim, &row->chans[0]);
  ObitInfoListGetTest(inData->info, "EChan",  &type, dim, &row->chans[1]);
  row->chans[0] =  MAX (1, row->chans[0]);
  strncpy (Reason, "                        ", 24);
  ObitInfoListGetTest(inData->info, "Reason",  &type, dim, Reason); Reason[24] = 0;
  strncpy (row->reason, Reason, 24);
  /* Stokes flags */
  strncpy (Stokes, "    ", 4);
  ObitInfoListGetTest(inData->info, "PFlag",  &type, dim, Stokes); Stokes[4] = 0;
  if      (!strncmp(Stokes, "    ", 4)) row->pFlags[0] = 1+2+4+8;     /* All */
  else if (!strncmp(Stokes, "I   ", 4)) row->pFlags[0] = 1+2;         /* I */
  else if (!strncmp(Stokes, "Q   ", 4)) row->pFlags[0] = 1+2+4+8;     /* Q */
  else if (!strncmp(Stokes, "U   ", 4)) row->pFlags[0] = 1+2+4+8;     /* U */
  else if (!strncmp(Stokes, "V   ", 4)) row->pFlags[0] = 1+2+4+8;     /* V */
  else if (!strncmp(Stokes, "RR  ", 4) || !strncmp(Stokes, "XX  ", 4)) row->pFlags[0] = 1;
  else if (!strncmp(Stokes, "LL  ", 4) || !strncmp(Stokes, "YY  ", 4)) row->pFlags[0] = 2;
  else if (!strncmp(Stokes, "RL  ", 4) || !strncmp(Stokes, "XY  ", 4)) row->pFlags[0] = 4;
  else if (!strncmp(Stokes, "LR  ", 4) || !strncmp(Stokes, "YY  ", 4)) row->pFlags[0] = 8;
  else if ((Stokes[0]=='1') || (Stokes[0]=='0')) {
    /* flag bits passed */
    row->pFlags[0] = 0;
    if (Stokes[0]=='1')  row->pFlags[0] |= 1;
    if (Stokes[1]=='1')  row->pFlags[0] |= 2;
    if (Stokes[2]=='1')  row->pFlags[0] |= 4;
    if (Stokes[3]=='1')  row->pFlags[0] |= 8;
  } else row->pFlags[0] = 1+2+4+8;   /* what the hell */
  
  /* Open and close data to get selection */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Source selection from UV data*/
  nsou  = inData->mySel->numberSourcesList;
  souId = inData->mySel->sources;
  if (souId) {
    count = 0;
    for (i=0; i<nsou; i++) {
      if (souId[i]<0) break;
      count++;
    }
    nsou = count;
    if (nsou<=0) souId = NULL;
  }

  /* Antenna selection */
  ObitInfoListGetP(inData->info, "Antennas",  &type, dim, (gpointer)&Antennas);
  nant = 0;
  if (Antennas!=NULL) { /* specified */
    for (i=0; i<dim[0]; i++) {
      if (Antennas[i]<=0) break;
      nant++;
    }
  } /* end Antennas specified */
  if (nant<=0) Antennas = NULL;  /* None given */

  /* Baseline selection */
  ObitInfoListGetP(inData->info, "Baseline",  &type, dim, (gpointer)&Baseline);
  nbl = 0;
  if (Baseline!=NULL) { /* specified */
    for (i=0; i<dim[0]; i++) {
      if (Baseline[i]<=0) break;
      nbl++;
    }
  } /* end Baselines specified */
  if (nbl<=0) Baseline = NULL;  /* None given */

  /* Loop over sources */
  for (isou=0; isou<MAX(1,nsou); isou++) {
    if (souId) row->SourID = souId[isou];
    else       row->SourID = 0;
    /* Loop over antennas */
    for (iant=0; iant<MAX(1,nant); iant++) {
      if (Antennas) row->ants[0] = Antennas[iant];
      else          row->ants[0] = 0;
      /* Loop over baselines */
      for (ibl=0; ibl<MAX(1,nbl); ibl++) {
	if (Baseline) row->ants[1] = Baseline[ibl];
	else          row->ants[1] = 0;
	/* Write flag */
	iFGRow = -1;
	ObitTableFGWriteRow (FlagTab, iFGRow, row, err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);
      } /* end loop over baseline */
    } /* end loop over antenna */
  } /* end loop over sources */

  /* Close output table */
  ObitTableFGClose (FlagTab, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Cleanup */
  FlagTab = ObitTableFGUnref(FlagTab);
  row     = ObitTableFGRowUnref(row);
} /* end FlagData */

/*----------------------------------------------------------------------- */
/*  Write History for UVFlag                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to update history                                */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVFlagHistory (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq", 
    "BChan", "EChan", "BIF", "EIF", "FreqID", "souCode", "Qual", 
    "Sources", "PFlag", "timeRange",  "subA", "Antennas", "Baseline", 
    "flagTab", "opCode", "minElev", "minShad", "minCross", "Reason",
    NULL};
  gchar *routine = "UVFlagHistory";

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
 
} /* end UVFlagHistory  */

