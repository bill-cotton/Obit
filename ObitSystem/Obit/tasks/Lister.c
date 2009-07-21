/* $Id:  $  */
/* Task to print the contents of various data files                   */
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

#include <unistd.h>
#include "ObitPrinter.h"
#include "ObitUV.h"
#include "ObitUVSel.h"
#include "ObitUVUtil.h"
#include "ObitTable.h"
#include "ObitTableUtil.h"
#include "ObitTableCL.h"
#include "ObitTableSN.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableNX.h"
#include "ObitTableFQ.h"
#include "ObitPosLabelUtil.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* ListerIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void ListerOut (ObitInfoList* outList, ObitErr *err);
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
/* Basic UV data display */
void doDATA (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Scan listing */
void doSCAN (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Gain listing */
void doGAIN (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Days to human string */
void day2dhms(ofloat time, gchar *timeString);
/* Extract data value from SN row */
ofloat getSNValue (ObitTableSNRow *SNRow, gboolean firstPol, 
		   olong iif, olong dtype);
/* Extract data value from CL row */
ofloat getCLValue (ObitTableCLRow *CLRow, gboolean firstPol, 
		   olong iif, olong dtype);
/* Set gain scaling for SN table */
ofloat getSNGainScale (ObitUVSel *sel, ObitTableSN *SNTable, ObitTableSNRow *SNRow, 
		       olong maxAnt, gboolean *doAnt, olong loAnt, olong hiAnt, 
		       olong iif, gboolean firstPol, olong dt,  olong ndig, 
		       ObitErr *err);
/* Set gain scaling for CL table */
ofloat getCLGainScale (ObitUVSel *sel, ObitTableCL *CLTable, ObitTableCLRow *CLRow, 
		       olong maxAnt, gboolean *doAnt, olong loAnt, olong hiAnt, 
		       olong iif, gboolean firstPol, olong dt, olong ndig, 
		       ObitErr *err);

/* Program globals */
gchar *pgmName = "Lister";       /* Program name */
gchar *infile  = "Lister.inp";   /* File with program inputs */
gchar *outfile = "Lister.out";   /* File to contain program outputs */
olong  pgmNumber;                /* Program number (like POPS no.) */
olong  AIPSuser;                 /* AIPS user number number  */
olong  nAIPS=0;                  /* Number of AIPS directories */
gchar **AIPSdirs=NULL;           /* List of AIPS data directories */
olong  nFITS=0;                  /* Number of FITS directories */
gchar **FITSdirs=NULL;           /* List of FITS data directories */
ObitInfoList *myInput  = NULL;   /* Input parameter list */
ObitInfoList *myOutput = NULL;   /* Output parameter list */


/* File info */
gchar  *Type;        /* File type, AIPS or FITS */
gchar  inFile[513];  /* FITS input file */
gchar  Aname[13];    /* AIPS Name */
gchar  Aclass[7];    /* AIPS class */
olong  Aseq;         /* AIPS sequence */
olong  disk;         /* AIPS or FITS disk */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Task to print the contents of various data files                     */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitErr      *err= NULL;
  ObitUV       *inData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        opType[24];

   /* Startup - parse command line */
  err = newObitErr();
  myInput = ListerIn (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get control parameters from myInput */
  ObitInfoListGet(myInput, "opType", &type, dim, opType, err);
  opType[dim[0]] = 0;  /* Be sure NULL terminated */

  /* operation by opType */
  if (!strncmp(opType, "DATA", 4)) {
    /* basic uv data display */
    doDATA (myInput, inData, err);
  } else if (!strncmp(opType, "SCAN", 4)) {
    /* Scan listing */
    doSCAN (myInput, inData, err);
  } else if (!strncmp(opType, "GAIN", 4)) {
    /* Gain listing */
    doGAIN (myInput, inData, err);
  } else {
    /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown opType %s", 
                   pgmName, opType);
  }

  /* error check */
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput = ObitInfoListUnref(myInput);    /* delete input list */
  inData  = ObitUnref(inData);
  
  /* Shutdown  */
 exit:
  /* Python kinda slow witted, wait for I/0 */
  usleep(500000); /* 0.5 sec */
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);                /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* ListerIn (int argc, char **argv, ObitErr *err)
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
  ObitInfoList* list;
  gchar *routine = "ListerIn";

  /* Make default inputs, outputs InfoList */
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
} /* end ListerIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Lister -input file -output ofile [args]\n");
    fprintf(stderr, "List UV related data in various forms\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Lister.in\n");
    fprintf(stderr, "  -output output result file, def Lister.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
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
/*     inDisk    Int        input AIPS or FITS image disk no  [def 1]     */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  ObitInfoList *out = newObitInfoList();
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
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
  strTemp = "GetJy.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "GetJyName";
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

 return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     Mean     Int        Image pixel mean  [0.0]                        */
/*     RMS      Int        Image pixel rms   [0.0]                        */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  ObitInfoList *out = newObitInfoList();
  /*  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
      ofloat ftemp;
      gchar *routine = "defaultOutputs";*/

  /* add parser items */
  
  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Digest inputs                                                         */
/*  Rewrites DataType as InDataType                                       */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void digestInputs(ObitInfoList *myInput, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  doCalSelect = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Rewrite DataType as "inDataType" */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  ObitInfoListAlwaysPut (myInput, "inDataType", type, dim, Type);

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
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *strTemp;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Qual", "calCode", "Stokes", "timeRange",  "FreqID", 
    "BIF", "EIF", "BChan", "EChan", "Antennas", "subA",
    "doCalib", "gainUse", "doPol", "flagVer", "doBand", "BPVer", "Smooth", 
    "doCalSelect",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Build basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Fetch input file definition */
  /* Disk number */
  disk = 0;
  ObitInfoListGetTest(myInput, "inDisk", &type, dim, &disk);
  /* AIPS Sequence */
  Aseq = 0;
  ObitInfoListGetTest(myInput, "inSeq",  &type, dim, &Aseq);
  /* AIPS Name */
  if (ObitInfoListGetP(myInput, "inName", &type, dim, (gpointer)&strTemp)) {
    strncpy (Aname, strTemp, 13);
  } else { /* Didn't find */
    strncpy (Aname, "No Name ", 13);
  } 
  /* AIPS Class */
  if (ObitInfoListGetP(myInput, "inClass", &type, dim, (gpointer)&strTemp)) {
    strncpy (Aclass, strTemp, 7);
  } else { /* Didn't find */
    strncpy (Aclass, "NoClas", 7);
  } 
  /* FITS name */
  if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (inFile, strTemp, 128);
  } else { 
    strncpy (inFile, "No_Filename_Given", 128);
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
/*  Basic UV data display                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input UV data                                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doDATA (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitPrinter  *myPrint = NULL;
  ObitIOCode   iretCode;
  ObitUVDesc   *inDesc;
  FILE         *outStream = NULL;
  gboolean     isInteractive = FALSE, quit = FALSE, first=TRUE, doReal=FALSE;
  gchar        line[1024], Title1[1024], Title2[1024];
  olong        LinesPerPage = 0;
  olong        i, ii, indx, jndx, count, bPrint, nPrint, doCrt=1, lenLine=0;
  ofloat       u, v, w, cbase, re, im, amp, phas, wt;
  ofloat       maxBL=0.0, maxWt=0.0, minWt=1.0e20, maxAmp=0.0;
  ofloat       blscale=1., wtscale=1., ampscale=1.;
  olong        start, ic, ncor, mcor, maxcor, ant1, ant2, souID, SubA, inc=2;
  olong        bif=1, bchan=1, ichan, doCalib, gainUse, flagVer;
  gchar        *prtFile=NULL, timeString[25];
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *cstokes[4] = {"RR", "LL", "RL", "LR"};
  gchar        *stokes[4]  = {"I ", "Q ", "U ", "V "};
  gchar        *routine = "doDATA";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get parameters */
  ObitInfoListGet(myInput, "bPrint",    &type, dim, &bPrint, err);
  ObitInfoListGet(myInput, "nPrint",    &type, dim, &nPrint, err);
  ObitInfoListGet(myInput, "doCalib",   &type, dim, &doCalib, err);
  ObitInfoListGet(myInput, "gainUse",   &type, dim, &gainUse, err);
  ObitInfoListGet(myInput, "flagVer",   &type, dim, &flagVer, err);
  ObitInfoListGetTest(myInput, "inc",   &type, dim, &inc);
  inc = MAX (1, inc);
  nPrint *= inc;   /* Correct for increment */
  ObitInfoListGetTest(myInput, "BIF",   &type, dim, &bif);
  ObitInfoListGetTest(myInput, "BChan", &type, dim, &bchan);
  ObitInfoListGet(myInput, "doReal",    &type, dim, &doReal, err);
  ObitInfoListGet(myInput, "doCrt",     &type, dim, &doCrt,  err);
  isInteractive = doCrt>0;
  if (isInteractive) { /* interactive, write to stdout */
    outStream = stdout;
  } else { /* write to file */
    ObitInfoListGetP(myInput, "prtFile",  &type, dim, (gpointer)&prtFile);
    if (prtFile) prtFile[MIN (47,dim[0])] = 0;
    ObitTrimTrail(prtFile);
    /* Make sure file named */
    Obit_return_if_fail(((strlen(prtFile)>2) && 
			 ((prtFile[0]!=' ') && (prtFile[1]!=' '))), err,
			"%s: Printer file not specified", routine);
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Writing output to %s",prtFile);
    ObitErrLog(err);
 }
  if (err->error) Obit_traceback_msg (err, routine, "myInput");

  /* How long is a line? */
  lenLine = MAX (72, MIN(doCrt,1023));
  if (!isInteractive) lenLine = 132;

  /* Open uv data */
  iretCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inData->name);
  count = 0;
  inDesc = inData->myDesc;  

   /* Titles */
  sprintf (Title1, "                                                   ");
  sprintf (Title2, "        Time         basel souID  U     V      W   ");

  /* How many vis to print in Stokes */
  start = strlen(Title2);
  maxcor = (lenLine-start) / 18;
  if (maxcor==3) maxcor = 2;  /* Only even number or 1 */
  ncor = MIN (inDesc->inaxes[inDesc->jlocs], maxcor);  /* Stokes */

  /* If more room available, loop in freq */
  if ((maxcor-ncor)>=2) {
    mcor = MIN (inDesc->inaxes[inDesc->jlocf], (maxcor/ncor));   /* Number freq */
  } else mcor = 1;

  /* Add visibility labels */
  for (ichan=0; ichan<mcor; ichan++) {  /* Freq loop */
    for (ic=0; ic<ncor; ic++) {  /* Stokes loop */
      start = strlen(Title2);
      if (doReal) 
	sprintf (&Title2[start], "   Real   Imag  Wt");
      else /* Amp/phase */
	sprintf (&Title2[start], "    Amp  Phase  Wt");
      start = strlen(Title1);
      ii = (olong)(fabs(inDesc->crval[inDesc->jlocs]) - 0.9);
      if (inDesc->crval[inDesc->jlocs]>0.0)
	sprintf (&Title1[start], "    %2d/%4d %s    ", 
		 bif,bchan+ichan,stokes[ic+ii]);
      else
	sprintf (&Title1[start], "    %2d/%4d %s    ", 
		 bif,bchan+ichan,cstokes[ic+ii]);
    }
  } /* end frequency loop */

 /* Create/Open Printer */
  myPrint = ObitPrinterCreate (pgmName, isInteractive, outStream, prtFile);
  ObitPrinterOpen (myPrint, LinesPerPage, Title1, Title2, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Write header */
  if (!strncmp(Type,"AIPS",4)) {
    sprintf (line, "AIPS %s %s %d disk %d user %d",Aname,Aclass,Aseq, disk, AIPSuser);
  } else if (!strncmp(Type,"FITS",4)) {
    sprintf (line, "FITS %s disk %d",inFile, disk);
  } else {
    sprintf (line, "UNKNOWN Data type ");
  }
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  sprintf(line,"Object: %s Telescope: %s Observed: %s  Freq: %8.3lf GHz", 
	  inDesc->object,inDesc->teles, inDesc->obsdat, inDesc->freq*1.0e-9);
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* calibration */
  if (doCalib>0) {
    sprintf(line,"   Calibrating data with gain table %d", gainUse);
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  }

  /* Editing */
  if (flagVer>0) {
    sprintf(line,"   Editing data with flag table %d", flagVer);
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  }

  /* Loop through data */
  while ((iretCode==OBIT_IO_OK) && (count<nPrint)) {
    /* read buffer full */
    iretCode = ObitUVReadSelect (inData, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* On first pass get scaling */
    if (first && (inDesc->numVisBuff>1)) {
      first = FALSE;
      for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
	indx = i*inDesc->lrec;
	jndx = i*inDesc->lrec + inDesc->nrparm;
	u   = inData->buffer[indx+inDesc->ilocu];
	maxBL = MAX (maxBL, fabs(u));
	v   = inData->buffer[indx+inDesc->ilocv];
	maxBL = MAX (maxBL, fabs(v));
	w   = inData->buffer[indx+inDesc->ilocw];
	maxBL = MAX (maxBL, fabs(w));
	re = inData->buffer[jndx];
	im = inData->buffer[jndx+1];
	amp = re*re+im*im;
	maxAmp = MAX(maxAmp, amp);
	wt  = inData->buffer[jndx+2];
	maxWt = MAX (wt, maxWt); minWt = MIN (wt, minWt);
      }
      /* Scalings */
      blscale = 1;
      if (maxBL>1.0e3) blscale = 1.0e-1;
      if (maxBL>1.0e4) blscale = 1.0e-2;
      if (maxBL>1.0e5) blscale = 1.0e-3;
      if (maxBL>1.0e6) blscale = 1.0e-4;
      if (maxBL>1.0e7) blscale = 1.0e-5;
      if (maxBL>1.0e7) blscale = 1.0e-6;
      if (maxBL>1.0e9) blscale = 1.0e-7;
      if (maxBL>1.0e10) blscale = 1.0e-8;
      if (maxBL>1.0e11) blscale = 1.0e-9;
      wtscale = 1.0;
      if (maxWt>1.0e2)  wtscale = 1.0e1;
      if (maxWt>1.0e3)  wtscale = 1.0;
      if (maxWt>1.0e4)  wtscale = 1.0e-1;
      if (maxWt>1.0e5)  wtscale = 1.0e-2;
      if (maxWt>1.0e6)  wtscale = 1.0e-3;
      if (maxWt<1.0e2)  wtscale = 1.0e3;
      if (maxWt<1.0e-1) wtscale = 1.0e4;
      if (maxWt<1.0e-2) wtscale = 1.0e5;
      if (maxWt<1.0e-3) wtscale = 1.0e6;
      if (maxWt<1.0e-4) wtscale = 1.0e7;
      if (maxWt<1.0e-5) wtscale = 1.0e7;
      if (maxWt<1.0e-6) wtscale = 1.0e9;
      if (maxWt<1.0e-7) wtscale = 1.0e10;
      if (maxWt<1.0e-8) wtscale = 1.0;
      sprintf( line, "   u,v,w scaled by %8.2g  weights scaled by  %8.2g",
	       blscale, wtscale);
      ObitPrinterWrite (myPrint, line,   &quit, err);
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

      ampscale = 1.0;
      maxAmp = sqrt(maxAmp);
      if (maxAmp>=1.3e0) ampscale = 1.0e1;
      if (maxAmp>=1.3e1) ampscale = 1.0e0;
      if (maxAmp>=1.3e2) ampscale = 1.0e-1;
      if (maxAmp>=1.3e3) ampscale = 1.0e-2;
      if (maxAmp>=1.3e4) ampscale = 1.0e-3;
      if (maxAmp>=1.3e5) ampscale = 1.0e-2;
      if (maxAmp>=1.3e6) ampscale = 1.0e-5;
      if (maxAmp<1.3e0)  ampscale = 1.0e2;
      if (maxAmp<=1.0e-1) ampscale = 1.0e3;
      if (maxAmp<=1.0e-2) ampscale = 1.0e4;
      if (maxAmp<=1.0e-3) ampscale = 1.0e5;
      if (maxAmp<=1.0e-4) ampscale = 1.0e6;
      if (maxAmp<=1.0e-5) ampscale = 1.0e7;
      if (maxAmp<=1.0e-6) ampscale = 1.0e8;
      if (ampscale!=1.0) {
	sprintf( line, "   Amplitudes scaled by %8.2g",ampscale);
	ObitPrinterWrite (myPrint, line,   &quit, err);
      }

      /* Titles */
      ObitPrinterWrite (myPrint, Title1, &quit, err);
      ObitPrinterWrite (myPrint, Title2, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    }
    
    /* loop through buffer */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Want this one% */
      if (((count%inc)!=0) || ((count+1)<bPrint)) goto skip;

      indx = i*inDesc->lrec;
      jndx = i*inDesc->lrec + inDesc->nrparm;

      /* Random parameters */
      u   = blscale * inData->buffer[indx+inDesc->ilocu];
      v   = blscale * inData->buffer[indx+inDesc->ilocv];
      w   = blscale * inData->buffer[indx+inDesc->ilocw];
      day2dhms (inData->buffer[indx+inDesc->iloct], timeString);
      cbase = inData->buffer[indx+inDesc->ilocb]; /* Baseline */
      ant1 = (cbase / 256.0) + 0.001;
      ant2 = (cbase - ant1 * 256) + 0.001;
      SubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
      if (inDesc->ilocsu>=0) 
	souID = (olong)(inData->buffer[indx+inDesc->ilocsu]+0.5);
      else souID = 0;

      /* Format line with random parms */
      sprintf(line,"%6d %s %2d-%2d %3d %6.1f %6.1f %6.1f",
	      count+1,timeString,ant1,ant2,souID, u,v,w);
      /* Add visibilities */
      for (ichan=0; ichan<mcor; ichan++) {  /* Freq loop */
	for (ic=0; ic<ncor; ic++) {  /* Stokes loop */
	  /* Visibilities */
	  re = inData->buffer[jndx  + ic*inDesc->incs + ichan*inDesc->incf];
	  im = inData->buffer[jndx+1+ ic*inDesc->incs + ichan*inDesc->incf];
	  wt   = wtscale * inData->buffer[jndx+2+ic*inDesc->incs + ichan*inDesc->incf];
	  amp  = ampscale*sqrt(re*re + im*im);
	  phas = RAD2DG * atan2(im, re+1.0e-20);
	  start = strlen(line);
	  if (doReal) 
	    sprintf(&line[start]," %6.2f %6.2f %3.0f", re*ampscale,im*ampscale,wt);
	  else /* Amp/phase */
	    sprintf(&line[start]," %7.3f %4.0f %4.0f", amp,phas,wt);

	}
      }

      /* Print line */
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
      
    skip:
      count++;
      if (count>=nPrint) break;
    } /* end loop over visibilities */

  } /* end loop reading buffer */
    
  /* Close Printer, data */
 Quit:
  ObitPrinterClose (myPrint, err);
  /* Close UV data */
  ObitUVClose(inData, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

} /* end doDATA */

/*----------------------------------------------------------------------- */
/*  Scan listing for multisource data set                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input UV data                                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doSCAN (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitPrinter  *myPrint = NULL;
  ObitTableNX  *NXTable=NULL;
  ObitTableNXRow *NXRow=NULL;
  ObitSourceList *SouList=NULL;
  ObitTableSU  *SUTable=NULL;
  ObitTableFQ  *FQTable=NULL;
  ObitTableFQRow *FQRow=NULL;
  ObitIOCode   iretCode;
  ObitUVDesc   *inDesc;
  FILE         *outStream = NULL;
  oint         numIF;
  gboolean     isInteractive = FALSE, quit = FALSE;
  gchar        line[1024], Title1[1024], Title2[1024];
  olong        i, iif, ivr, ivd, lenLine, LinesPerPage = 0, maxSUID, *noVis=NULL;
  olong        ver, iRow, count, start, souID, lastSouID, doCrt=1, qual = 0;
  gchar        *prtFile=NULL, btimeString[25], etimeString[25];
  gchar        source[17]="source          ", calcode[5]="    ";
  gchar        *VelReference[3] = {"LSR", "Helio", "Observer"};
  gchar        *VelDef[2] = {"optical","radio"};
  gchar        RAString[16], DecString[16];
  odouble      freq;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *routine = "doSCAN";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  Obit_return_if_fail((inData->myDesc->ilocsu>=0), err,
		      "I can only work on multisource files");

  /* Get parameters */
  ObitInfoListGet(myInput, "doCrt",     &type, dim, &doCrt,  err);
  isInteractive = doCrt>0;
  if (isInteractive) { /* interactive, write to stdout */
    outStream = stdout;
  } else { /* write to file */
    ObitInfoListGetP(myInput, "prtFile",  &type, dim, (gpointer)&prtFile);
    if (prtFile) prtFile[MIN (47,dim[0])] = 0;
    ObitTrimTrail(prtFile);
    /* Make sure file named */
    Obit_return_if_fail(((strlen(prtFile)>2) && 
			 ((prtFile[0]!=' ') && (prtFile[1]!=' '))), err,
			"%s: Printer file not specified", routine);
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Writing output to %s",prtFile);
    ObitErrLog(err);
}
  if (err->error) Obit_traceback_msg (err, routine, "myInput");

  /* How long is a line? */
  lenLine = MAX (72, MIN(doCrt,1023));
  if (!isInteractive) lenLine = 132;
  count = 0;
  inDesc = inData->myDesc;  

  /* Convert SU table into Source List */
  ver = 1;
  if (inDesc->jlocif>=0) numIF = inDesc->inaxes[inDesc->jlocif];
  else numIF = 1;
  SUTable = newObitTableSUValue (inData->name, (ObitData*)inData, &ver, numIF, 
				 OBIT_IO_ReadOnly, err);
  if (SUTable) SouList = ObitTableSUGetList (SUTable, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  SUTable = ObitTableSUUnref(SUTable);

  /* Init vis counting */
  maxSUID = 1;
  for (i=0; i<SouList->number; i++) maxSUID = MAX (maxSUID, SouList->SUlist[i]->SourID);
  noVis = g_malloc0(maxSUID*sizeof(olong));

  /* Open NX Table */
  ver = 1;
  NXTable = newObitTableNXValue (inData->name, (ObitData*)inData, &ver, OBIT_IO_ReadOnly, err);
  /* Should be there */
  Obit_return_if_fail((NXTable!=NULL), err,
		      "%s: iNdeX table does not exist - use UV.PUtilIndex to create", routine);
  /* Open Index table  */
  iretCode = ObitTableNXOpen (NXTable, OBIT_IO_ReadOnly, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inData->name);
  /* Create row structure */
  NXRow = newObitTableNXRow(NXTable);

   /* Titles */
  sprintf (Title1, "                                                   ");
  sprintf (Title2, "Scan      Source      Qual  Calcode Sub         Timerange          FrqID  ");
  if (lenLine>95) {
    start = strlen(Title2);
    sprintf (&Title2[start], "Start Vis    End Vis");
  }

  /* Create/Open Printer */
  myPrint = ObitPrinterCreate (pgmName, isInteractive, outStream, prtFile);
  ObitPrinterOpen (myPrint, LinesPerPage, Title1, Title2, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Write header */
  if (!strncmp(Type,"AIPS",4)) {
    sprintf (line, "AIPS %s %s %d disk %d user %d",Aname,Aclass,Aseq, disk, AIPSuser);
  } else if (!strncmp(Type,"FITS",4)) {
    sprintf (line, "FITS %s disk %d",inFile, disk);
  } else {
    sprintf (line, "UNKNOWN Data type ");
  }
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  sprintf(line,"Object: %s Telescope: %s Observed: %s  Freq: %8.3lf GHz", 
	  inDesc->object,inDesc->teles, inDesc->obsdat, inDesc->freq*1.0e-9);
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  sprintf(line,"Scan summary listing");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  sprintf(line,"    ");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Titles */
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  
  /* Loop through Table */
  lastSouID = -999; 
  for (iRow=1; iRow<=NXTable->myDesc->nrow; iRow++) {
    iretCode = ObitTableNXReadRow (NXTable, iRow, NXRow, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Timerange */
    day2dhms (NXRow->Time-0.5*NXRow->TimeI, btimeString);
    day2dhms (NXRow->Time+0.5*NXRow->TimeI, etimeString);

    /* Source */
    souID = (olong)NXRow->SourID;
    /* Look up source */
    if (souID!=lastSouID) {
      lastSouID = souID;
      for (i=0; i<SouList->number; i++) {
	if (souID==SouList->SUlist[i]->SourID) {
	  qual = SouList->SUlist[i]->Qual;
	  strncpy (source, SouList->SUlist[i]->SourceName, 16);
	  strncpy (calcode, SouList->SUlist[i]->CalCode, 4);
	  break;
	}
      }
    }

    /* Sum visibilities */
    noVis[souID-1] += NXRow->EndVis-NXRow->StartVis+1;

    /* Format line  */
    sprintf(line,"%4d %16s : %4.4d %4s %4d %s - %s %3d",
	    count+1, source, qual, calcode, (olong)NXRow->SubA, btimeString, etimeString,
	    (olong)NXRow->FreqID);
    /* Add range of vis if room */
    if (lenLine>95) {
      start = strlen(line);
      sprintf (&line[start], " %10d %10d", (olong)NXRow->StartVis, (olong)NXRow->EndVis);
    }

    /* Print line */
    ObitPrinterWrite (myPrint, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    
    count++;
  } /* end loop reading NX Table */

  /* Close index table */
  iretCode = ObitTableNXClose (NXTable, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inData->name);

  /* Source info - new page */
  ObitPrinterSetTitle (myPrint, NULL, NULL, err);  /* No page titles at top */
  ObitPrinterNewPage (myPrint, &quit, err);
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    
  sprintf(line,"Source summary");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  /* Velocity type */
  ivr = MAX(0, MIN((inDesc->VelReference-1),2));
  ivd = MAX(0, MIN(inDesc->VelDef,1));
  sprintf(line,"Velocity type = '%s'    Definition = '%s'", VelDef[ivd], VelReference[ivd]);
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Titles */
  sprintf (Title1, "                                                   ");
  sprintf (Title2, "  ID Source           Qual  Calcode RA(2000.0)   Dec(2000.0)      ");
  if (lenLine>105) {
    start = strlen(Title2);
    sprintf (&Title2[start], "IFlux   QFlux   UFlux   VFlux  No. vis");
  }
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  ObitPrinterSetTitle (myPrint, Title1, Title2, err);  /* No page titles at top */
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Loop over sources */
  numIF = SouList->SUlist[0]->numIF;
  for (i=0; i<SouList->number; i++) {
    ObitPosLabelUtilRA2HMS  (SouList->SUlist[i]->RAMean,  inDesc->ctype[inDesc->jlocr], RAString);
    ObitPosLabelUtilDec2DMS (SouList->SUlist[i]->DecMean, inDesc->ctype[inDesc->jlocd], DecString);
    /* Loop over IFs */
    for (iif=0; iif<SouList->SUlist[i]->numIF; iif++) {
      /* Format line  */
      if (iif==0) { /* Only 1st IF */
	sprintf(line,"%4d %16s : %4.4d %4s %13s %13s ",
		SouList->SUlist[i]->SourID, SouList->SUlist[i]->SourceName, SouList->SUlist[i]->Qual,  
		SouList->SUlist[i]->CalCode, RAString, DecString);
      } else { /* Subsequent IFs */
	sprintf(line,"                                                               ");
      }
      /* Add Flux densities, no vis if room */
      if (lenLine>105) {
	start = strlen(line);
	if (iif==0)
	  sprintf (&line[start], "  %7.3f %7.3f %7.3f %7.3f%9d", 
		   SouList->SUlist[i]->IFlux[iif], SouList->SUlist[i]->QFlux[iif], 
		   SouList->SUlist[i]->UFlux[iif], SouList->SUlist[i]->VFlux[iif], 
		   noVis[SouList->SUlist[i]->SourID-1]);
	else
 	  sprintf (&line[start], " %7.3f %7.3f %7.3f %7.3f", 
		   SouList->SUlist[i]->IFlux[iif], SouList->SUlist[i]->QFlux[iif], 
		   SouList->SUlist[i]->UFlux[iif], SouList->SUlist[i]->VFlux[iif]);
     }
      /* Print line */
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    } /* end loop over IFs */
 } /* end loop over sources */

  /* Source frequency info */
  /* Titles */
  sprintf (Title1, "                                                   ");
  sprintf (Title2, "  ID Source            Freq(GHz) Velocity(Km/s) Rest freq (GHz)");
  ObitPrinterSetTitle (myPrint, Title1, Title2, err);  /* titles at top */
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  if (quit) goto Quit;
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Loop over sources */
  for (i=0; i<SouList->number; i++) {
     /* Loop over IFs */
    for (iif=0; iif<SouList->SUlist[i]->numIF; iif++) {
      if (inDesc->jlocif>=0) 
	freq = inDesc->freqIF[iif] + SouList->SUlist[i]->FreqOff[iif];
      else
	freq = inDesc->freq + SouList->SUlist[i]->FreqOff[iif];
      /* Format line  */
      if (iif==0) { /* Only 1st IF */
	sprintf(line,"%4d %16s %10.4lf %10.4lf %14.4lf",
		SouList->SUlist[i]->SourID, SouList->SUlist[i]->SourceName, freq*1.0e-9, 
		SouList->SUlist[i]->LSRVel[iif]*1.0e-3, 
		SouList->SUlist[i]->RestFreq[iif]*1.0e-9);
      } else { /* Subsequent IFs */
	sprintf(line,"     IF(%3d)          %10.4lf %10.4lf %14.4lf",
		iif+1, freq*1.0e-9, 
		SouList->SUlist[i]->LSRVel[iif]*1.0e-3, 
		SouList->SUlist[i]->RestFreq[iif]*1.0e-9);
      }
      /* Print line */
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    } /* end loop over IFs */
  } /* End loop over sources */

  /* Frequency info */
  /* Source info - new page */
  ObitPrinterSetTitle (myPrint, NULL, NULL, err);  /* No page titles at top */
  ObitPrinterNewPage (myPrint, &quit, err);
  if (quit) goto Quit;
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    
  sprintf(line,"Frequency Table summary");
  ObitPrinterWrite (myPrint, line, &quit, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
  ObitPrinterSetTitle (myPrint, NULL, NULL, err);  /* No page titles at top */
  /* Titles */
  sprintf (Title1, "                                                   ");
  sprintf (Title2, "FQID IF#      Freq(GHz)      BW(kHz)   Ch.Sep(kHz)  Sideband");
  ObitPrinterWrite (myPrint, Title1, &quit, err);
  ObitPrinterWrite (myPrint, Title2, &quit, err);
  ObitPrinterSetTitle (myPrint, Title1, Title2, err);  /* Page titles at top */
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Open FQ Table */
  ver = 1;
  FQTable = newObitTableFQValue (inData->name, (ObitData*)inData, &ver, numIF, 
				 OBIT_IO_ReadOnly, err);
  /* Should be there */
  Obit_return_if_fail((FQTable!=NULL), err,
		      "%s: FQ table does not exist", routine);
  /* Open/close Index table  */
  iretCode = ObitTableFQOpen (FQTable, OBIT_IO_ReadOnly, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, FQTable->name);
  /* Create row structure */
  FQRow = newObitTableFQRow(FQTable);

  /* Loop over table */
  for (iRow=1; iRow<=FQTable->myDesc->nrow; iRow++) {
     iretCode = ObitTableFQReadRow (FQTable, iRow, FQRow, err);
    if (err->error) Obit_traceback_msg (err, routine, FQTable->name);
     /* Loop over IFs */
    for (iif=0; iif<FQTable->numIF; iif++) {
      /* Frequency */
      if (inDesc->jlocif>=0) 
	freq = inDesc->freqIF[iif]; /* + FQRow->freqOff[iif];*/
      else
	freq = inDesc->freq + FQRow->freqOff[iif];
      /* Format line  */
      if (iif==0) { /* Only 1st IF */
	sprintf(line,"%4d %3d %14.8lf %14.4f %10.4f %6d",
		FQRow->fqid, iif+1, freq*1.0e-9, 
		FQRow->totBW[iif]*1.0e-3,  FQRow->chWidth[iif]*1.0e-3, 
		(olong)FQRow->sideBand[iif]);
      } else { /* Subsequent IFs */
	sprintf(line,"     %3d %14.8lf %14.4f %10.4f %6d",
		iif+1, freq*1.0e-9, 
		FQRow->totBW[iif]*1.0e-3,  FQRow->chWidth[iif]*1.0e-3, 
		(olong)FQRow->sideBand[iif]);
      }
      /* Print line */
      ObitPrinterWrite (myPrint, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
    } /* end loop over IFs */
  } /* End loop over table */

   /* Close FQ */
  iretCode = ObitTableFQClose (FQTable, err);

  /* Done - Close Printer */
 Quit:
  ObitPrinterClose (myPrint, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Cleanup */
  NXTable = ObitTableNXUnref(NXTable);
  NXRow   = ObitTableNXRowUnref(NXRow);
  FQTable = ObitTableFQUnref(FQTable);
  FQRow   = ObitTableFQRowUnref(FQRow);
  SouList = ObitSourceListUnref(SouList);
  if (noVis) g_free(noVis);

} /* end doSCAN */

/*----------------------------------------------------------------------- */
/*  Gain table listing                                                    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input UV data                                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doGAIN (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitPrinter     *myPrint = NULL;
  ObitSourceList  *SouList = NULL;
  ObitTableSU     *SUTable = NULL;
  ObitAntennaList *AntList = NULL;
  ObitTableAN     *ANTable = NULL;
  ObitTableCL     *CLTable = NULL;
  ObitTableCLRow  *CLRow   = NULL;
  ObitTableSN     *SNTable = NULL;
  ObitTableSNRow  *SNRow   = NULL;
  ObitUVSel       *sel;
  ObitIOCode   iretCode;
  ObitUVDesc   *inDesc;
  FILE         *outStream = NULL;
  gboolean     isInteractive = FALSE, quit = FALSE;
  gboolean     doSN, done=FALSE, souChange=FALSE, firstPol, *doAnt=NULL;
  gchar        line[1024], Title1[1024], Title2[1024];
  oint         numPol, numIF, numTerm, numPCal, numOrb;
  olong        BIF, i, ia, dt, ipass, npass, ndig, nfit, lenLine, LinesPerPage = 0;
  olong        ver, iRow, maxAnt, nAnt, mAnt, iif, ipol, start, gainVer, doCrt=1;
  olong        loAnt, hiAnt, nrow, souID=0, SubA, antNo=0, freqID, lastSouID, lastSou;
  ofloat       time=0., lasttime, value=0., *valueArr=NULL, fblank = ObitMagicF();
  ofloat       scale = 1.0;
  gchar        *prtFile=NULL, timeString[25], inTab[28], dispType[10], source[20];
  gchar        *dTypes[] = {"AMP     ","PHASE   ","WT      ","DELAY   ","RATE    "};
  gchar        *dLabel[] = {"Amplitude","Phase","Weight/SNR","Delay","Rate"};
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *routine = "doSCAN";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get parameters */
  ObitInfoListGet(myInput, "doCrt",     &type, dim, &doCrt,  err);
  isInteractive = doCrt>0;
  if (isInteractive) { /* interactive, write to stdout */
    outStream = stdout;
  } else { /* write to file */
    ObitInfoListGetP(myInput, "prtFile",  &type, dim, (gpointer)&prtFile);
    if (prtFile) prtFile[MIN (47,dim[0])] = 0;
    ObitTrimTrail(prtFile);
    /* Make sure file named */
    Obit_return_if_fail(((strlen(prtFile)>2) && 
			 ((prtFile[0]!=' ') && (prtFile[1]!=' '))), err,
			"%s: Printer file not specified", routine);
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Writing output to %s",prtFile);
    ObitErrLog(err);
  }
  ObitInfoListGet(myInput, "gainVer",   &type, dim, &gainVer, err);
  ObitInfoListGet(myInput, "BIF",       &type, dim, &BIF, err);
  BIF = MAX (1, BIF);
  ObitInfoListGet(myInput, "ndig",      &type, dim, &ndig, err);
  ndig = MAX(4, MIN(8, ndig));
  ObitInfoListGet(myInput, "inTab",     &type, dim, inTab, err);
  inTab[dim[0]] = 0;
  /* Default "AIPS CL" */
  if (!strncmp (inTab, "    ", 4)) sprintf (inTab, "AIPS CL");
  ObitTrimTrail(inTab);
  ObitInfoListGet(myInput, "dispType",  &type, dim, dispType, err);
  dispType[dim[0]] = 0;
  if (err->error) Obit_traceback_msg (err, routine, "myInput");

  /* Digest */
  doSN = (!strncmp (inTab, "AIPS SN", 7));  /* Use SN? */
  dt = -1;       /* Display type code */
  /* Default is "AMP" */
  if (!strncmp(dispType, "    ",4)) sprintf (dispType,"%s",dTypes[0]);
  for (i=0; i<5; i++) {
    if (!strcmp (dispType, dTypes[i])) {
      dt = i;
      break;
    }
  }
  Obit_return_if_fail((dt>=0), err, "Unknown dispType %s",dispType );

  /* How long is a line? */
  lenLine = MAX (72, MIN(doCrt,1023));
  if (!isInteractive) lenLine = 132;

  /* Open uv data - to use for data selection */
  iretCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inData->name);
  inDesc = inData->myDesc;  
  sel    = inData->mySel;
  /* Convert SU table into Source List */
  ver = 1;
  if (inDesc->jlocif>=0) numIF = inDesc->inaxes[inDesc->jlocif];
  else numIF = 1;
  BIF = MIN (BIF, numIF);
  SUTable = newObitTableSUValue (inData->name, (ObitData*)inData, &ver, numIF, 
				 OBIT_IO_ReadOnly, err);
  if (SUTable) SouList = ObitTableSUGetList (SUTable, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  SUTable = ObitTableSUUnref(SUTable);

  /* Convert AN table into AntennaList */
  ver = 1;
  numPCal  = 0;
  numOrb   = 0;
  ANTable = newObitTableANValue (inData->name, (ObitData*)inData, &ver, 
				 numOrb, numPCal, OBIT_IO_ReadOnly, err);
  if (ANTable) AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ANTable = ObitTableANUnref(ANTable);

  /* List of antennas to display */
  maxAnt = 0;
  for (i=0; i<AntList->number; i++) maxAnt = MAX (maxAnt, AntList->ANlist[i]->AntID);
  doAnt    = g_malloc0((maxAnt+5)*sizeof(gboolean));
  valueArr = g_malloc0((maxAnt+5)*sizeof(ofloat));
  for (i=0; i<maxAnt; i++) doAnt[i]    = FALSE;
  for (i=0; i<maxAnt; i++) valueArr[i] = fblank;
  nAnt = 0;
  for (i=0; i<AntList->number; i++) {
    doAnt[i] = ObitUVSelWantAnt (sel, AntList->ANlist[i]->AntID);
    if (doAnt[i]) nAnt++;   /* How many to display? */
  }
  
  /* Open CL/SN Table */
  ver = gainVer;
  numPol  = inDesc->inaxes[inDesc->jlocs];
  numTerm = 0;
  if (doSN) {  /* SN table */
    SNTable = newObitTableSNValue (inData->name, (ObitData*)inData, &ver, OBIT_IO_ReadOnly, 
				   numPol, numIF, err);
    gainVer = ver;
    /* Should be there */
    Obit_return_if_fail((SNTable!=NULL), err, "SN table %d does not exist", gainVer);
    /* Sort input SN if needed */
    ObitTableUtilSort ((ObitTable*)SNTable, "TIME    ", FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    /* Open  table  */
    iretCode = ObitTableSNOpen (SNTable, OBIT_IO_ReadOnly, err);
    if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_msg (err, routine, inData->name);
    /* Create row structure */
    SNRow = newObitTableSNRow(SNTable);
    nrow = SNTable->myDesc->nrow;   /* How many rows */
  } else { /* CL table */
    CLTable = newObitTableCLValue (inData->name, (ObitData*)inData, &ver, OBIT_IO_ReadOnly, 
				   numPol, numIF, numTerm, err);
    gainVer = ver;
    /* Should be there */
    Obit_return_if_fail((CLTable!=NULL), err, "CL table %d does not exist", gainVer);
    /* Open  table  */
    iretCode = ObitTableCLOpen (CLTable, OBIT_IO_ReadOnly, err);
    if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback, return */
      Obit_traceback_msg (err, routine, inData->name);
    /* Create row structure */
    CLRow = newObitTableCLRow(CLTable);
    nrow = CLTable->myDesc->nrow;   /* How many rows */
  }

  /* How many will fit on a row per pass? */
  nfit = (lenLine - 20)/ndig;

  /* How many passes? */
  npass = nAnt/nfit;
  if ((nAnt%nfit)!=0) npass++;
  
  /* Create/Open Printer */
  myPrint = ObitPrinterCreate (pgmName, isInteractive, outStream, prtFile);
  ObitPrinterOpen (myPrint, LinesPerPage, Title1, Title2, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	
  /* Loop in IF */
  if (inDesc->jlocif>=0) numIF = inDesc->inaxes[inDesc->jlocif];
  else numIF = 1;
  for (iif=BIF-1; iif<BIF+numIF-1; iif++) {
    
    /* Loop in poln */
    numPol  = inDesc->inaxes[inDesc->jlocs];
    for (ipol=0; ipol<numPol; ipol++) {
      /* Is this the first poln in the table? */
      firstPol = (ipol==0) && (inDesc->crval[inDesc->jlocs]>-1.5);
      
      /* Loop making passes */
      loAnt = 0; hiAnt = -1;
      for (ipass=0; ipass<npass; ipass++) {

	/* May need to start on a new page */
	if ((iif>BIF-1) || (ipol>0) || (ipass>0)) {
	  /* New page */
	  ObitPrinterSetTitle (myPrint, NULL, NULL, err);  /* No page titles at top */
	  ObitPrinterNewPage (myPrint, &quit, err);
	  if (quit) goto Quit;
	  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	}


	/* Write header */
	if (!strncmp(Type,"AIPS",4)) {
	  sprintf (line, "AIPS %s %s %d disk %d user %d",Aname,Aclass,Aseq, disk, AIPSuser);
	} else if (!strncmp(Type,"FITS",4)) {
	  sprintf (line, "FITS %s disk %d",inFile, disk);
	} else {
	  sprintf (line, "UNKNOWN Data type ");
	}
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	
	sprintf(line,"Object: %s Telescope: %s Observed: %s  Freq: %8.3lf GHz", 
		inDesc->object,inDesc->teles, inDesc->obsdat, inDesc->freq*1.0e-9);
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	
	sprintf(line,"Gain table listing %s version %d", inTab, gainVer);
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	
	sprintf(line,"    ");
	ObitPrinterWrite (myPrint, line, &quit, err);
	if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	
	/* Set antenna range this pass */
	mAnt = 0;
	loAnt = -1;
	for (ia=hiAnt+1; ia<=maxAnt; ia++) {
	  if (doAnt[ia]) { /* Good one */
	    if (loAnt<0) loAnt = ia;
	    hiAnt = ia;
	    mAnt++;
	    if (mAnt>=nfit) break;
	  }
	} /* end loop finding antenna range */

	/* Set scaling */
	if (doSN)
	  scale = getSNGainScale (sel, SNTable, SNRow, maxAnt, doAnt, loAnt, hiAnt, 
				  iif, firstPol, dt, ndig, err);
	else 
	  scale = getCLGainScale (sel, CLTable, CLRow, maxAnt, doAnt, loAnt, hiAnt, 
				  iif, firstPol, dt, ndig, err);
	
	/* Titles */
	if (firstPol)
	  sprintf (Title1, "IF %d RCP: %s scaled by %g", iif+1,dLabel[dt], scale);
	else
	  sprintf (Title1, "IF %d LCP: %s scaled by %g", iif+1,dLabel[dt], scale);
	sprintf (Title2, "    Time        Source         ");
	for (ia=loAnt; ia<=hiAnt; ia++) {
	  if (doAnt[ia]) {
	    start = strlen(Title2);
	    sprintf (&Title2[start], "-%*d",ndig,ia+1);
	  }
	}
	ObitPrinterSetTitle (myPrint, Title1, Title2, err);  /* Page titles at top */
	ObitPrinterWrite (myPrint, Title1, &quit, err);
	ObitPrinterWrite (myPrint, Title2, &quit, err);
	if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

	/* Loop through Table */
	lastSouID = -999; 
	lastSou   = -999; 
	lasttime  = -1.0e20;
	/* Extra pass through loop for write final */
	for (iRow=1; iRow<=nrow+1; iRow++) {
	  done = (iRow>nrow);
	  /* Read */
	  if (!done) {
	    if (doSN) {   /* SN table */
	      iretCode = ObitTableSNReadRow (SNTable, iRow, SNRow, err);
	      time   = SNRow->Time;
	      souID  = SNRow->SourID;
	      freqID = SNRow->FreqID;
	      antNo  = SNRow->antNo-1;
	      SubA   = SNRow->SubA;
	      value  = getSNValue(SNRow, firstPol, iif, dt);
	    } else {  /* CL table */
	      iretCode = ObitTableCLReadRow (CLTable, iRow, CLRow, err);
	      time   = CLRow->Time;
	      souID  = CLRow->SourID;
	      freqID = CLRow->FreqID;
	      antNo  = CLRow->antNo-1;
	      SubA   = CLRow->SubA;
	      value  = getCLValue(CLRow, firstPol, iif, dt);
	    } /* End CL table */
	    if (err->error) Obit_traceback_msg (err, routine, inData->name);
	  
	    /* Initialize lasttime, lastSou */
	    if (lasttime<-1.0e20) lasttime = time;
	    if (lastSou<0)        lastSou  = souID;
	    
	    /* Want this antenna? */
	    if (!doAnt[antNo] || (antNo<loAnt) || (antNo>hiAnt)) continue;
	    /* Want this source? */
	    if (!ObitUVSelWantSour (sel, souID)) continue;
	    /* Want this Subarray? */
	    if ((sel->SubA>0) && (SubA!=sel->SubA)) continue;
	    /* Want this FreqID? */
	    if ((sel->FreqID>0) && (freqID!=sel->FreqID)) continue;
	    /* Valid value? */
	    if (value==fblank) continue;
	    
	  } /* end read if not done */

	  /* If time greater than last write current line */
	  if ((time>lasttime) || done) {
	    day2dhms(lasttime, timeString);
	    /* Format line  */
	    sprintf(line,"%13s %16s ", timeString, source);
	    /* Add values */
	    for (ia=loAnt; ia<=hiAnt; ia++) {
	      if (!doAnt[ia]) continue;
	      start = strlen(line);
	      if (valueArr[ia]!=fblank) {
		/* Encode as scaled value with ndig digits */
		sprintf (&line[start], " %*.0f", ndig, scale*valueArr[ia]);
		valueArr[ia] = fblank;
	      } else {
		/* Blank fill if in value */
		for (i=start; i<=start+ndig; i++) line[i] = ' '; line[i] = 0;
	      }
	    }
	    
	    /* Blank line if changing source */
	    if (souChange) {
	      souChange = FALSE;
	      sprintf (timeString, "    ");
	      ObitPrinterWrite (myPrint, timeString, &quit, err);
	      if (quit) goto Quit;
	      if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	      lastSou = souID;
	    }
	    
	    /* Print line */
	    ObitPrinterWrite (myPrint, line, &quit, err);
	    if (quit) goto Quit;
	    if (err->error) Obit_traceback_msg (err, routine, myPrint->name);
	    lasttime = time;  /* Reset last time */
	    if (done) break;
	  } /* End write */

	  /* Accumulate */
	  valueArr[antNo] = value;
	  /* Look up source */
	  if (souID!=lastSouID) {
	    if(lastSouID>-99) souChange = TRUE;
	    lastSouID = souID;
	    sprintf (source,"Unspecified     ");
	    for (i=0; i<SouList->number; i++) {
	      if (souID==SouList->SUlist[i]->SourID) {
		strncpy (source, SouList->SUlist[i]->SourceName, 16);
		break;
	      }
	    }
	  } /* end accumulate */
	  
	} /* end loop reading gain Table */
      } /* End passes loop */
    } /* end poln loop */
  } /* End IF loop */
    /* Done - Close Printer */
 Quit:
  ObitPrinterClose (myPrint, err);
  ObitUVClose(inData, err);
  if (err->error) Obit_traceback_msg (err, routine, myPrint->name);

  /* Cleanup */
  CLTable = ObitTableCLUnref(CLTable);
  CLRow   = ObitTableCLRowUnref(CLRow);
  SNTable = ObitTableSNUnref(SNTable);
  SNRow   = ObitTableSNRowUnref(SNRow);
  SouList = ObitSourceListUnref(SouList);
  AntList = ObitAntennaListUnref(AntList);
  if (doAnt)    g_free(doAnt);
  if (valueArr) g_free(valueArr);

} /* end doGAIN */

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

/** 
 * Extract value from an SN table row
 * \param SNRow     SN Table row
 * \param firstPol  If true, first poln, else second
 * \param iif       0-rel IF number
 * \param dtype     Desired data type code
 *                  0=amp, 1=phase, 3=wt, 4=Delay, 5=Rate
 * \return value, blanked if necessary
 */
ofloat getSNValue (ObitTableSNRow *SNRow, gboolean firstPol, 
		   olong iif, olong dtype)
{
  ofloat       value, fblank = ObitMagicF();

  value = fblank;  /* Initial value */
  /* Data by type */
  switch (dtype) {
  case 0:  /* "AMP"   */
    if (firstPol) {
      if ((SNRow->Real1[iif]!=fblank) && (SNRow->Imag1[iif]!=fblank) && (SNRow->Weight1[iif]>0.0))
	value = sqrt (SNRow->Real1[iif]*SNRow->Real1[iif] + 
		      SNRow->Imag1[iif]*SNRow->Imag1[iif]);
    } else {
      if ((SNRow->Real2[iif]!=fblank) && (SNRow->Imag2[iif]!=fblank) && (SNRow->Weight2[iif]>0.0))
	value = sqrt (SNRow->Real2[iif]*SNRow->Real2[iif] + 
		      SNRow->Imag2[iif]*SNRow->Imag2[iif]);
    }
    break;
  case 1:  /* "PHASE" */
    if (firstPol) {
      if ((SNRow->Real1[iif]!=fblank) && (SNRow->Imag1[iif]!=fblank) && (SNRow->Weight1[iif]>0.0))
	value = RAD2DG * atan2 (SNRow->Imag1[iif], SNRow->Real1[iif]+1.0e-20);
    } else {
      if ((SNRow->Real2[iif]!=fblank) && (SNRow->Imag2[iif]!=fblank) && (SNRow->Weight2[iif]>0.0))
	value = RAD2DG * atan2 (SNRow->Imag2[iif], SNRow->Real2[iif]+1.0e-20);
    }
  case 2:  /* "WT"    */
    if (firstPol) {
      value = SNRow->Weight1[iif];
    } else {
      value = SNRow->Weight2[iif];
    }
    break;
  case 3:  /* "DELAY" */
    if (firstPol) {
      if ((SNRow->Real1[iif]!=fblank) && (SNRow->Imag1[iif]!=fblank) && (SNRow->Weight1[iif]>0.0))
	value = SNRow->Delay1[iif];
      
    } else {
      if ((SNRow->Real2[iif]!=fblank) && (SNRow->Imag2[iif]!=fblank) && (SNRow->Weight2[iif]>0.0))
	value = SNRow->Delay2[iif];
    }
    break;
  case 4:  /* "RATE"  */
    if (firstPol) {
      if ((SNRow->Real1[iif]!=fblank) && (SNRow->Imag1[iif]!=fblank) && (SNRow->Weight1[iif]>0.0))
	value = SNRow->Rate1[iif];
      
    } else {
      if ((SNRow->Real2[iif]!=fblank) && (SNRow->Imag2[iif]!=fblank) && (SNRow->Weight2[iif]>0.0))
	value = SNRow->Rate2[iif];
    }
    break;
  default:
    value = fblank;
  }; /* end switch by type */
  return value;
} /* end getSNValue */

/** 
 * Extract value from an CL table row
 * \param CLRow     CL Table row
 * \param firstPol  If true, first poln, else second
 * \param iif       0-rel IF number
 * \param dtype     Desired data type code
 *                  0=amp, 1=phase, 3=wt, 4=Delay, 5=Rate
 * \return value, blanked if necessary
 */
ofloat getCLValue (ObitTableCLRow *CLRow, gboolean firstPol, 
		   olong iif, olong dtype)
{
  ofloat       value, fblank = ObitMagicF();

  value = fblank;  /* Initial value */
  /* Data by type */
  switch (dtype) {
  case 0:  /* "AMP"   */
    if (firstPol) {
      if ((CLRow->Real1[iif]!=fblank) && (CLRow->Imag1[iif]!=fblank) && (CLRow->Weight1[iif]>0.0))
	value = sqrt (CLRow->Real1[iif]*CLRow->Real1[iif] + 
		      CLRow->Imag1[iif]*CLRow->Imag1[iif]);
    } else {
      if ((CLRow->Real2[iif]!=fblank) && (CLRow->Imag2[iif]!=fblank) && (CLRow->Weight2[iif]>0.0))
	value = sqrt (CLRow->Real2[iif]*CLRow->Real2[iif] + 
		      CLRow->Imag2[iif]*CLRow->Imag2[iif]);
    }
    break;
  case 1:  /* "PHASE" */
    if (firstPol) {
      if ((CLRow->Real1[iif]!=fblank) && (CLRow->Imag1[iif]!=fblank) && (CLRow->Weight1[iif]>0.0))
	value = RAD2DG * atan2 (CLRow->Imag1[iif], CLRow->Real1[iif]+1.0e-20);
    } else {
      if ((CLRow->Real2[iif]!=fblank) && (CLRow->Imag2[iif]!=fblank) && (CLRow->Weight2[iif]>0.0))
	value = RAD2DG * atan2 (CLRow->Imag2[iif], CLRow->Real2[iif]+1.0e-20);
    }
    break;
  case 2:  /* "WT"    */
    if (firstPol) {
      value = CLRow->Weight1[iif];
    } else {
      value = CLRow->Weight2[iif];
    }
    break;
  case 3:  /* "DELAY" */
    if (firstPol) {
      if ((CLRow->Real1[iif]!=fblank) && (CLRow->Imag1[iif]!=fblank) && (CLRow->Weight1[iif]>0.0))
	value = CLRow->Delay1[iif];
      
    } else {
      if ((CLRow->Real2[iif]!=fblank) && (CLRow->Imag2[iif]!=fblank) && (CLRow->Weight2[iif]>0.0))
	value = CLRow->Delay2[iif];
    }
    break;
  case 4:  /* "RATE"  */
    if (firstPol) {
      if ((CLRow->Real1[iif]!=fblank) && (CLRow->Imag1[iif]!=fblank) && (CLRow->Weight1[iif]>0.0))
	value = CLRow->Rate1[iif];
      
    } else {
      if ((CLRow->Real2[iif]!=fblank) && (CLRow->Imag2[iif]!=fblank) && (CLRow->Weight2[iif]>0.0))
	value = CLRow->Rate2[iif];
    }
    break;
  default:
    value = fblank;
  }; /* end switch by type */
  return value;
} /* end getCLValue */

/** 
 * Get scaling for data in an SN table
 * \param SNTable   SN Table
 * \param SNRow     SN Table row
 * \param maxAnt    Maximum antenna number
 * \param doAnt     0-rel indexed got antenna array
 * \param loAnt     0-rel lowest antenna number to check
 * \param hiAnt     0-rel highest antenna number to check
 * \param iif       0-rel IF number
 * \param firstPol  If true, first poln, else second
 * \param dtype     Desired data type code
 *                  0=amp, 1=phase, 3=wt, 4=Delay, 5=Rate
 * \param ndig      Number of digits in display
 * \param err       Obit Error stack 
 * \return scaling value
 */
ofloat getSNGainScale (ObitUVSel *sel, ObitTableSN *SNTable, ObitTableSNRow *SNRow, 
		       olong maxAnt, gboolean *doAnt, olong loAnt, olong hiAnt, 
		       olong iif, gboolean firstPol, olong dtype, olong ndig, 
		       ObitErr *err)
{
  ofloat out = 1.0;
  ofloat value, maxValue=-1.0e20, ftest, ratio, fblank = ObitMagicF();
  olong nrow, iRow, count=0, antNo, souID, SubA, freqID, ilog;
  gboolean isNeg = FALSE;
  ObitIOCode   iretCode;
  gchar *routine = "getSNGainScale";

  nrow = SNTable->myDesc->nrow;   /* How many rows */
  /* Loop over table or first 1000 hits */
  for (iRow=1; iRow<=nrow; iRow++) {
    iretCode = ObitTableSNReadRow (SNTable, iRow, SNRow, err);
    if (err->error) Obit_traceback_val (err, routine, SNTable->name, out);

    /* Want this antenna? */
    antNo  = SNRow->antNo-1;
    if (!doAnt[antNo] || (antNo<loAnt) || (antNo>hiAnt)) continue;
	
    /* Want this source? */
    souID  = SNRow->SourID;
    if (!ObitUVSelWantSour (sel, souID)) continue;

    /* Want this Subarray? */
    SubA   = SNRow->SubA;
    if ((sel->SubA>0) && (SubA!=sel->SubA)) continue;

    /* Want this FreqID? */
    freqID = SNRow->FreqID;
    if ((sel->FreqID>0) && (freqID!=sel->FreqID)) continue;

    /* Valid value? */
    value  = getSNValue(SNRow, firstPol, iif, dtype);
    if (value==fblank) continue;

    /* Keep max abs value */
    if (maxValue==fblank) {maxValue = fabs(value); isNeg = (value<0.0);}
    else if (fabs(value)>maxValue) 
      {maxValue = fabs(value); isNeg = (value<0.0);}
    /* Done enough? */
    count++;
    if (count>1000) break;
  } /* end loop */

  /* Get scaling - first target max */
  if (isNeg) {  /* Negative */
    ftest = pow(10.0, (odouble)(ndig-2.0));
  } else {      /* Positive */
    ftest = pow(10.0, (odouble)(ndig-1.0));
  }
   if (maxValue<1.0e-20) {
    out = 1.0;
  } else {
    ratio = ftest/maxValue;
    ilog = (olong)log10(ratio);
    out = pow(10.0, (odouble)ilog);
  }
  return out;
} /*  end getSNGainScale */

/** 
 * Get scaling for data in an CL table
 * \param sel       ObitUVSel to decide which entries wanted
 * \param CLTable   CL Table
 * \param CLRow     CL Table row
 * \param maxAnt    Maximum antenna number
 * \param doAnt     0-rel indexed got antenna array
 * \param loAnt     0-rel lowest antenna number to check
 * \param hiAnt     0-rel highest antenna number to check
 * \param iif       0-rel IF number
 * \param firstPol  If true, first poln, else second
 * \param dtype     Desired data type code
 *                  0=amp, 1=phase, 3=wt, 4=Delay, 5=Rate
 * \param ndig      Number of digits in display
 * \param err       Obit Error stack 
 * \return scaling value
 */
ofloat getCLGainScale (ObitUVSel *sel, ObitTableCL *CLTable, ObitTableCLRow *CLRow, 
		       olong maxAnt, gboolean *doAnt, olong loAnt, olong hiAnt, 
		       olong iif, gboolean firstPol, olong dtype, olong ndig, 
		       ObitErr *err)
{
  ofloat out = 1.0;
  ofloat value, maxValue=-1.0e20, ftest, ratio, fblank = ObitMagicF();
  olong nrow, iRow, count=0, antNo, souID, SubA, freqID, ilog;
  gboolean isNeg = FALSE;
  ObitIOCode   iretCode;
  gchar *routine = "getCLGainScale";

  nrow = CLTable->myDesc->nrow;   /* How many rows */
  /* Loop over table or first 1000 hits */
  for (iRow=1; iRow<=nrow; iRow++) {
    iretCode = ObitTableCLReadRow (CLTable, iRow, CLRow, err);
    if (err->error) Obit_traceback_val (err, routine, CLTable->name, out);

    /* Want this antenna? */
    antNo  = CLRow->antNo-1;
    if (!doAnt[antNo] || (antNo<loAnt) || (antNo>hiAnt)) continue;
	
    /* Want this source? */
    souID  = CLRow->SourID;
    if (!ObitUVSelWantSour (sel, souID)) continue;

    /* Want this Subarray? */
    SubA   = CLRow->SubA;
    if ((sel->SubA>0) && (SubA!=sel->SubA)) continue;

    /* Want this FreqID? */
    freqID = CLRow->FreqID;
    if ((sel->FreqID>0) && (freqID!=sel->FreqID)) continue;

    /* Valid value? */
    value  = getCLValue(CLRow, firstPol, iif, dtype);
    if (value==fblank) continue;

    /* Keep max abs value */
    if (maxValue==fblank) {maxValue = fabs(value); isNeg = (value<0.0);}
    else if (fabs(value)>maxValue) 
      {maxValue = fabs(value); isNeg = (value<0.0);}
    /* Done enough? */
    count++;
    if (count>1000) break;
  } /* end loop */

  /* Get scaling - first target max */
  if (isNeg) {  /* Negative */
    ftest = pow(10.0, (odouble)(ndig-2.0));
  } else {      /* Positive */
    ftest = pow(10.0, (odouble)(ndig-1.0));
  }
  if (maxValue<1.0e-20) {
    out = 1.0;
  } else {
    ratio = ftest/maxValue;
    ilog = (olong)log10(ratio);
    out = pow(10.0, (odouble)ilog);
  }

  return out;
} /*  end getCLGainScale */
