/* $Id$  */
/* Differential instrumental polarization calibration       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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
/*; Correspondence about this software should be addressed as follows:*/
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
#include "ObitTableCCUtil.h"
#include "ObitHistory.h"
#include "ObitData.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/* internal prototypes */
/* Get inputs */
ObitInfoList* DPCalIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void DPCalOut (ObitInfoList* outList, ObitErr *err);
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
ObitUV* setOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Get calibration solutions */
void GetCalSoln (ObitInfoList *myInput, ObitUV* inData, ObitErr *err);
/* Apply calibrator solutions */
void ApplyCalSoln (ObitInfoList *myInput, ObitUV* inData, ObitUV* outData, ObitErr *err);
/* Write history */
void DPCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err);

/* Program globals */
gchar *pgmName = "DPCal";       /* Program name */
gchar *infile  = "DPCal.in" ;   /* File with program inputs */
gchar *outfile = "/tmp/DPCal.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
olong  nIF;                    /* Number of IFs in calData */
olong  nTime;                  /* Number of times in calData */
ofloat *calData        = NULL; /* Calibrator data */
ofloat *timeData       = NULL; /* Times of calibrator data */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Differential instrumental polarization calibration                   */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line, read inputs */
  err = newObitErr();
  myInput = DPCalIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest inputs */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get output uvdata */
  outData = setOutputData (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get calibration solutions */
  GetCalSoln (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Apply calibrator solutions */
  ApplyCalSoln (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* History */
  DPCalHistory (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
  if (calData)  g_free(calData);
  if (timeData) g_free(timeData);
 
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* DPCalIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "DPCalIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
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
      
     } else if (strcmp(arg, "-DataType") == 0) { /* Data type AIPS or FITS */
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
      
    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS output UV sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outDisk") == 0) { /* output UV disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-outName") == 0) { /* AIPS UV outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* AIPS UV outClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

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
} /* end DPCalIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: DPCal -input file -output ofile [args]\n");
    fprintf(stderr, "DPCal Differential inst. poln cal.\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def DPCal.in\n");
    fprintf(stderr, "  -output output result file, def DPCal.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output uv FITS  file\n");  
    fprintf(stderr, "  -outName output uv AIPS file name\n");
    fprintf(stderr, "  -outClass output uv AIPS file class\n");
    fprintf(stderr, "  -outSeq output uv AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output uv (AIPS or FITS) disk number (1-rel) \n");
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
/*     inFile    Str [?]    input FITS uv file name [no def]              */
/*     inDisk    Int        input AIPS or FITS uv disk no  [def 1]        */
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     BIF       Int [1]    first IF to process, 0=>1                     */
/*     EIF       Int [1]    highest IF to process, 0=>all                 */
/*     outDisk   Int        output AIPS or FITS image disk no  [def 1]    */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     outName   Str [12]   output AIPS image name  [no def]              */
/*     outClass  Str [6]    output AIPS image class  [no def]             */
/*     outSeq    Int        output AIPS image sequence no  [no def]       */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
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
  strTemp = "DPCal.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "DPCalName";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS input uv sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* AIPS or FITS input uv disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* BIF  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* EIF  */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "EIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input FITS UV file name root*/
  strTemp = "DPCalModel.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2File", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS UV file name */
  strTemp = "UV";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Name", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS UV file class root */
  strTemp = "IMAGER";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "in2Class", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Input AIPS UV sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "in2Seq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS UV file name */
  strTemp = "DPCalOut.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file name */
  strTemp = "DPCalOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file class */
  strTemp = "Class ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS UV disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
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
  /*ObitInfoType type; */
    gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    gboolean btemp;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /*  Output file not expected to exist */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (myInput, "outExist", OBIT_bool, dim, &btemp);

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
  olong        nvis=1000;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "Stokes", "timeRange", "BChan", "EChan",  "BIF", "EIF", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "Smooth", "Antennas",  "Sources",  "souCode", "Qual", 
    "FreqID",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
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

  /* Ensure inData fully instantiated and OK and selector set */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  return inData;
} /* end getInputData */
  
/*----------------------------------------------------------------------- */
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV    *outData = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong     nvis=1000;
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outData;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Build basic output UV data Object */
  outData = ObitUVFromFileInfo ("out", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  
  /* Set buffer size */
  nvis = 1000;
  ObitInfoListAlwaysPut (outData->info, "nVisPIO",  OBIT_long, dim,  &nvis);
 
 /* Clone from input */
  ObitUVClone (inData, outData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);  

  /* Ensure outData fully instantiated and OK and selector set */
  ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  ObitUVClose (outData, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outData);
  
  ObitErrLog(err); /* Show messages */
  return outData;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Write History for DPCal                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void DPCalHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "BIF", "EIF", "Sources",  "Qual", "FreqID", "souCode", "calSour",
    "doCalSelect", "doCalib", "gainUse", "doPol", "PDVer", "flagVer", 
    "doBand", "BPVer", "Smooth", 
    "outFile",  "outDisk",  "outName", "outClass", "outSeq",
    "solInt", "RLPhase", "PPol", "noIFs",
    NULL};
  gchar *routine = "DPCalHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outData, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
  } else { /* simply copy history */
    ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end DPCalHistory  */

/* Get calibration solutions */
/*----------------------------------------------------------------------- */
/*  Get calibration solutions                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to get data from                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void GetCalSoln (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  olong nchan, i, first, middle, last, iTime=1;
  olong ivis, nvis, ifreq, nfreq, iif, nstok, nif, indx;
  ObitUVDesc *desc;
  ofloat RLPhase[]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
  ofloat RM[]     ={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
  ofloat PPol[]   ={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; 
  ofloat fblank = ObitMagicF();
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat *vis, *ifvis, *fvis, *svis, *rec=NULL;
  ofloat *Isum=NULL, *Iwt=NULL, *Qsum=NULL, *Usum=NULL, *Qwt=NULL, *Uwt=NULL;
  ofloat lastTime=-1.0e3, solInt = 1.0, pdif, pflux;
  odouble freq, reffreq, lambda, reflambda;
  gchar *Stokes="IQUV";
#define MAXTIME 1000   /* Maximum number of times in calData */
  gchar *dataParms[] = {"calSour", NULL};
  gchar *newParms[]  = {"Sources", NULL};
  gchar *routine = "GetCalSoln";

  /* error checks */
  g_assert(ObitErrIsA(err));
  g_assert(ObitUVIsA(inData));
  if (err->error) return;

  /* Get descriptor */
  desc  = inData->myDesc;

  /* Reference wavelength */
  reffreq   = desc->freq;
  reflambda = VELIGHT/reffreq;

  /* Number of channels/poln */
  nchan = desc->inaxes[desc->jlocf];
  if (desc->jlocif>=0)
    nIF   = desc->inaxes[desc->jlocif];
  else
    nIF = 1; 

  /* Allocate calData (global) array [Qfract,Ufract][IF][time] */
  calData  = g_malloc0(2*nIF*MAXTIME*sizeof(ofloat));
  timeData = g_malloc0((MAXTIME+1)*sizeof(ofloat));

  /* Calibrator poln */
  dim[0] = 1;
  ObitInfoListGetTest (myInput, "RM",      &type, dim, RM); 
  ObitInfoListGetTest (myInput, "PPol", &type, dim, PPol); 
  ObitInfoListGetTest (myInput, "RLPhase", &type, dim, RLPhase); 
  for (i=0; i<dim[0]; i++) RLPhase[i] *= DG2RAD;    /* to radians */
  ObitInfoListGetTest (myInput, "solInt", &type, dim, &solInt); 
  if (solInt<=0.0) solInt = 1.0;   /* 1 min */
  solInt /= 1440.0;                /* To days */

  /* Accumulator arrays */
  Isum = g_malloc0(nIF*sizeof(ofloat));
  Iwt  = g_malloc0(nIF*sizeof(ofloat));
  Qsum = g_malloc0(nIF*sizeof(ofloat));
  Qwt  = g_malloc0(nIF*sizeof(ofloat));
  Usum = g_malloc0(nIF*sizeof(ofloat));
  Uwt  = g_malloc0(nIF*sizeof(ofloat));

  /* Want to convert to IQUV */
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut(inData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Select Calibrators */
  ObitInfoListCopyListRename (myInput, inData->info, dataParms, newParms);

  /* Open input */
  retCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if (err->error) goto cleanup;
  desc  = inData->myDesc;
  
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  nstok = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  
  /* loop over blocks of data */
  while (retCode == OBIT_IO_OK) {
    
    /* read buffer */
    retCode = ObitUVReadSelect (inData, NULL, err);
    if (retCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) goto cleanup;
    
    /* how much data? */
    nvis  = desc->numVisBuff;
    
    /* Loop over visibilities in buffer */
    vis = inData->buffer+desc->nrparm;
    rec = inData->buffer;
    /* Trap first time */
    if (lastTime<-1.0e2) {
      lastTime    = rec[desc->iloct];
      timeData[0] = lastTime;
    }
    for (ivis=0; ivis<nvis; ivis++) {
      
      /* Is the previous solution ended? */
      if (rec[desc->iloct]>=(lastTime+solInt)) {
	  timeData[iTime] = lastTime+solInt;   /* End time */
	  lastTime = rec[desc->iloct];
	  /* Get/save averages */
	  indx = (iTime-1)*nIF*2;
	  for (iif=0; iif<nIF; iif++) {
	    if (Iwt[iif]>0.0) Isum[iif] /= Iwt[iif];
	    else              Isum[iif]  = fblank;
	    if (Qwt[iif]>0.0) Qsum[iif] /= Qwt[iif];
	    else              Qsum[iif]  = fblank;
	    if (Uwt[iif]>0.0) Usum[iif] /= Uwt[iif];
	    else              Usum[iif]  = fblank;
	    /* Get model */
	    if ((Isum[iif]!=fblank) && (Qsum[iif]!=fblank) && (Usum[iif]!=fblank)) {
	      pflux = Isum[iif]*PPol[0];  /* Model polarized flux density */
	      /* Model phase difference 
		 IF wavelength */
	      freq   = reffreq + desc->freqIF[iif] + (nfreq/2) * desc->cdelt[desc->jlocf];
	      lambda = VELIGHT/freq;
	      pdif = RLPhase[0] + 2.0 * (lambda*lambda - reflambda*reflambda) * RM[0];

	      /*PARALLACTIC ANGLE CORRECTION?; **************************????????????; */

	      Qsum[iif] -= pflux * cos(pdif);
	      Usum[iif] -= pflux * sin(pdif);
	      Qsum[iif] /= Isum[iif];
	      Usum[iif] /= Isum[iif];
	    } else {
	      Qsum[iif] = Usum[iif] = fblank;
	    }
	    calData[indx++] = Qsum[iif];
	    calData[indx++] = Usum[iif];
	    Qsum[iif] = Qwt[iif] = Usum[iif] = Uwt[iif] = Isum[iif] = Iwt[iif] = 0.0;
	  }
	  iTime++;
	  /* Bounds check */
	  Obit_return_if_fail ((iTime<MAXTIME), err,
			       "%s: Too many time segments, > %d", 
			       routine, MAXTIME);
      }
      /* loop over IFs */
      ifvis = vis;
      for (iif=0; iif<nIF; iif++) {
	
	/* loop over frequency */
	fvis = ifvis;
	for (ifreq = 0; ifreq<nfreq; ifreq++) {
	  
	  /*  Stokes */
	  svis = fvis;
	  /* I */
	  if (*(svis+2)>0.0) {
	    Isum[iif] += *svis * *(svis+2);
	    Iwt[iif]  += *(svis+2);
	  }
	  svis += desc->incs; /* visibility pointer */
	  /* Q */
	  if (*(svis+2)>0.0) {
	    Qsum[iif] += *svis * *(svis+2);
	    Qwt[iif]  += *(svis+2);
	  }
	  svis += desc->incs; 
	  /* U */
	  if (*(svis+2)>0.0) {
	    Usum[iif] += *svis * *(svis+2);
	    Uwt[iif]  += *(svis+2);
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
  
  /* Final solution*/
  timeData[iTime] = rec[desc->iloct];   /* End time */
  /* Get/save averages */
  indx = (iTime-1)*nIF*2;
  for (iif=0; iif<nIF; iif++) {
    if (Iwt[iif]>0.0) Isum[iif] /= Iwt[iif];
    else              Isum[iif]  = fblank;
    if (Qwt[iif]>0.0) Qsum[iif] /= Qwt[iif];
    else              Qsum[iif]  = fblank;
    if (Uwt[iif]>0.0) Usum[iif] /= Uwt[iif];
    else              Usum[iif]  = fblank;
    /* Get model */
    if (Isum[iif]!=fblank) {
      pflux = Isum[iif]*PPol[0];  /* Model polarized flux density */
      /* Model phase difference 
	 IF wavelength */
      freq   = reffreq + desc->freqIF[iif] + (nfreq/2) * desc->cdelt[desc->jlocf];
      lambda = VELIGHT/freq;
      pdif = RLPhase[0] + 2.0 * (lambda*lambda - reflambda*reflambda) * RM[0];
      Qsum[iif] -= pflux * cos(pdif);
      Usum[iif] -= pflux * sin(pdif);
      Qsum[iif] /= Isum[iif];
      Usum[iif] /= Isum[iif];
    } else {
      Qsum[iif] = Usum[iif] = fblank;
    }
    calData[indx++] = Qsum[iif];
    calData[indx++] = Usum[iif];
  }
  nTime = iTime-1;   /* Number of times */

  /* Diagnostics? first. middle, last IFs */
  if (err->prtLv>=2) {
    first = 0; last = nIF-1; middle = (last-first-1)/2;
    Obit_log_error(err, OBIT_InfoErr,
		   "time       Q  IF(%2d) U     Q  IF(%d) U     Q  IF(%d) U     Q    Avgd  U",
		   first+1, middle+1, last+1);
    for (iTime=0; iTime<nTime; iTime++) {
      indx = iTime*2*nIF;
      /* Average */
      Qsum[0] = Usum[0] = Qwt[0] = Uwt[0] = 0.0;
      for (iif=0; iif<nIF; iif++) {
	if (calData[indx+iif*2]!=fblank) {
	  Qsum[0] += calData[indx+iif*2]; Qwt[0] += 1.0;
	}
	if (calData[indx+iif*2+1]!=fblank) {
	  Usum[0] += calData[indx+iif*2+1]; Uwt[0] += 1.0;
	}
      }
      if (Qwt[0]>0.0) Qsum[0] /= Qwt[0];
      else            Qsum[0] = -1.0;
      if (Uwt[0]>0.0) Usum[0] /= Uwt[1];
      else            Usum[0] = -1.0;
      Obit_log_error(err, OBIT_InfoErr,
		     "%8.6f  %6.3f %6.3f   %6.3f %6.3f   %6.3f %6.3f   %7.4f %7.4f",
		     timeData[iTime], calData[indx], calData[indx+1],
		     calData[indx+middle*2], calData[indx+middle*2+1],
		     calData[indx+last*2], calData[indx+last*2+1],
		     Qsum[0], Usum[0]);
    }
    ObitErrLog(err); /* show  messages on err */
  }

  /* Close data */
  retCode = ObitUVClose (inData, err);
  if (err->error) goto cleanup;
 cleanup:
  if (Isum) g_free(Isum);
  if (Qsum) g_free(Qsum);
  if (Usum) g_free(Usum);
  if (Iwt)  g_free(Iwt);
  if (Qwt)  g_free(Qwt);
  if (Uwt)  g_free(Uwt);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
} /* end GetCalSoln */

/*----------------------------------------------------------------------- */
/* Apply calibrator solutions                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to calibrate                                     */
/*      outData   ObitUV to write                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void ApplyCalSoln (ObitInfoList *myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr *err)
{
  ObitIOCode iretCode = OBIT_IO_OK, oretCode = OBIT_IO_OK;
  olong iTime=1, ivis, nvis, ifreq, nfreq, iif, nstok, indx;
  ObitUVDesc *inDesc, *outDesc;
  ofloat fblank = ObitMagicF();
  ObitInfoType type; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat *vis, *ifvis, *fvis, *svis, *rec, Isum, Iwt, IPol, RLPol, LRPol, wt;
  gboolean noIFs=FALSE;
  gchar *Stokes="    ", *today=NULL;
  gchar *dataParms[] = {"Sources", NULL};
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PD","AIPS BP",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  gchar *routine = "ApplyCalSoln";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inData));
  g_assert (ObitUVIsA(outData));

  /* Get descriptors */
  inDesc  = inData->myDesc;
  outDesc = outData->myDesc;

   /* Reference wavelength
  reffreq   = inDesc->freq;
  reflambda = VELIGHT/reffreq; */

  /* Set Stokes */
  dim[0] = strlen(Stokes);
  ObitInfoListAlwaysPut(inData->info, "Stokes", OBIT_string, dim, Stokes);

  /* Select Sources */
  ObitInfoListCopyList (myInput, inData->info, dataParms);

  /* Open input */
  iretCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if (err->error) goto cleanup;

  /* copy Descriptor */
  outData->myDesc = ObitUVDescCopy(inData->myDesc, outData->myDesc, err);

  /* drop IFs in output header? */
  ObitInfoListGetTest (myInput, "noIFs", &type, dim, &noIFs); 
  if (noIFs) {
    outData->myDesc->inaxes[outData->myDesc->jlocf] *= 
      outData->myDesc->inaxes[outData->myDesc->jlocif];
    outData->myDesc->inaxes[outData->myDesc->jlocif] = 1;
  }

  /* Creation date today */
  today = ObitToday();
  strncpy (outData->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outData->buffer) ObitIOFreeBuffer(outData->buffer); /* free existing */
  outData->buffer = NULL;
  outData->bufferSize = -1;

  /* test open output */
  oretCode = ObitUVOpen (outData, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outData->buffer = NULL;
    outData->bufferSize = 0;
    Obit_traceback_msg (err, routine, outData->name);
  }

  /* iretCode = ObitUVClose (inData, err); DEBUG */
  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inData, outData, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inData->mySel->numberSourcesList>1) || (!inData->mySel->selectSources))
  iretCode = ObitUVCopyTables (inData, outData, NULL, sourceInclude, err);
  if (err->error) {
    outData->buffer = NULL;
    outData->bufferSize = 0;
    Obit_traceback_msg (err, routine, inData->name);
  }

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inData->myIO,  inData->info, err);
  oretCode = ObitIOSet (outData->myIO, outData->info, err);
  if (err->error) Obit_traceback_msg (err, routine,inData->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inData, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inData->name);

  iretCode = ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inData->name);
  outData->bufferSize = -1;

  /* Open output */
  oretCode = ObitUVOpen (outData, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  
  nfreq = inDesc->inaxes[inDesc->jlocf];
  nIF = 1;
  if (inDesc->jlocif>=0) nIF = inDesc->inaxes[inDesc->jlocif];
  nstok = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  
  /* loop over blocks of data */
  iTime = 1;
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    
    /* read buffer */
    iretCode = ObitUVReadSelect (inData, NULL, err);
    if (iretCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) goto cleanup;
    
    /* how much data? */
    nvis  = inDesc->numVisBuff;
    outDesc->numVisBuff = inDesc->numVisBuff;
    
    /* Loop over visibilities in buffer */
    vis = inData->buffer+inDesc->nrparm;
    rec = inData->buffer;
    for (ivis=0; ivis<nvis; ivis++) {
      
      /* Is the previous solution ended? */
      if ((rec[inDesc->iloct]>=timeData[iTime]) && (iTime<nTime)) iTime++;

      /* loop over IFs */
      ifvis = vis;
      for (iif=0; iif<nIF; iif++) {
	
	/* loop over frequency */
	fvis = ifvis;
	Isum = Iwt = 0.0;
	for (ifreq = 0; ifreq<nfreq; ifreq++) {
	  /* Average IF Stokes I */
	  svis = fvis;
	  /* I */
	  IPol = wt = 0.0;
	  if (*(svis+2)>0.0) {wt = *(svis+2);  IPol = (*svis)*wt;}           /* RR */
	  svis += inDesc->incs; /* visibility pointer */
	  if (*(svis+2)>0.0) {wt += *(svis+2); IPol += (*svis)*(*(svis+2));} /* LL */
	  if (wt>0.0) {
	    Isum += IPol;   /* Ipol times weight */
	    Iwt  += wt;
	  }
	  fvis += inDesc->incf; /* visibility pointer */
	} /* end average loop over frequencies */

	/* Apply correction */
	/*PARALLACTIC ANGLE CORRECTION; **************************????????????; */
	indx = (iTime-1)*nIF*2 + iif*2;
	if ((Iwt>0.0) && (calData[indx]!=fblank)) {
	  IPol = Isum/Iwt;
	  RLPol = IPol * (calData[indx] - calData[indx+1]);
	  LRPol = IPol * (calData[indx] + calData[indx+1]);
	} else IPol = fblank;
	fvis = ifvis;
	for (ifreq = 0; ifreq<nfreq; ifreq++) {
	  /* IF average IPol and if good cal */
	  svis = fvis + 2*inDesc->incs;
	  if (IPol!=fblank) {
	    *(svis)   -= IPol * calData[indx];  /* real part of RL */
	    *(svis+1) -= IPol * calData[indx+1];    /* imag part of RL */
	    svis += inDesc->incs;
	    *(svis)   -= IPol * calData[indx];  /* real part of LR */
	    *(svis+1) += IPol * calData[indx+1];    /* imag part of LR */
	  } else {
	    /* Bad IPol or cal, flag RL,LR */
	    *(svis+2) = 0.0;
	    svis += inDesc->incs;
	    *(svis+2) = 0.0;
	  }
	  fvis += inDesc->incf; /* visibility pointer */
	} /* end application loop over frequency */
	ifvis += inDesc->incif; /* visibility pointer */
      } /* Loop over IFs */
      
	/* update data pointers */
      vis += inDesc->lrec;
      rec += inDesc->lrec;
    } /* end loop over buffer */
    
    /* Write */
    oretCode = ObitUVWrite (outData, inData->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop over data */
  /* Close data */
 
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outData->buffer = NULL;
  outData->bufferSize = 0;
  
 cleanup:
  iretCode = ObitUVClose (inData, err);
  oretCode = ObitUVClose (outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
} /* end ApplyCalSoln */
