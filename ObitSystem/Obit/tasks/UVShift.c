/* Need: rotate U,V,W to new position */
/* $Id$  */
/*  Shift the phase tracking center of a UV data.    .                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2023                                               */
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

#include "ObitUVDesc.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUV.h"
#include "ObitUVUtil.h"
#include "ObitTableAN.h"
#include "ObitAIPSDir.h"
#include "ObitPosLabelUtil.h"
#include <math.h>
void sincos(double x, double *sin, double *cos); /* Fooey */

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVShiftIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void UVShiftOut (ObitInfoList* outList, ObitErr *err);
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
/* Write history */
void UVShiftHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err);
/* Average data */
ObitUV* ObitUVShift (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/* Program globals */
gchar *pgmName = "UVShift";       /* Program name */
gchar *infile  = "UVShift.in" ;   /* File with program inputs */
gchar *outfile = "UVShift.out";   /* File to contain program outputs */
olong  pgmNumber;               /* Program number (like POPS no.) */
olong  AIPSuser;                /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;                 /* Number of AIPS directories */
gchar **AIPSdirs=NULL;          /* List of AIPS data directories */
olong  nFITS=0;                 /* Number of FITS directories */
gchar **FITSdirs=NULL;          /* List of FITS data directories */
ObitInfoList *myInput  = NULL;  /* Input parameter list */
ObitInfoList *myOutput = NULL;  /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Average data on multiple baselines                                   */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line, read inputs */
  err = newObitErr();

  myInput = UVShiftIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Digest inputs */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get output uvdata */
  outData = setOutputData (myInput, inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Copy Averaging */
  ObitUVShift(inData, outData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Make sure there is no SU table on output */
  ObitUVZapTable (outData, "AIPS SU", -1, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Make sure there is no PD table on output */
  ObitUVZapTable (outData, "AIPS PD", -1, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* History */
  UVShiftHistory (myInput, inData, outData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outData   = ObitUnref(outData);
 
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* UVShiftIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "UVShiftIn";

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

    } else if (strcmp(arg, "-BChan") == 0) { /* BChan */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "BChan", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-EChan") == 0) { /* EChan */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "EChan", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-BIF") == 0) { /* BIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "BIF", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-EIF") == 0) { /* EIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "EIF", OBIT_oint, dim, &itemp, err);
      
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
} /* end UVShiftIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVShift -input file -output ofile [args]\n");
    fprintf(stderr, "UVShift Average multiple baselines\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVShift.in\n");
    fprintf(stderr, "  -output output result file, def UVShift.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BChan first channel to copy\n");
    fprintf(stderr, "  -EChan highest channel to copy\n");
    fprintf(stderr, "  -BIF first IF to copy\n");
    fprintf(stderr, "  -EIF highest IF to copy\n");
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
/*     BChan     Int [1]    channel number, 0=>all                        */
/*     EChan     Int [1]    channel number, 0=>all                        */
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
  ObitInfoListPut (out, "outDType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "UVShift.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "UVShiftName";
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

  /* channels  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BChan", OBIT_oint, dim, &itemp, err);
  itemp = 0; 
  ObitInfoListPut (out, "EChan", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* IFs  */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BIF", OBIT_oint, dim, &itemp, err);
  itemp = 0; 
  ObitInfoListPut (out, "EIF", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS UV file name */
  strTemp = "UVShiftOut.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

 /* correlation type */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListAlwaysPut (out, "corrType", OBIT_oint, dim, &itemp);
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

  /* do3D, Always True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "do3D", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

   /* Stokes always "   " */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Stokes", OBIT_string, dim, strTemp, err);
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
  /*gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    gchar *strTemp;
    gchar *routine = "defaultOutputs";*/
  
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[8], nameStr[13], inName[13];
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Output AIPS UV file class */
  strncpy(strTemp, "     ", 6);
  ObitInfoListGetTest(myInput, "outClass", &type, dim, strTemp);
  if (!strncmp(strTemp, "     ", 5)) {
      strncpy(strTemp, "UVShif", 7);
      dim[0] = strlen (strTemp); dim[1] = 1;
      ObitInfoListAlwaysPut (myInput, "outClass", OBIT_string, dim, strTemp);
  }

  /* Output name, default inName */
  ObitInfoListGetTest(myInput, "outName", &type, dim, nameStr);
  if (!strncmp(nameStr, "     ", 5)) {
    ObitInfoListGetTest(myInput, "inName", &type, dim, inName);
    dim[0] = 12; dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, inName);
  }

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

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
  gboolean     doCalSelect;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "Stokes", "timeRange", "UVRange",
    "BChan", "EChan", "chanInc", "BIF", "EIF", "IFInc", "FreqID", "corrType", 
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "keepLin", "Smooth", "Antennas",  "subA", "Sources", "souCode", "Qual",
    "Shift",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input UV data Object */
  inData = ObitUVFromFileInfo ("in", myInput, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
 
  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
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
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      nvis;
  gboolean  btemp;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     *outParms[] = {  /* Parameters for output data */
    "Compress", NULL};
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /*  Output file not expected to exist */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListAlwaysPut (myInput, "outExist", OBIT_bool, dim, &btemp);

  /* Defaults set in digestInputs */
  /* Create basic output UV Object */
  outUV = ObitUVFromFileInfo ("out", myInput, err);
    
  /* Set buffer size */
  nvis = 1000; type = OBIT_long;
  ObitInfoListGetTest(inData->info, "nVisPIO", &type, dim, &nvis);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO",  type, dim,  &nvis);
    
  /* Clone from input */
  ObitUVClone (inData, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
  
  /* Ensure outUV fully instantiated and OK */
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  /* Get input parameters from myInput, copy to outUV */
  ObitInfoListCopyList (myInput, outUV->info, outParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  ObitErrLog(err); /* Show messages */
  return outUV;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/*  Write History for UVShift                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVShiftHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile", "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "BChan", "EChan", "BIF", "EIF", 
    "Sources", "Qual", "souCode", "subA", "Antennas", 
    "doCalSelect", "doCalib", "gainUse", "doPol", "PDVer", "keepLin", 
    "flagVer", "doBand", "BPVer", "Smooth", "Shift",
    "outFile", "outDisk", "outName", "outClass", "outSeq", "Compress",
    NULL};
  gchar *routine = "UVShiftHistory";

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
 
} /* end UVShiftHistory  */

/**
 * Average baselines, output labeled baseline 1-2
 * \param inUV     Input uv data to copy. 
 *  Control parameter on info
 * \li Shift   OBIT_float (2) shift in asec (NYI)
 * \param err      Error stack, returns if not empty.
 * \return the zeroed ObitUV.
 */
ObitUV* ObitUVShift (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect, doShift=FALSE, uvwneed=FALSE;
  olong i, ii, j, indx, jndx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  olong ilocu, ilocv, ilocw;
  ofloat Shift[]={0.0,0.0}, dShift[]={0.0,0.0}, URot3D[3][3], PRot3D[3][3];
  ofloat visRe, visIm, wt, dxyzc[3]={0.,0.,0.};
  odouble ra, dec, ra0, dec0, phase, cp, sp, u, v, w, uu, vv, ww;
  gchar *rach="RA  ", rast[15], *decch="DEC ", decst[21];
  gchar *routine = "ObitUVShift";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));

  /* Position shift - Full 3D */
  ObitInfoListGetTest(inUV->info, "Shift", &type, dim, Shift);

  /* Clone from input */
  ObitUVClone (inUV, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer = NULL;
  outUV->bufferSize = -1;

  /* Open Files  */
  iretCode = ObitUVOpen (inUV, access, err);
  ObitUVSelSetDesc(inUV->myDesc, inUV->mySel, outUV->myDesc,err);
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  outUV->buffer = inUV->buffer;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;
  ilocu = inDesc->ilocu; ilocv = inDesc->ilocv; ilocw = inDesc->ilocw;
   
  /* Shift */
  doShift = ((Shift[0]!=0.0) || (Shift[1]!=0.0));
  if (doShift) {
    dShift[0] = Shift[0]/3600.; dShift[1] = Shift[1]/3600.; /* In degrees */
    ObitUVDescShiftPosn (inDesc, dShift[0], dShift[1], dxyzc, err);
    /* Shifted position */
    if (err->error) goto cleanup;
    ra0  = inDesc->crval[inDesc->jlocr];
    dec0 = inDesc->crval[inDesc->jlocd];
    ObitSkyGeomXYShift (ra0, dec0, dShift[0], dShift[1], ObitUVDescRotate(inDesc), &ra, &dec);

    /* Set shifted position in output */
    outDesc->crval[outDesc->jlocr] = ra;
    outDesc->crval[outDesc->jlocd] = dec;

    /* uvw Rotation matrix */
    uvwneed = ObitUVDescShift3DPos (inDesc, dShift, 0.0, TRUE, URot3D, PRot3D);

    /* Tell where the shift is to */
    ObitPosLabelUtilRA2HMS (ra, rach, rast);
    ObitPosLabelUtilDec2DMS(dec, decch, decst);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Shifted to %s = %s, %s = %s", rach, rast, decch, decst);
    ObitErrLog(err); 
  } /* end do Shift */
 
  /* we're in business, copy, zero data, set weight to 1 */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* Shift buffer */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec;
      /* U,V,W */
      u = inUV->buffer[indx+ilocu]; v = inUV->buffer[indx+ilocv]; w = inUV->buffer[indx+ilocw];
      jndx = indx + inDesc->nrparm;
      for (j=0; j<inDesc->ncorr; j++) { /* loop over correlations */
	wt = inUV->buffer[jndx+2];
	if (wt>0.0) {
	  /* Shift? */
	  if (doShift) {
	    ii = j/inDesc->inaxes[inDesc->jlocs];
	    phase =  -(+dxyzc[0]*u + dxyzc[1]*v + dxyzc[2]*w)*inDesc->fscale[ii];  /* Scale to freq */
	    sincos(phase, &sp, &cp);
	    visRe = (ofloat)(inUV->buffer[jndx]*cp - inUV->buffer[jndx+1]*sp);
	    visIm = (ofloat)(inUV->buffer[jndx]*sp + inUV->buffer[jndx+1]*cp);
	  } else {
	    visRe = inUV->buffer[jndx]; visIm = inUV->buffer[jndx+1];
	  }
	  inUV->buffer[jndx] = visRe;  inUV->buffer[jndx+1] = visIm;  /* Replace with shifted */
	} /* end if valid */
	jndx += 3;
      } /* end loop over correlations */
      if (uvwneed) {
	uu = u*URot3D[0][0] + v*URot3D[0][1] + w*URot3D[0][2];
	vv = u*URot3D[1][0] + v*URot3D[1][1] + w*URot3D[1][2];
	ww = u*URot3D[2][0] + v*URot3D[2][1] + w*URot3D[2][2];
	inUV->buffer[indx+ilocu] = uu; inUV->buffer[indx+ilocv]= vv; inUV->buffer[indx+ilocw] = ww;
      } /* end rotate u,v,w */
    } /* end loop over visibilities */
    /* Write */
    outDesc->numVisBuff = inDesc->numVisBuff;
    oretCode = ObitUVWrite (outUV, inUV->buffer , err);
    if (err->error) goto cleanup;
  } /* end loop processing data */
  
  
   /* Give report */  
  Obit_log_error(err, OBIT_InfoErr, 
		 "Wrote %d shifted visibilities", outDesc->nvis);
  ObitErrLog(err); 
 /* unset input buffer (may be multiply deallocated ;'{ ) */
 cleanup:
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  
   /* close files */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, inUV->name, outUV);
  
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);

 return outUV;
} /* end ObitUVShift */
