/* $Id$  */
/* Obit Task to append/concatenate uv data           .                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUV.h"
#include "ObitAIPSDir.h"
#include "ObitFITS.h"
#include "ObitTableAN.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* UVAppendIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void UVAppendOut (ObitInfoList* outList, ObitErr *err);
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
void UVAppendHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err);
/* Copy AN table */
olong CopyAN (ObitUV* inData, ObitUV* outData, ObitErr* err);
/* Append data */
void ObitUVAppend(ObitUV *inUV, ObitUV *outUV, ObitErr *err);
/* Copy AN Tables */
ObitIOCode ObitTableANSelect2 (ObitUV *inUV, olong ncopy, ObitUV *outUV, 
			       ObitErr *err);

/* Program globals */
gchar *pgmName = "UVAppend";       /* Program name */
gchar *infile  = "UVAppend.in" ;   /* File with program inputs */
gchar *outfile = "UVAppend.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
gboolean newOutput=TRUE;       /* Did output previous exist? */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit Task to append the data from one dataset to a compatable one    */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  olong        newSubA;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitErr      *err= NULL;

   /* Startup - parse command line, read inputs */
  err = newObitErr();

  myInput = UVAppendIn (argc, argv, err);
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

  /* Copy input AN table, get new subarray number */
  newSubA = CopyAN (inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Set new subarray number */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(inData->info, "newSubA", OBIT_long, dim, &newSubA);

  /* Append */
  ObitUVAppend(inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* History */
  UVAppendHistory (myInput, inData, outData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
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

ObitInfoList* UVAppendIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "UVAppendIn";

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
} /* end UVAppendIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: UVAppend -input file -output ofile [args]\n");
    fprintf(stderr, "UVAppend Obit task to subtract model from UV data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def UVAppend.in\n");
    fprintf(stderr, "  -output output result file, def UVAppend.out\n");
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
  strTemp = "UVAppend.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "UVAppendName";
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
  strTemp = "UVAppendOut.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file name */
  strTemp = "UVAppendOut";
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
  /*  ObitInfoType type;
      gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
      gchar *strTemp, *opcode=NULL; */
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

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
  ObitInfoType type;
  olong         Aseq, disk, cno, nvis=1000;
  gchar        *Type, *strTemp, inFile[129];
  oint         doCalib;
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Sources", "souCode", "Qual", "Stokes", "timeRange", 
    "BChan", "EChan", "chanInc", "BIF", "EIF", "IFInc", "FreqID", "corrType", 
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "Smooth", "Antennas",  "subA", "Sources", "souCode", "Qual",
     NULL};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
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
    ObitUVSetFITS (inData, nvis, disk, inFile,  err); 
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inData;
  }

  /* Make sure doCalSelect set properly */
  doCalSelect = TRUE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  /* Always */
  doCalSelect = TRUE;
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
/*  Create output uv data                                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   In global                                                            */
/*      newOutput TRUE if output just created                             */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputData (ObitInfoList *myInput, ObitUV* inData, ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  olong      i, n, Aseq, disk, cno;
  gchar     *Type, *strTemp, outFile[129];
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129];
  gchar     *outFull=NULL, *FITS = "FITS";
  gchar     *outParms[] = {  /* Parameters for output data */
    "Compress", NULL};
  gchar     *routine = "setOutputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* Output File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */

    /* outName given? */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
    for (i=0; i<12; i++) Aname[i] = ' ';  Aname[i] = 0;
    for (i=0; i<MIN(12,dim[0]); i++) Aname[i] = strTemp[i];
    /* Save any defaulting on myInput */
    dim[0] = 12;
    ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, Aname);

      
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    /* Default out class is "UVCop" */
    if (!strncmp(Aclass, "      ", 6)) strncpy (Aclass, "UVCop", 7);

    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    newOutput = !exist;  /* Is this new? */
    
    /* define object */
    nvis = 1000;
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS UV data %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */

    /* outFile given? */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&strTemp);
    /* if not use inName */
    if ((strTemp==NULL) || (!strncmp(strTemp, "            ", 12)))
      ObitInfoListGetP (myInput, "inFile", &type, dim, (gpointer)&strTemp);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) outFile[i] = strTemp[i]; outFile[i] = 0;
    ObitTrimTrail(outFile);  /* remove trailing blanks */

    /* Save any defaulting on myInput */
    dim[0] = strlen(outFile);
    ObitInfoListAlwaysPut (myInput, "outFile", OBIT_string, dim, outFile);

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

    /* Does output exist? */
    outFull = ObitFITSFilename (disk, outFile, err);  /* Full path */
    /* Does new filename exist? */
    newOutput = !(ObitFileExist (outFull, err) || err->error);
    if (outFull) g_free(outFull);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
   
    /* define object */
    nvis = 1000;
    ObitUVSetFITS (outUV, nvis, disk, outFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS UV data %s on disk %d", outFile, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  /* Get input parameters from myInput, copy to outUV */
  ObitInfoListCopyList (myInput, outUV->info, outParms);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  /* Clone from input if new */
  if (newOutput) {
    ObitUVClone (inData, outUV, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* Ensure outUV fully instantiated and OK */
    ObitUVFullInstantiate (outUV, FALSE, err);
    Obit_log_error(err, OBIT_InfoErr, "Create new output dataset");
  } else {
    ObitUVFullInstantiate (outUV, TRUE, err);
    Obit_log_error(err, OBIT_InfoErr, "Appending to existing dataset");
  }
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  ObitErrLog(err); /* Show messages */
  return outUV;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/*  Write History for UVAppend                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void UVAppendHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "BChan", "EChan", "chanInc", "BIF", "EIF", "IFInc", "Stokes", 
    "Sources",  "Qual", "souCode", "subA", "Antennas", 
    "doCalSelect", "doCalib", "gainUse", "doPol", "PDVer", "flagVer", 
    "doBand", "BPVer", "Smooth",  "corrType", 
    "outFile",  "outDisk",  "outName", "outClass", "outSeq", "Compress",
    NULL};
  gchar *routine = "UVAppendHistory";

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
 
} /* end UVAppendHistory  */

/**
 * Copy AN tables from inData to outData with selection
 * \param inData   Input UV 
 * \param outData  Output UV, must already be defined
 * \param err      Error stack
 * \return the number of the first new subarray
 */
olong CopyAN (ObitUV* inData, ObitUV* outData, ObitErr* err)
{
  olong newSubA=0;
  olong oldANHi, newANHi;
  gchar *routine = "CopyAN";

  /* error checks */
  if (err->error) return newSubA;

  /* Check global newOutput */
  if (newOutput) return newSubA;

  /* Highest input AN */
  newANHi = ObitTableListGetHigh (inData->tableList, "AIPS AN");

  /* Highest extant output AN */
  oldANHi = ObitTableListGetHigh (outData->tableList, "AIPS AN");
  newSubA = oldANHi + 1;   /* New output AN = subarray */
  
  /* Copy AN Tables */
  ObitTableANSelect2 (inData, newANHi, outData, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, newSubA);

  Obit_log_error(err, OBIT_InfoErr, 
		 "Copying %d subarrays starting with %d",newANHi,newSubA);
  ObitErrLog(err); 

  return newSubA;
} /* end CopyAN */

/**
 * Append the contents of one UV onto the end of another
 * \param inUV   Input UV 
 *      Control info on info member in addition to data selection/calibration
 * \li newSubA OBIT_long scalar subarray number of new data [def 0]
 *     This affects the baseline abd time random parameters of the output.
 * \param outUV  Output UV, must already be defined
 * \param err    Error stack
 */
void ObitUVAppend(ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode;
  gboolean doCalSelect, done;
  ObitInfoType type;
  ObitIOAccess access;
  gboolean incompatible;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *inDesc, *outDesc;
  olong inNPIO, outNPIO, NPIO, indx, i, newSubA=0;
  ofloat timeAdd, blAdd;
  gchar *routine = "ObitUVAppend";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* New subarray number */
  ObitInfoListGetTest(inUV->info, "newSubA", &type, dim, &newSubA);
  newSubA = MAX (1, newSubA);
  /* 5 day offset for each subarray */
  timeAdd = 5.0 * (newSubA-1);
  /* 0.01 offset per subarray to baseline */
  blAdd = 0.01 * (newSubA-1);

  /* Get input descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Check compatability between inUV, outUV */
  incompatible = (inDesc->ncorr!=outDesc->ncorr);
  incompatible = incompatible || (inDesc->jlocs!=outDesc->jlocs);
  incompatible = incompatible || (inDesc->jlocf!=outDesc->jlocf);
  incompatible = incompatible || (inDesc->jlocif!=outDesc->jlocif);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible structures",
		   routine);
      return ;
 }

  /* Check position - same equinox and within 1 mas */
  incompatible = fabs(inDesc->equinox-outDesc->equinox) > 0.01;
  incompatible = incompatible || 
    fabs(inDesc->crval[inDesc->jlocr]-outDesc->crval[inDesc->jlocr]) > 0.001/3600;
  incompatible = incompatible || 
    fabs(inDesc->crval[inDesc->jlocd]-outDesc->crval[inDesc->jlocd]) > 0.001/3600;
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible positions",
		   routine);
      return ;
 }

  /* Check frequency - within 1 Hz */
  incompatible = fabs(inDesc->freq-outDesc->freq) > 1.0;
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible frequencies",
		   routine);
      return ;
 }

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Set number of vis per I/O */
  inNPIO = 1000;
  ObitInfoListGetTest (inUV->info, "nVisPIO", &type, dim, &inNPIO);
  outNPIO = 1000;
  ObitInfoListGetTest (outUV->info, "nVisPIO", &type, dim, &outNPIO);
  NPIO = 1000; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (inUV->info,  "nVisPIO", OBIT_long, dim,  &NPIO);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim,  &NPIO);

  /* Open Input Data */
  retCode = ObitUVOpen (inUV, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inUV->name);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer     = inUV->buffer;
  outUV->bufferSize = inUV->bufferSize;

  /* Open Output Data */
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outUV->buffer = NULL; /* remove pointer to inUV buffer */
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }
  outDesc->firstVis = outDesc->nvis+1; /* Write to end */

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

    /* How many? */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify input time and baseline */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec;
      inUV->buffer[indx+inDesc->iloct] += timeAdd;
      inUV->buffer[indx+inDesc->ilocb] += blAdd;
    } /* end visibility loop */

    /* Write buffer */
    retCode = ObitUVWrite (outUV, NULL, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, outUV->name);
    }
  } /* end loop over data */
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* Close input */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  /* Close output */
  retCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Reset number of vis per I/O */
  ObitInfoListAlwaysPut (inUV->info,  "nVisPIO", OBIT_long, dim,  &inNPIO);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim,  &outNPIO);

} /* end ObitUVAppend */

/**
 * Copies AN tables from inUV to outUV with selection in inUV
 * AN tables copies are added at higher version numbers than 
 * any extant tables (old one kept).
 * If poln calibration is selected on inUV, polarization 
 * calibration info is removed.
 * \param inUV     Input UV to copy from
 * \param ncopy    Number of Table to copy starting from version 1
 * \param outUV    Output UV to copy to
 * \param *err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableANSelect2 (ObitUV *inUV, olong ncopy, ObitUV *outUV, 
			       ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableAN    *inTab=NULL, *outTab=NULL;
  ObitTableANRow *inRow=NULL, *outRow=NULL;
  ObitInfoType type;
  olong iif, oif, i, polRefAnt;
  olong iANver, oANver, inANRow, outANRow;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  oint numOrb, numIF, numPCal;
  odouble dtemp;
  gboolean wanted, doPol;
  gchar tempName[MAXKEYCHARTABLEAN+4];
  gchar *CopyList[] = {"ARRAYX", "ARRAYY", "ARRAYZ", "GSTIA0",
		       "DEGPDY", "FREQ",   "RDATE",  "POLARX", "POLARY",
		       "UT1UTC", "DATUTC", "TIMSYS", "ARRNAM",
		       "NUMORB", "NOPCAL", "FREQID", "IATUTC", 
		       "POLTYPE", "P_REFANT", 
		       "P_DIFF01", "P_DIFF02", "P_DIFF03", "P_DIFF04", 
		       "P_DIFF05", "P_DIFF06", "P_DIFF07", "P_DIFF08", 
		       NULL};
  gchar *routine = "ObitTableANSelect2";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Poln calibration selected? */
  doPol = FALSE;
  ObitInfoListGetTest(inUV->info, "doPol", &type, (gint32*)dim, &doPol);

  /* Fully instantiate UV files */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error )Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Open Output Data */
  ObitUVOpen (outUV, OBIT_IO_ReadWrite, err) ;
  if (err->error )Obit_traceback_val (err, routine, outUV->name, retCode);

  /* Loop over AN tables */
  for (iANver=1; iANver<=ncopy; iANver++) {

    /* Get input table */
    numOrb  = 0;
    numPCal = 0;
    numIF   = 0;
    inTab = 
      newObitTableANValue (inUV->name, (ObitData*)inUV, &iANver, OBIT_IO_ReadOnly, 
			   numIF, numOrb, numPCal, err);
    if (err->error) Obit_traceback_val (err, routine, inTab->name, retCode);
    /* Find it */
    if (inTab==NULL) continue;  /* No keep looping */

    /* Open input table */
    retCode = ObitTableANOpen (inTab, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);

    /* Create output table */
    numOrb  = inTab->numOrb;
    numIF   = inUV->mySel->numberIF;
    if (inTab->numPCal>0) numPCal = 2;
    else numPCal = 0;
    oANver = -1;   /* Force new AN table */
    outTab = 
      newObitTableANValue (outUV->name, (ObitData*)outUV, &oANver, OBIT_IO_WriteOnly, 
			    numIF, numOrb, numPCal, err);
    if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);
    /* Create it? */
    Obit_retval_if_fail((outTab!=NULL), err, retCode,
			"%s: Could not create AN table %d for %s", 
			routine, iANver, outTab->name);

   /* Open output table */
  retCode = ObitTableANOpen (outTab, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inTab->name, retCode);
  
   /* Update header info */
    outTab->ArrayX  = inTab->ArrayX;
    outTab->ArrayY  = inTab->ArrayY;
    outTab->ArrayZ  = inTab->ArrayZ;
    outTab->GSTiat0 = inTab->GSTiat0;
    outTab->DegDay  = inTab->DegDay;
    outTab->Freq    = inTab->Freq;
    outTab->PolarX  = inTab->PolarX;
    outTab->PolarY  = inTab->PolarY;
    outTab->ut1Utc  = inTab->ut1Utc;
    outTab->dataUtc = inTab->dataUtc;
    outTab->FreqID  = inTab->FreqID;
    outTab->iatUtc  = inTab->iatUtc;
    outTab->P_Refant  = inTab->P_Refant;
    outTab->P_Diff01  = inTab->P_Diff01;
    outTab->P_Diff02  = inTab->P_Diff02;
    outTab->P_Diff03  = inTab->P_Diff03;
    outTab->P_Diff04  = inTab->P_Diff04;
    outTab->P_Diff05  = inTab->P_Diff05;
    outTab->P_Diff06  = inTab->P_Diff06;
    outTab->P_Diff07  = inTab->P_Diff07;
    outTab->P_Diff08  = inTab->P_Diff08;
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->RefDate[i] = inTab->RefDate[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->TimeSys[i] = inTab->TimeSys[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->ArrName[i] = inTab->ArrName[i];
    for (i=0; i<MAXKEYCHARTABLEAN; i++)  
      outTab->polType[i] = inTab->polType[i];

   /* Copy InfoList stuff */
    ObitInfoListCopyList (inTab->myDesc->info, outTab->myDesc->info, CopyList);

    /* Poln cal info in info member */
    g_snprintf (tempName, 7, "       ");
    dim[0] = 7; dim[1] = dim[2] = dim[3] = 1;
    type = OBIT_string;
    ObitInfoListGetTest(inTab->myDesc->info, "POLTYPE", &type, dim, &tempName);
    tempName[dim[0]] = 0;
    if (doPol) {  /* If cal, blank */
      g_snprintf (tempName, 7, "       ");
    }
    ObitInfoListAlwaysPut(outTab->myDesc->info, "POLTYPE", type, dim, tempName);
    
    /* Polarization reference antenna */
    dim[0] = dim[1] = 1; type = OBIT_long; polRefAnt = 0;
    ObitInfoListGetTest(inTab->myDesc->info, "P_REFANT", &type, dim, &polRefAnt);
    ObitInfoListAlwaysPut(outTab->myDesc->info, "P_REFANT", type, dim, &polRefAnt);

    /* R-L Phase differences */
    if (!doPol) {  
      oif = 1;
      for (iif=inUV->mySel->startIF; 
	   iif<=inUV->mySel->startIF+inUV->mySel->numberIF-1;
	   iif++) {
	g_snprintf (tempName, 9, "P_DIFF%2.2d",iif);
	dim[0] = dim[1] = 1; type = OBIT_double; dtemp = 0.0;
	ObitInfoListGetTest(inTab->myDesc->info, tempName, &type, dim, &dtemp);
	g_snprintf (tempName, 9, "P_DIFF%2.2d",oif);
	ObitInfoListAlwaysPut(outTab->myDesc->info, "tempName", type, dim, &polRefAnt);
	oif++;
      }
    }

    /* Set rows */
    inRow  = newObitTableANRow (inTab);
    outRow = newObitTableANRow (outTab);
    ObitTableANSetRow (outTab, outRow, err);
    if (err->error) Obit_traceback_val (err, routine, outTab->name, retCode);

    /* Loop over table copying selected data */
    outANRow = -1;
 
    for (inANRow=1; inANRow<=inTab->myDesc->nrow; inANRow++) {
      retCode = ObitTableANReadRow (inTab, inANRow, inRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
      if (inRow->status==-1) continue;
  
      /* Want this one? */
      wanted = ObitUVSelWantAnt(inUV->mySel, inRow->noSta);
      if (!wanted) continue;

      /* Copy selected data */
      outRow->noSta   = inRow->noSta;
      outRow->staXof  = inRow->staXof;
      outRow->PolAngA = inRow->PolAngA;
      outRow->PolAngB = inRow->PolAngB;
      for (i=0; i<8; i++) outRow->AntName[i] = inRow->AntName[i];
      for (i=0; i<3; i++) outRow->StaXYZ[i]  = inRow->StaXYZ[i];
      for (i=0; i<numOrb; i++) outRow->OrbParm[i] = inRow->OrbParm[i];
      outRow->polTypeA[0] = inRow->polTypeA[0];
      outRow->polTypeB[0] = inRow->polTypeB[0];
      /* IF dependent poln cal */
      oif = 0; 
      for (iif=inUV->mySel->startIF-1;  
	   iif<inUV->mySel->startIF+inUV->mySel->numberIF-1; 
	   iif++) { 
	if (doPol && (numPCal>0)) { /* zero cal */
	  outRow->PolCalA[2*oif]   = 0.0; 
	  outRow->PolCalA[2*oif+1] = 0.0; 
	  outRow->PolCalB[2*oif]   = 0.0; 
	  outRow->PolCalB[2*oif+1] = 0.0; 
	} else if (numPCal>0) {  /* Copy */
	  outRow->PolCalA[2*oif]   = inRow->PolCalA[2*iif]; 
	  outRow->PolCalA[2*oif+1] = inRow->PolCalA[2*iif+1];
	  outRow->PolCalB[2*oif]   = inRow->PolCalB[2*iif]; 
	  outRow->PolCalB[2*oif+1] = inRow->PolCalB[2*iif+1];
	} 
      	oif++; 
      } 

      retCode = ObitTableANWriteRow (outTab, outANRow, outRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, inUV->name, retCode);
    } /* end loop over rows */
    
    /* Close tables */
    retCode = ObitTableANClose (inTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, inTab->name, retCode);
    retCode = ObitTableANClose (outTab, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, outTab->name, retCode);
 
    /* release table objects */
    inTab  = ObitTableANUnref(inTab);
    outTab = ObitTableANUnref(outTab);

    /* release row objects */
    inRow  = ObitTableANRowUnref(inRow);
    outRow = ObitTableANRowUnref(outRow);
  } /* end loop over tables */

  /* Close output */
  ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, retCode);

  return retCode;
} /* end ObitTableANSelect2 */
