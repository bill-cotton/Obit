/* $Id$  */
/*  Convert linear feed basis to circular .                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2018,2022                                          */
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
#include "ObitTableANUtil.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitAIPSDir.h"
#include "ObitSinCos.h"
#include "ObitPrecess.h"
#include "ObitAntennaList.h"
#include "ObitSourceList.h"
/* internal prototypes */
/* Get inputs */
ObitInfoList* Lin2CirIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void Lin2CirOut (ObitInfoList* outList, ObitErr *err);
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
void Lin2CirHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err);
/* Convert data */
ObitUV* ObitUVLin2Cir (ObitUV *inUV, ObitUV *outUV, ObitErr *err);

/* Program globals */
gchar *pgmName = "Lin2Cir";       /* Program name */
gchar *infile  = "Lin2Cir.in" ;   /* File with program inputs */
gchar *outfile = "Lin2Cir.out";   /* File to contain program outputs */
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
/*  Convert linear feed basis to circular .                               */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem=NULL;
  ObitUV       *inData=NULL, *outData=NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line, read inputs */
  err = newObitErr();

  myInput = Lin2CirIn (argc, argv, err);
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

  /* Convert */
  ObitUVLin2Cir(inData, outData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* History */
  Lin2CirHistory (myInput, inData, outData, err);
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

ObitInfoList* Lin2CirIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "Lin2CirIn";

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
} /* end Lin2CirIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Lin2Cir -input file -output ofile [args]\n");
    fprintf(stderr, "Lin2Cir Convert XX,YY... to RR,LL...\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Lin2Cir.in\n");
    fprintf(stderr, "  -output output result file, def Lin2Cir.out\n");
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
  strTemp = "Lin2Cir.intab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS input uv file name */
  strTemp = "Lin2CirName";
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
  strTemp = "Lin2CirOut.uvtab";
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar strTemp[8], nameStr[13], inName[13];
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Output AIPS UV file class */
  strncpy(strTemp, "Ln2Cr", 6);
  ObitInfoListGetTest(myInput, "outClass", &type, dim, strTemp);
  if (!strncmp(strTemp, "     ", 5)) {
      strncpy(strTemp, "Ln2Cr", 6);
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
    "doPol", "PDVer", "Smooth", "Antennas",  "subA", "Sources", "souCode", "Qual",
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

  /* MUST be Linear with 4 corelations */
  Obit_retval_if_fail ((inData->myDesc->crval[inData->myDesc->jlocs]==-5.0), err,inData,
		       "%s: input MUST be linear basis data", routine);  
  Obit_retval_if_fail ((inData->myDesc->inaxes[inData->myDesc->jlocs]>=4), err,inData,
		       "%s: input MUST have at least 4 Stokes, have  %d",  
		       routine, inData->myDesc->inaxes[inData->myDesc->jlocs]);  

 
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

  /* Change to circular basis (RR, LL,RL,LR) */
  outUV->myDesc->crval[outUV->myDesc->jlocs] = -1.0;
  
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
/*  Write History for Lin2Cir                                             */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   ObitUV to write history to                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void Lin2CirHistory (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
		   ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", 
    "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "BChan", "EChan", "BIF", "EIF", 
    "Sources",  "Qual", "souCode", "subA", "Antennas", 
    "doCalSelect", "doCalib", "gainUse", "doPol", "PDVer", "flagVer", 
    "doBand", "BPVer", "Smooth",  "corrType", 
    "outFile",  "outDisk",  "outName", "outClass", "outSeq", "Compress",
    NULL};
  gchar *routine = "Lin2CirHistory";

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
 
} /* end Lin2CirHistory  */

  /** Private: 4x4 complex matrix * 4x1 complex vector multiply */
  static void MatxVec4Mult(ofloat in1[32], ofloat in2[8], ofloat out[8])
{
  olong ic, ii, n=4;
  ofloat sumr, sumi;
  for (ic=0; ic<4; ic++) {
    sumr = sumi = 0.0;
    for (ii=0; ii<4; ii++) {
      sumr += in1[(ic*n+ii)*2]*in2[ii*2]   - in1[(ic*n+ii)*2+1]*in2[ii*2+1];
      sumi += in1[(ic*n+ii)*2]*in2[ii*2+1] + in1[(ic*n+ii)*2+1]*in2[ii*2];
    }
    out[ic*2]   = sumr;
    out[ic*2+1] = sumi;
  }
} /* end MatxVec4Mult */

/** Private: Muller matrix from outer product of Jones matrices 
   Supports blanking */
static void MatxOuter(ofloat in1[8], ofloat in2[8], ofloat out[32])
{
  /* out = in1 (outer product) conjg(in2) */
  out[0]  =  in1[0] * in2[0] + in1[1] * in2[1];  out[1]  =  in1[1] * in2[0] - in1[0] * in2[1];
  out[2]  =  in1[0] * in2[2] + in1[1] * in2[3];  out[3]  =  in1[1] * in2[2] - in1[0] * in2[3];
  out[4]  =  in1[2] * in2[0] + in1[3] * in2[1];  out[5]  =  in1[3] * in2[0] - in1[2] * in2[1];
  out[6]  =  in1[2] * in2[2] + in1[3] * in2[3];  out[7]  =  in1[3] * in2[2] - in1[2] * in2[3];
  
  out[8]  =  in1[0] * in2[4] + in1[1] * in2[5];  out[9]  =  in1[1] * in2[4] - in1[0] * in2[5];
  out[10] =  in1[0] * in2[6] + in1[1] * in2[7];  out[11] =  in1[1] * in2[6] - in1[0] * in2[7];
  out[12] =  in1[2] * in2[4] + in1[3] * in2[5];  out[13] =  in1[3] * in2[4] - in1[2] * in2[5];
  out[14] =  in1[2] * in2[6] + in1[3] * in2[7];  out[15] =  in1[3] * in2[6] - in1[2] * in2[7];
  
  out[16] =  in1[4] * in2[0] + in1[5] * in2[1];  out[17] =  in1[5] * in2[0] - in1[4] * in2[1];
  out[18] =  in1[4] * in2[2] + in1[5] * in2[3];  out[19] =  in1[5] * in2[2] - in1[4] * in2[3];
  out[20] =  in1[6] * in2[0] + in1[7] * in2[1];  out[21] =  in1[7] * in2[0] - in1[6] * in2[1];
  out[22] =  in1[6] * in2[2] + in1[7] * in2[3];  out[23] =  in1[7] * in2[2] - in1[6] * in2[3];
  
  out[24] =  in1[4] * in2[4] + in1[5] * in2[5];  out[25] =  in1[5] * in2[4] - in1[4] * in2[5];
  out[26] =  in1[4] * in2[6] + in1[5] * in2[7];  out[27] =  in1[5] * in2[6] - in1[4] * in2[7];
  out[28] =  in1[6] * in2[4] + in1[7] * in2[5];  out[29] =  in1[7] * in2[4] - in1[6] * in2[5];
  out[30] =  in1[6] * in2[6] + in1[7] * in2[7];  out[31] =  in1[7] * in2[6] - in1[6] * in2[7];
  } /* end  MatxOuter */

/**
 * Convert linear feed basis data to circular for vis record
 * \param Muller  Inverse Muller matrix
 * \param inData  input array of r, i, w for XX,YY,XY,YX for ncorr vis
 * \param outData output array of r, i, w for XX,YY,XY,YX for ncorr vis
 * \param ncorr   Number of correlations
 */
void ConvertL2C(ofloat Muller[32], ofloat *inData, ofloat *outData, olong ncorr)
{
  olong i, j;
  ofloat *in=inData, *out=outData;
  ofloat iV[8], oV[8],sum;
  gboolean flag;

  /* Loop over correlations */
  for (i=0; i<ncorr-1; i++) {
    /* anything flagged? -> zero all */
    flag = FALSE;  for (j=0; j<12; j+=3) if (in[j+2]<=0.0) {flag = TRUE; break;}
    if (flag) {
      for (j=0; j<12; j++) out[j] = 0.0;
      in += 12; out += 12;  /* update pointers */
      continue;
    }
    /* Load input vector [XX,XY,YX,YY]*/
    sum = 0.0;
    for (j=0; j<4; j++) sum += in[j*3+2];
    iV[0] = in[0]; iV[1] = in[1];
    iV[2] = in[6]; iV[3] = in[7];
    iV[4] = in[9]; iV[5] = in[10];
    iV[6] = in[3]; iV[7] = in[4];

    /* Multiply by inverse Muller matrix */
    MatxVec4Mult(Muller, iV, oV);

    /* unload output vector from [RR,RL,LR,LL] */
    sum *= 0.25;
    out[0] = oV[0]; out[1]  = oV[1]; out[2]  = sum;
    out[3] = oV[6]; out[4]  = oV[7]; out[5]  = sum;
    out[6] = oV[2]; out[7]  = oV[3]; out[8]  = sum;
    out[9] = oV[4]; out[10] = oV[5]; out[11] = sum;
    in += 12; out += 12;  /* update pointers */
  } /* end loop over correlations */
} /* end ConvertL2C */

/* Nominal inverse of Jones matrix for perfect inear feed */
void IJones(ofloat J[8])
{
  ofloat elp_x=0.0, elp_y=0.0, ori_x=0.0, ori_y=G_PI*0.5;
  ofloat angle[4], sina[4], cosa[4], Jones[8] ,Det[2], d;

  angle[0] = G_PI*0.25+elp_x; angle[1] = G_PI*0.25-elp_y;
  angle[2] = ori_x;           angle[3] = ori_y;
  ObitSinCosVec(4, angle, sina, cosa);
  Jones[0] =  cosa[0]*cosa[2]; Jones[1] = -cosa[0]*sina[2];
  Jones[2] =  sina[0]*cosa[2]; Jones[3] =  sina[0]*sina[2];
  Jones[4] =  sina[1]*cosa[3]; Jones[5] = -sina[1]*sina[3];
  Jones[6] =  cosa[1]*cosa[3]; Jones[7] =  cosa[1]*sina[3];
  /* inverse of determinant */
  Det[0] = (Jones[0]*Jones[6] - Jones[1]*Jones[7]) - (Jones[2]*Jones[4] - Jones[3]*Jones[5]);
  Det[1] = (Jones[0]*Jones[7] + Jones[1]*Jones[6]) - (Jones[2]*Jones[5] + Jones[3]*Jones[4]);
  /* Inverse of determinant */
  d = Det[0]*Det[0] + Det[1]*Det[1];
  if (d!=0.0) d = 1.0 / d;
  else d = 1.0;
  Det[0] *=  d;
  Det[1] *= -d;
  /* Inverse matrix out */
  J[6] =   Jones[0] * Det[0] - Jones[1] * Det[1];
  J[7] =   Jones[0] * Det[1] + Jones[1] * Det[0];
  J[2] = -(Jones[2] * Det[0] - Jones[3] * Det[1]);
  J[3] = -(Jones[2] * Det[1] + Jones[3] * Det[0]);
  J[4] = -(Jones[4] * Det[0] - Jones[5] * Det[1]);
  J[5] = -(Jones[4] * Det[1] + Jones[5] * Det[0]);
  J[0] =   Jones[6] * Det[0] - Jones[7] * Det[1];
  J[1] =   Jones[6] * Det[1] + Jones[7] * Det[0];
} /* end IJones */

/**
 *  Get Source List
 * \param   inUV       ObitUV with SU Table
 * \param   err        Obit Error stack
 * \return SourceList, should be Unrefed when done
 */
ObitSourceList* GetSourceList (ObitUV* inUV, ObitErr* err)
{
  olong iver = 1;
  ObitTableSU *SUTable=NULL;
  ObitSourceList  *SList=NULL;
  gchar *routine = "GetSourceList";

  if (err->error) return SList;

  SUTable = newObitTableSUValue ("SUTable", (ObitData*)inUV, &iver, 
			       OBIT_IO_ReadOnly, 0, err);
  if (SUTable) {
    SList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, SUTable->name, SList);
    SUTable = ObitTableSUUnref(SUTable);   /* Done with table */
  } else {  /* Use position from header */
    SList = ObitSourceListCreate ("SList", 1);
    SList->SUlist[0]->equinox = inUV->myDesc->equinox;
    SList->SUlist[0]->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    SList->SUlist[0]->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, SList->SUlist[0]);
  }
  return SList;
} /* end GetSourceList */

/**
 * Convert linear feed basis data to circular for data file
 * \param inUV     Input uv data to convery
 *  Control parameter on info
 * \param outUV    Output uv dat
 * \param err      Error stack, returns if not empty.
 * \return the modified ObitUV.
 */
ObitUV* ObitUVLin2Cir (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  olong i, j, indx, incs, icor, jndx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  ofloat J[8], Muller[32], tr, ti, gr=1.0, gi=0.0, PA, cosPA, sinPA;
  ObitTableAN *ANTable=NULL;
  olong ver, numIF=1, numOrb,  numPCal, ncorr;
  olong suId=-1, lastSu=0;
  ofloat lastTime=-1.0e-20, time;
  ObitAntennaList *AntList=NULL;   /* Antenna list*/
  ObitSourceList *SList=NULL;
  gchar *routine = "ObitUVLin2Cir";
  /* DEBUG 
  ofloat ipol=1.234,qpol=0.01,upol=-0.2,vpol=0.05;*/
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));

  /* Generate inverse Muller Matrix for conversion */
  /* First, nominal inverse Jones */
  IJones (J);
  MatxOuter (J, J, Muller);

  /* Clone output from input */
  ObitUVClone (inUV, outUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open Files  */
  iretCode = ObitUVOpen (inUV, access, err);
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  /* Change to circular basis (RR, LL,RL,LR) */
  outUV->myDesc->crval[outUV->myDesc->jlocs] = -1.0;
  
 
  /* Get descriptors */
  inDesc  = inUV->myDesc;
  incs    = inDesc->incs;
  ncorr   = inDesc->ncorr/inDesc->inaxes[inDesc->jlocs];
  outDesc = outUV->myDesc;

  /* Get antenna list */
  /* Create output Antenna table object */
  ver      = 1;   numOrb   = 0;   numPCal  = 0;
  ANTable = newObitTableANValue ("AN table", (ObitData*)outUV, 
				 &ver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
  if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
  AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);

  /* Get Source List */
  SList = GetSourceList (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);
  if (inDesc->ilocsu<0) {
    /* May have selected one source */
    if (inUV->mySel->selectSources && (inUV->mySel->numberSourcesList==1)) 
         suId = inUV->mySel->sources[0];
    else suId = 1;
  }
 
  /* we're in business, convert */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* Loop over buffer */
    outDesc->numVisBuff = 0 /* no. vis in buffer */;
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over buffer */
      indx = i*inDesc->lrec;  /* Data offset */
      /* Copy random parameters */
      for (j=0; j<inDesc->nrparm; j++) outUV->buffer[indx+j] = inUV->buffer[indx+j];

      /* source/time */
      if (inDesc->ilocsu>0) suId = (olong)(0.5 + inUV->buffer[indx+inDesc->ilocsu]);
      time = inUV->buffer[indx+inDesc->iloct];
      indx += inDesc->nrparm;

      /* New source/time? */
      if ((time>lastTime) || (suId!=lastSu)) {
	/* Parallactic angle  - assume EVLA */
	PA = ObitAntennaListParAng (AntList, 1, time, SList->SUlist[suId-1]);
	lastTime = time; lastSu = suId;
	cosPA = cos(PA); sinPA = sin(PA);
	/* cos(2PA), sin(2PA) */
	gr = cosPA * cosPA - sinPA * sinPA;
	gi = cosPA * sinPA + sinPA * cosPA; 
	/* DEBUG - place to put break point */
	if (time>(2.5/24.)) {
	  PA = 0.0;
	} /* End DEBUG */
	/* DEBUG 
	gr = 1.0;
	gi = 0.0;*/
     }
      /* DEBUG
      inUV->buffer[indx+0] = ipol + qpol*cos(2*PA) + upol*sin(2*PA);
      inUV->buffer[indx+1] = 0.0;
      inUV->buffer[indx+3] = ipol - qpol*cos(2*PA) - upol*sin(2*PA);
      inUV->buffer[indx+4] = 0.0;
      inUV->buffer[indx+6] = -qpol*sin(2*PA) + upol*cos(2*PA);
      inUV->buffer[indx+7] = vpol;
      inUV->buffer[indx+9]  = -qpol*sin(2*PA) + upol*cos(2*PA);
      inUV->buffer[indx+10] = -vpol; */
      /* end DEBUG */
      /* Convert vis to circular */
      ConvertL2C(Muller, &inUV->buffer[indx], &outUV->buffer[indx], ncorr);
      
      /* Parallactic angle correction */
      jndx = indx;
      for (icor=0; icor<ncorr; icor++) {
	/* Correct RL */
	tr = outUV->buffer[jndx+2*incs];  
	ti = outUV->buffer[jndx+2*incs+1];
	outUV->buffer[jndx+2*incs]   = tr * gr - ti * gi;
	outUV->buffer[jndx+2*incs+1] = ti * gr + tr * gi;
	
	/* Correct LR */
	tr = outUV->buffer[jndx+3*incs];
	ti = outUV->buffer[jndx+3*incs+1];
	outUV->buffer[jndx+3*incs]   = tr * gr + ti * gi;
	outUV->buffer[jndx+3*incs+1] = ti * gr - tr * gi;
	jndx += 12;
      } /* end loop over correlation */
      outDesc->numVisBuff += 1;  /* count in buffer */
    } /* end loop over buffer */
    /* Write */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop over file */
  ObitErrLog(err); 
 cleanup:
  /* Cleanup */
  AntList  = ObitAntennaListUnref(AntList);
  ANTable = ObitTableANUnref(ANTable);  
   
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
  } /* end ObitUVLin2Cir */
