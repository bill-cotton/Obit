/* $Id$  */
/* Obit task to make ionospheric movies from an SN table             */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007,2010                                          */
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
#include "ObitImage.h"
#include "ObitUV.h"
#include "ObitTableSN.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableUtil.h"
#include "ObitUVSoln.h"
#include "ObitFArray.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* IonMovieIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void IonMovieOut (ObitInfoList* outList, ObitErr *err);
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
/* Create output image */
ObitImage* setOutput (ObitInfoList *myInput, olong ntime, ObitUV* inData, 
		      ObitErr *err);
/* Count times/Convert SN table into movie */
olong doMovie (ObitInfoList* myInput, ObitUV* inData, ObitImage *outImage, 
	       ObitErr* err);
/* Write history */
void IonMovieHistory (ObitInfoList* myInput, ObitUV* inData, 
		      ObitImage* outImage, ObitErr* err);
/* Convert accumulations to degrees */
void GetDeg (ObitFArray *real, ObitFArray *imag, ObitFArray *deg);


/* Program globals */
gchar *pgmName = "IonMovie";       /* Program name */
gchar *infile  = "IonMovie.in" ;   /* File with program inputs */
gchar *outfile = "IonMovie.out";   /* File to contain program outputs */
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
/*   Obit task to make ionospheric movie from an SN table                 */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  olong        ntime = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitImage    *outImage = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = IonMovieIn (argc, argv, err);
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

  /* Count times in SN table */
  ntime = doMovie (myInput, inData, NULL, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Create output image */
  outImage = setOutput (myInput, ntime, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Process */
  doMovie (myInput, inData, outImage, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* History */
  IonMovieHistory (myInput, inData, outImage, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  outImage  = ObitUnref(outImage);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* IonMovieIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "IonMovieIn";

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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
    } else if (strcmp(arg, "-BIF") == 0) { /* BIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "BIF", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-EIF") == 0) { /* EIF */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "EIF", OBIT_oint, dim, &itemp, err);
      
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
      
    } else if (strcmp(arg, "-outSeq") == 0) { /* AIPS image sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outSeq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-outDisk") == 0) { /* output image disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "outDisk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-outName") == 0) { /* AIPS image outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outName", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outClass") == 0) { /* AIPS image outClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outClass", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

     } else if (strcmp(arg, "-AIPSdir") == 0) { /* Single AIPS Directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "AIPSdirs", OBIT_string, dim, strTemp);
      /* Only one AIPS Directory */
      dim[0] = 1;dim[1] = 1;
      itemp = 1; /* number of AIPS directories (1) */
      ObitInfoListPut (list, "nAIPS", OBIT_oint, dim, &itemp, err);

     } else if (strcmp(arg, "-FITSdir") == 0) { /* Single FITS Directory */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "FITSdirs", OBIT_string, dim, strTemp);
      /* Only one FITS Directory */
      dim[0] = 1;dim[1] = 1;
      itemp = 1; /* number of FITS directories (1) */
      ObitInfoListPut (list, "nFITS", OBIT_oint, dim, &itemp, err);

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
} /* end IonMovieIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IonMovie -input file -output ofile [args]\n");
    fprintf(stderr, "IonMovie Obit task make movie from SN table \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IonMovie.in\n");
    fprintf(stderr, "  -output output result file, def IonMovie.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -AIPSuser AIPS user number, def 2 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input \n");
    fprintf(stderr, "  -outDType AIPS or FITS type for output\n");
    fprintf(stderr, "  -inFile input FITS UV file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input image (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BIF first IF to copy\n");
    fprintf(stderr, "  -EIF highest IF to copy\n");
    fprintf(stderr, "  -outFile output image (FITS Image file\n");  
    fprintf(stderr, "  -outName output image (AIPS file name\n");
    fprintf(stderr, "  -outClass output image (AIPS file class\n");
    fprintf(stderr, "  -outSeq output image (AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -AIPSdir single AIPS data directory\n");
    fprintf(stderr, "  -FITSdir single FITS data directory\n");
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
/*     inName    Str [12]   input AIPS uv name  [no def]                  */
/*     inClass   Str [6]    input AIPS uv class  [no def]                 */
/*     inSeq     Int        input AIPS uv sequence no  [no def]           */
/*     outDisk   Int        output AIPS or FITS image disk no  [def 1]    */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     outName   Str [12]   output AIPS image name  [no def]              */
/*     outClass  Str [6]    output AIPS image class  ["Movie"]            */
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     solnVer   Int (1)    SN table to use, 0=> highest                  */
/*     subA      Int (1)    Subarray, def=1                               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[3];
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
  ObitInfoListPut (out, "outDType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "IonMovie.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "IonMovieName";
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

  /* output FITS Image file name */
  strTemp = "IonMovieOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "IonMovieOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file class */
  strTemp = "Movie ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "outSeq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS Image disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Sources selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Sources", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* Stokes parameter to image */
  strTemp = "I   ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Stokes", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  SN table to use, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "solnVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* maxSCLoop */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "maxSCLoop", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Subarray */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "subA", OBIT_oint, dim, &itemp, err);
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
  /* ObitInfoType type;
     gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
     gchar *strTemp;
     ofloat ftemp;
     olong itemp;
     gchar *routine = "digestInputs"; */

  /* error checks */
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
/*       ObitUV with input data                                           */
/*----------------------------------------------------------------------- */
ObitUV* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitUV       *inData = NULL;
  ObitInfoType type;
  olong         Aseq, disk, cno, nvis=1000;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
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
  
  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Create output image                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      ntime     Number of positions along time axis                     */
/*      inData    ObitUV with input data                                  */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Return:                                                              */
/*      outImage  Output image                                            */
/*----------------------------------------------------------------------- */
ObitImage* setOutput (ObitInfoList *myInput, olong ntime, ObitUV *inData, 
		      ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  ObitIOType IOType;
  olong      i, n, Aseq, disk, cno, nx, ny;
  ofloat    xCells, yCells, solInt;
  gchar     Type[10], *strTemp, outFile[129], *outName, *outF;
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean  exist;
  gchar     tname[129], Sources[50], *today = NULL;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutput";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image");
  outImage = newObitImage(tname);

  /* Source name? */
  g_snprintf (Sources, 10, "        ");
  ObitInfoListGetTest (myInput, "Sources", &type, dim, &Sources);
  ObitTrimTrail (Sources);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetTest (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetTest (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) strcpy (Type, FITS);
  dim[0] = strlen (Type); dim[1] = 1;
  ObitInfoListAlwaysPut (myInput, "outDType", OBIT_string, dim, Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    /* Something in source name? */
    if ((outName==NULL) || (outName[0]==' ') || (outName[0]==0)) 
      g_snprintf (tname, 100, "%s%s", Sources, outName);
    else g_snprintf (tname, 100, "%s", outName);
    /* If no name use input name */
    if ((tname[0]==' ') || (tname[0]==0)) {
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
      g_snprintf (tname, 100, "%s", strTemp);
    }
      
    IOType = OBIT_IO_AIPS;  /* Save file type */
    /* output AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    strncpy (Aname, tname, 13); Aname[12] = 0;
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Movie", 7);
    }
    if ((Aclass[0]==' ') && (Aclass[1]==' ') && (Aclass[2]==' ')) 
	strncpy (Aclass, "Movie", 7);
    Aclass[6] = 0;

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput",outImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* Tell about it */
    Obit_log_error(err, OBIT_InfoErr, "Output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);

    /* define object */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Sources, outName */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&outF);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = outF[i]; tname[i] = 0;
    /* If blank use ".fits" */
    if ((tname[0]==' ') || (tname[0]==0)) g_snprintf (tname, 128, ".fits");
    /* Something in source name? */
    if ((Sources==NULL) || (Sources[0]==' ') || (Sources[0]==0)) 
      g_snprintf (outFile, 128, "%s", tname);
    else g_snprintf (outFile, 128, "%s%s", Sources, tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);

    /* Write 32 bit float image */
    outImage->myDesc->bitpix = -32;
    
    /* Give output Image name */
    Obit_log_error(err, OBIT_InfoErr, "Output FITS image %s on disk %d ",
		    outFile, disk);

    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS Image %s on disk %d", outFile, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outImage;
  }

  /* Set image size, etc */
  nx = 230;
  ObitInfoListGetTest(myInput, "nx",  &type, dim, &nx);
  if (nx<1) nx = 230;
  ny = 230;
  ObitInfoListGetTest(myInput, "ny",  &type, dim, &ny);
  if (ny<1) ny = 230;
  xCells = 200.0;
  ObitInfoListGetTest(myInput, "xCells",  &type, dim, &xCells);
  if (xCells<1.0) xCells = 200.0;
  yCells = 200.0;
  ObitInfoListGetTest(myInput, "yCells",  &type, dim, &yCells);
  if (yCells<1.0) yCells = 200.0;
  solInt = 10.0 / 60.0;
  ObitInfoListGetTest(myInput, "solInt",  &type, dim, &solInt);
  if (solInt<1.0e-3) solInt = 10.0 / 60.0;
 
  outImage->myDesc->naxis = 3;
  /* E-W axis */
  outImage->myDesc->inaxes[0] = nx;
  g_snprintf (outImage->myDesc->ctype[0], 9, "East");
  outImage->myDesc->crval[0]  = 0.0;
  outImage->myDesc->cdelt[0]  = xCells;
  outImage->myDesc->crpix[0]  = outImage->myDesc->inaxes[0]/2;
  outImage->myDesc->crota[0]  = 0.0;
  /* N-S axis */
  outImage->myDesc->inaxes[1] = ny;
  g_snprintf (outImage->myDesc->ctype[1], 9, "North");
  outImage->myDesc->crval[1]  = 0.0;
  outImage->myDesc->cdelt[1]  = yCells;
  outImage->myDesc->crpix[1]  = outImage->myDesc->inaxes[1]/2;
  outImage->myDesc->crota[1]  = 0.0;
  /* Time axis */
  outImage->myDesc->inaxes[2] = ntime;
  g_snprintf (outImage->myDesc->ctype[2], 9, "Time");
  outImage->myDesc->crval[2]  = 0.0;
  outImage->myDesc->cdelt[2]  = solInt/1440.0;
  outImage->myDesc->crpix[2]  = 1.0;
  outImage->myDesc->crota[2]  = 0.0;

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE-1);
  if (today) g_free(today);
  /* Other stuff */
  strncpy (outImage->myDesc->object,     inData->myDesc->object,     IMLEN_VALUE-1);
  strncpy (outImage->myDesc->teles,      inData->myDesc->teles,      IMLEN_VALUE-1);
  strncpy (outImage->myDesc->instrument, inData->myDesc->instrument, IMLEN_VALUE-1);
  strncpy (outImage->myDesc->observer,   inData->myDesc->observer,   IMLEN_VALUE-1);
  strncpy (outImage->myDesc->obsdat,     inData->myDesc->obsdat,     IMLEN_VALUE-1);
  strncpy (outImage->myDesc->origin,     pgmName,                    IMLEN_VALUE-1);
  strncpy (outImage->myDesc->bunit,      "DEGREE  ",                 IMLEN_VALUE-1);
  outImage->myDesc->epoch  = inData->myDesc->epoch;
  outImage->myDesc->obsra  = inData->myDesc->obsra;
  outImage->myDesc->obsdec = inData->myDesc->obsdec;

  return outImage;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/*  Count times/Convert SN table into a Movie (3-D image)                 */
/*  Selection by Sources, FreqID, subA, timeRange, Antennas               */
/*  The two functions of this routine are to keep the logic consistent    */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*      outImage  Output image, if NULL only count times                  */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return:                                                              */
/*   Number of planes in output image, 0 on failure.                      */
/*----------------------------------------------------------------------- */
olong doMovie  (ObitInfoList* myInput, ObitUV* inData, ObitImage *outImage,
	       ObitErr* err)
{
  olong           iplane = 0;
  ObitInfoType    type;
  gint32          dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong           SNver, ANver, highVer, loop;
  olong            itemp, xcell, ycell, subA, iif, refAnt, numif, numpol;
  olong            bif, eif;
  odouble         antX, antY, timeEnd=-9999.0;
  ofloat          solInt, fblank = ObitMagicF(), *realP, *imagP;
  gboolean        doImage, doRCP, doLCP, wanted;
  ObitFArray      *accReal=NULL, *accImag=NULL, *outPlane=NULL;
  ObitTableSN     *SNTable=NULL;
  ObitTableSNRow  *SNRow=NULL;
  ObitTableAN     *ANTable=NULL;
  ObitAntennaList *antList=NULL;
  ObitImageDesc   *imdesc=NULL;
  ObitUVSel       *sel=NULL;
  ObitIOSize      IOBy = OBIT_IO_byPlane;
  gchar           Stokes[6];
  olong            blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong            trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar           *dataParms[] = {  /* Source selection*/
    "Sources", "Qual", "timeRange", "Antennas", "FreqID",
    NULL
  };
  gchar *routine = "doMovie";

  /* error checks */
  if (err->error) return iplane;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_val (err, routine, inData->name, iplane );

  /* Make sure selector set on inData */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);

  /* local pointers for structures */
  sel  = inData->mySel;

  /* Get AN table info */
  subA = 0;
  ObitInfoListGetTest(myInput, "subA", &type, dim, &subA);
  ANver = MAX (1, subA);
  ANTable = newObitTableANValue (inData->name, (ObitData*)inData, &ANver, 
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  antList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_val (err, routine, inData->name, iplane);
  ANTable = ObitTableANUnref (ANTable);

  /* Get SN table */
  itemp = 0;
  ObitInfoListGetTest(myInput, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  /* Which SN table? */
  highVer = ObitTableListGetHigh (inData->tableList, "AIPS SN");
  if (SNver==0) SNver = highVer;
  SNTable = newObitTableSNValue (inData->name, (ObitData*)inData, &SNver, 
				 OBIT_IO_ReadOnly, 0, 0, err);

  /* Get descriptive info - IFs */
  numif  = SNTable->numIF;
  bif = 1;
  ObitInfoListGetTest(myInput, "BIF", &type, dim, &bif);
  bif = MAX (0, bif-1);  /* 0-rel */
  eif = numif;
  ObitInfoListGetTest(myInput, "EIF", &type, dim, &eif);
  if (eif<=0) eif = numif;
  eif = MIN (numif-1, eif-1);  /* 0-rel */
  
  /* Polarization */
  numpol = SNTable->numPol;
  Stokes[0] = Stokes[1] = Stokes[2] = Stokes[3] = ' '; Stokes[4] = 0;
  ObitInfoListGetTest(myInput, "Stokes", &type, dim, Stokes);
  /* Do RCP or only poln? */
  doRCP = (Stokes[0]=='R') || (Stokes[0]==' ') || (numpol==1);
  doLCP = (Stokes[0]=='L') && (numpol>1); /* LCP? */

  /* averaging time */
  solInt = 10.0 / 60.0;
  ObitInfoListGetTest(myInput, "solInt",  &type, dim, &solInt);
  if (solInt<1.0e-3) solInt = 10.0 / 60.0;
  solInt /= 1440.0;  /* to days */

  /* Reference table if needed - once */
  refAnt = 0;
  ObitInfoListGetTest(myInput, "refAnt", &type, dim, &refAnt);
  if ((refAnt>0) && (outImage==NULL)) {
    ObitUVSolnRefAnt (SNTable, subA, &refAnt, err);
    if (err->error) Obit_traceback_val (err, routine, SNTable->name, iplane);
  }

  /* Open table */
  ObitTableSNOpen (SNTable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, SNTable->name, iplane);
 
  /* row to read SN table */  
  SNRow = newObitTableSNRow(SNTable);

  /* Output image? */
  doImage = ObitImageIsA(outImage);
  if (doImage) {
    IOBy = OBIT_IO_byPlane; /* Do I/O by plane */
    dim[0] = 1;
    ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, &IOBy, err);
    dim[0] = 7;
    ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, trc, err);
    outImage->extBuffer = FALSE;  /* Make sure it has buffer */

    /* Open Image, image on member image, an ObitFArray */
    ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
    if (err->error) Obit_traceback_val (err, routine, outImage->name, iplane);
    imdesc = outImage->myDesc;

    /* Clone accumulators */
    outPlane = outImage->image;
    accReal  = ObitFArrayCreate("accReal", outPlane->ndim, outPlane->naxis);
    accImag  = ObitFArrayCreate("accImag", outPlane->ndim, outPlane->naxis);
  } /* end initialize image */

  /* Loop over table */
  iplane = 0;
  for (loop=1; loop<=SNTable->myDesc->nrow; loop++) { 

    ObitTableSNReadRow (SNTable, loop, SNRow, err);
    if (err->error) Obit_traceback_val (err, routine, SNTable->name, iplane);
    if (SNRow->status<0) continue;  /* Skip deselected record */

    /* Initial time to header */
    if (doImage && (loop==1)) imdesc->crval[2] = SNRow->Time;

    /* Check if this record wanted */
    wanted = (SNRow->Time>=sel->timeRange[0]) && (SNRow->Time<=sel->timeRange[1]);
    wanted = wanted && ObitUVSelWantSour (sel, SNRow->SourID);
    wanted = wanted && ObitUVSelWantAnt (sel, SNRow->antNo);
    wanted = wanted && ((SNRow->FreqID==sel->FreqID) || (sel->FreqID<=0));
    wanted = wanted && ((SNRow->SubA==sel->SubA) || (sel->SubA<=0));
    if (!wanted) continue;

    /* End of this accumulation*/
    if (timeEnd<-1000.0) timeEnd = SNRow->Time + solInt;

    /* Accumulation finished? */
    if (SNRow->Time>timeEnd) {
      /* Write or count plane */
      iplane++;
      if (doImage) {
	GetDeg (accReal, accImag, outPlane);
	/* Add end time in pixel 0 */
	outPlane->array[0] = (ofloat)timeEnd;
	ObitImageWrite (outImage, NULL, err);
	if (err->error) Obit_traceback_val (err, routine, outImage->name, iplane);
      }

      /* If gap bigger than 2 solInt add 2 empty planes */
      if ((SNRow->Time-timeEnd)>=2.0*solInt) {
        iplane +=2;
	if (doImage) {
	  ObitFArrayFill (outPlane, fblank);
	  ObitImageWrite (outImage, NULL, err);
	  ObitImageWrite (outImage, NULL, err);
	  if (err->error) Obit_traceback_val (err, routine, outImage->name, iplane);
	}
      }
      timeEnd = SNRow->Time + solInt;  /* reset for next integration */
    } /* end accumumation done */

    /* Accumulate if writing image */
    if (doImage) {
      /* If VLA */
      if ( antList->isVLA) {
	/* Reproject Array onto Plains of San Augustin (lat=34.07874889 deg)
	   new X=E, Y=N */
	antX = antList->ANlist[(SNRow->antNo)-1]->AntXYZ[1];
	antY = antList->ANlist[(SNRow->antNo)-1]->AntXYZ[2] * 0.828268220 -
	  antList->ANlist[(SNRow->antNo)-1]->AntXYZ[0] * 0.560331827;
      } else {
	/*  don't know how to project */
	antX = antList->ANlist[(SNRow->antNo)-1]->AntXYZ[0];
	antY = antList->ANlist[(SNRow->antNo)-1]->AntXYZ[1];
      }

      /* Which cell (0-rel) */
      xcell = (olong)(-0.5 + imdesc->crpix[0] + (antX / imdesc->cdelt[0]));
      xcell = MAX (0, MIN (xcell, imdesc->inaxes[0]-1));
      ycell = (olong)(-0.5 + imdesc->crpix[1] + (antY / imdesc->cdelt[1]));
      ycell = MAX (0, MIN (ycell, imdesc->inaxes[1]-1));
 
     /* Pointers in accumulators */
      realP = accReal->array + xcell + ycell * imdesc->inaxes[0];
      imagP = accImag->array + xcell + ycell * imdesc->inaxes[0];

      /* Loop over selected IFs/poln */
      for (iif=bif; iif<=eif; iif++) { 
	if (doRCP && (SNRow->Weight1[iif]>0.0)  &&  (SNRow->Real1[iif]!=fblank)) {
	  *realP += SNRow->Real1[iif];
	  *imagP += SNRow->Imag1[iif];
	}
	if (doLCP && (SNRow->Weight2[iif]>0.0)  &&  (SNRow->Real2[iif]!=fblank)) {
	  *realP += SNRow->Real2[iif];
	  *imagP += SNRow->Imag2[iif];
	}
      }

    } /* end writing image */
  } /* end loop over table */  

  /* final accumulation */
  if (timeEnd>-1000.0) {
    iplane++;
    if (doImage) {
      GetDeg (accReal, accImag, outPlane);
      /* Add end time in pixel 0 */
      outPlane->array[0] = (ofloat)timeEnd;
      ObitImageWrite (outImage, NULL, err);
      if (err->error) Obit_traceback_val (err, routine, outImage->name, iplane);
    }
  }

 /* Close output image if writing */
  if (doImage) {
    ObitImageClose (outImage, err); /* Close */
    if (err->error) Obit_traceback_val (err, routine, outImage->name, iplane);
    /* Cleanup */
    outPlane = outImage->image = ObitFArrayUnref(outImage->image);
    accReal  = ObitFArrayUnref(accReal);
    accImag  = ObitFArrayUnref(accImag);
  }

  /* Cleanup */
  ObitTableSNClose (SNTable, err);
  if (err->error) Obit_traceback_val (err, routine, SNTable->name, iplane);
  antList = ObitAntennaListUnref( antList);
  SNTable = ObitTableSNUnref(SNTable);
  SNRow   = ObitTableSNRowUnref(SNRow);

  return iplane;
}  /* end doMovie */

/*----------------------------------------------------------------------- */
/*  Write History for IonMovie                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IonMovieHistory (ObitInfoList* myInput, ObitUV* inData, 
		      ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "Sources", "Qual", "Stokes", "subA", "Antennas",
    "FreqID", "BIF", "EIF", "timeRange", "solnVer",
    "xCells", "yCells", "nx", "ny", "SolInt", "refAnt", 
    NULL};
  gchar *routine = "IonMovieHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitImageIsA(outImage));

  /* Do Image history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
  } else if (outHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopy2Header (inHistory, outHistory, err);
  } else { /* simply copy history */
    ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* If output FITS - copy to header */
  if (outHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopy2Header (outHistory, outHistory, err);
  }

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IonMovieHistory  */

/**
 *  Convert accumulations to degrees
 *  out = RAD2DG*atan2(in2, in1),  
 *  if either is blanked or both zero, the result is blanked
 *  Quantizes at 0.1 deg
 * \param real  Input object with real part, zeroed on return
 * \param imag  Input object with imaginary part, zeroed on return
 * \param deg   Output array 
 */
void GetDeg (ObitFArray *real, ObitFArray *imag, ObitFArray *deg)
{
  olong i, itemp;
  ofloat ftemp, fblank = ObitMagicF();
   /* error checks */
  g_assert (ObitFArrayIsA(real));
  g_assert (ObitFArrayIsA(imag));
  g_assert (ObitFArrayIsCompatable(real, imag));
  g_assert (ObitFArrayIsCompatable(real, deg));

  for (i=0; i<real->arraySize; i++) {
    if ((real->array[i]!=fblank) && (imag->array[i]!=fblank) && 
	((real->array[i]!=0.0)&&(imag->array[i]!=0.0))) {
      ftemp = RAD2DG * atan2 (imag->array[i], real->array[i]);
      /* quantize at 0.1 deg */
      itemp = ftemp * 10.0;
      deg->array[i] = 0.1 * itemp;
    }
    else deg->array[i] = fblank;
  }

  /* Zero accumulators */
  ObitFArrayFill(real, 0.0);
  ObitFArrayFill(imag, 0.0);
 
} /* end ObitFArrayAdd */

