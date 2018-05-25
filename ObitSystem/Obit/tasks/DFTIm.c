/* $Id$  */
/* Obit task to DFT image a uv data set                               */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2018                                               */
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

#include "ObitThread.h"
#include "ObitFArray.h"
#include "ObitImageUtil.h"
#include "ObitUVUtil.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitUVUtil.h"
#include "ObitFITS.h"
#include "ObitSinCos.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* DFTImIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void DFTImOut (ObitInfoList* outList, ObitErr *err);
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
 ObitImage* setOutputImage (gchar *Source, olong iStoke, ObitInfoList *myInput, 
		    ObitUV* inData, ObitErr *err);

/* Loop over sources */
void doSources (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Loop over Poln */
void doChanPoln (gchar *Source, ObitInfoList* myInput, ObitUV* inData, 
		 ObitErr* err);

/* Write history */
void DFTImHistory (gchar *Source, gchar Stok, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitErr* err);
/* Image */
void doImage ( ObitUV* inData, ObitImage *outImage, ObitErr* err);

/* Program globals */
gchar *pgmName = "DFTIm";       /* Program name */
gchar *infile  = "DFTIm.in" ;   /* File with program inputs */
gchar *outfile = "DFTIm.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

/* Threaded function argument */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* Data descriptor */
  ObitUVDesc *myDesc;
  /* Data buffer */
  ofloat *buffer;
  /* Frequency scaling */
  ofloat *fscale;
  /* Array of X values per cell */
  ofloat* X;
  /* Array of Y values per cell */
  ofloat* Y;
  /* work array 1 per cell */
  ofloat *work1;
  /* work array 2 per cell */
  ofloat *work2;
  /* Cell accumulation */
  odouble *accum;
  /* Sum of weights */
  odouble sumWt;
  /* Number of cells */
  olong ncell;
  /* First vis, 1 rel */
  olong loVis;
  /* Highest vis, 1 rel */
  olong hiVis;
  /* Making a beam? */
  gboolean doBeam;
  /* thread number  */
  olong        ithread;
} DFTImFuncArg;
static olong MakeDFTImFuncArgs (ObitThread *thread, ObitUVDesc *myDesc,
				ofloat *buffer, ofloat *fscale,
				ofloat *X, ofloat *Y, olong ncell,
				gboolean doBeam,
				DFTImFuncArg ***ThreadArgs);
static void KillDFTImFuncArgs (olong nargs, DFTImFuncArg **ThreadArgs);
static gpointer ThreadDFTIm (gpointer arg);

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Obit task to image a uv data set                                     */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;
 
   /* Startup - parse command line */
  err = newObitErr();
  myInput = DFTImIn (argc, argv, err);
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

  /* Process */
  doSources (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* DFTImIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "DFTImIn";

  /* error checks */
  if (err->error) return list;

  /* Make default inputs/outputs InfoList */
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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "DataType", OBIT_string, dim, strTemp);
      
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
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end DFTImIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: DFTIm -input file -output ofile [args]\n");
    fprintf(stderr, "DFTIm Obit task to dft image data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def DFTIm.in\n");
    fprintf(stderr, "  -output output result file, def DFTIm.out\n");
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
    fprintf(stderr, "  -BChan first channel to copy\n");
    fprintf(stderr, "  -EChan highest channel to copy\n");
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
/*     outClass  Str [6]    output AIPS image class  [sClean]             */
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     UVRange   Flt (2)    Range n uv plane in klambda, def=all          */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*     subA      Int (1)    Subarray, def=1                               */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat farray[10];
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
  ObitInfoListPut (out, "outDType", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "DFTIm.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "DFTImName";
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
  strTemp = "DFTImOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "DFTImOut";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file class */
  strTemp = "Class ";
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

  /* Range in uv plane in klambda, 0=>all */
  dim[0] = 2;dim[1] = 1;
  farray[0] = 0.0; farray[1] = 0.0;
  ObitInfoListPut (out, "UVRange", OBIT_float, dim, farray, err);
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
  gchar Stokes[10];
  /*ofloat ftemp;
    gboolean *booTemp, btemp, do3D;
    olong itemp;*/
    gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Default Stokes is 'I' */
  strcpy (Stokes, "I   ");
  ObitInfoListGetTest(myInput, "Stokes",  &type, dim, Stokes);
  if (Stokes[0]==' ') Stokes[0] = 'I';
  ObitInfoListAlwaysPut(myInput, "Stokes",  type, dim, Stokes);

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

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
  olong         Aseq, disk, cno, nvis;
  gchar        *Type, *strTemp, inFile[129];
  oint         doCalib;
  gchar        Aname[13], Aclass[7], *Atype = "UV";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doCalSelect;
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
    
    /* define object  */
    nvis = 1;
    ObitUVSetAIPS (inData, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

    /* define object */
    nvis = 1;
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
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);
 
  /* Ensure inData fully instantiated and OK */
  ObitUVFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  /* Set number of vis per IO */
  nvis = 1000;  /* How many vis per I/O? */
  nvis =  ObitUVDescSetNVis (inData->myDesc, myInput, nvis);
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (inData->info, "nVisPIO", OBIT_long, dim,  &nvis);

  return inData;
} /* end getInputData */
/*----------------------------------------------------------------------- */
/*  Set output info on uvdata, create output image                        */
/*  Sets AIPS like Name,CLASS,seq info even for FITS files                */
/*  One output image per requested Stokes                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      iStoke    Stokes number (1-rel), I, Q, U, V, R, L                 */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to accept parameters defining output             */
/*      err       Obit Error stack                                        */
/*   Return:                                                              */
/*      outImage  Output image depending on Stokes request                */
/*-----------------------------------------------------------------------  */
 ObitImage* setOutputImage (gchar *Source, olong iStoke, ObitInfoList *myInput, 
		    ObitUV* inData, ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  ObitIOType IOType;
  olong      i, n, Aseq, disk, cno, nx, ny;
  gchar     *Type, *strTemp, outFile[129], *outName, *outF;
  ofloat    xCells, yCells;
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong     blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong     trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean  exist;
  gchar     tname[129], *chStokes="IQUVRL";
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputImage";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Image info */
  xCells = yCells = 1.0;
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &yCells);
  if (yCells <= 0.0) yCells = xCells;
  nx = ny = 10;
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  if (ny <= 0) ny = nx;

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image %cPol",chStokes[iStoke-1]);
  outImage = newObitImage(tname);
  ObitImageUtilUV2ImageDesc(inData->myDesc, outImage->myDesc, FALSE, 1); /* Init with UV Info */
  /* Set image info */
  outImage->myDesc->naxis = 4;
  outImage->myDesc->inaxes[0] = nx;
  outImage->myDesc->inaxes[1] = ny;
  outImage->myDesc->inaxes[2] = 1;
  outImage->myDesc->inaxes[3] = 1;
  outImage->myDesc->cdelt[0]  = -xCells/3600.;
  outImage->myDesc->cdelt[1]  =  yCells/3600.;
  outImage->myDesc->crval[0]  = inData->myDesc->crval[inData->myDesc->jlocr];
  outImage->myDesc->crval[1]  = inData->myDesc->crval[inData->myDesc->jlocd];
  outImage->myDesc->crpix[0]  = nx/2.0;
  outImage->myDesc->crpix[1]  = ny/2.0;
  strncpy(outImage->myDesc->ctype[0], "RA---SIN", 8);
  strncpy(outImage->myDesc->ctype[1], "DEC--SIN", 8);
  strncpy(outImage->myDesc->object, Source, 8);
    
  /* FILE type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) g_snprintf (tname, 100, "%s", outName);
    else g_snprintf (tname, 100, "%s%s", Source, outName);
    /* If no name use input name */
    if ((tname[0]==' ') || (tname[0]==0)) {
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
      g_snprintf (tname, 100, "%s", strTemp);
    }
      
    IOType = OBIT_IO_AIPS;  /* Save file type */
    /* input AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    strncpy (Aname, tname, 13); 
    Aname[12] = 0;
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "xClass", 7);
    }
    Aclass[0] = chStokes[iStoke-1];
    Aclass[6] = 0;

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
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
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&outF);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = outF[i]; tname[i] = 0;
    /* If blank use ".uvtab" */
    if ((tname[0]==' ') || (tname[0]==0)) g_snprintf (tname, 128, ".uvtab");
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (outFile, 128, "%c%s", chStokes[iStoke-1], tname);
    else g_snprintf (outFile, 128, "%s%c%s", Source, chStokes[iStoke-1], tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* Give output Image name */
    Obit_log_error(err, OBIT_InfoErr, "Output FITS image %s on disk %d ",
		    outFile, disk);

    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS Image %s on disk %d", outFile, disk);
    
    /* Make up AIPS-like name, class...  */
    if (strncmp (Source, "    ", 4))
      strncpy (Aname, Source, 13);
    else
      strncpy (Aname, outFile, 13);
    strncpy (Aclass, "XMap", 7);
    Aclass[0] = chStokes[iStoke-1];  /* Stokes type as first char */
    Aseq = 1;
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outImage;
  }
  
  /* Open/close to fully instantiate input and see if it's OK */
  ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);
  ObitImageClose (outImage, err);
  if (err->error>0) Obit_traceback_val (err, routine, outImage->name, outImage);

  /* Copy Field info to inData InfoList */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(inData->info, "imSeq",  OBIT_long, dim, &Aseq);
  ObitInfoListAlwaysPut(inData->info, "imDisk", OBIT_long, dim, &disk);
  ObitInfoListAlwaysPut(inData->info, "imFileType", OBIT_long, dim, &IOType);
  dim[0] = strlen (Aname);
  ObitInfoListAlwaysPut(inData->info, "imName", OBIT_string, dim, Aname);
  dim[0] = strlen (Aclass);
  ObitInfoListAlwaysPut(inData->info, "imClass", OBIT_string, dim, Aclass);

  return outImage;
} /* end setOutputImage */

/*----------------------------------------------------------------------- */
/*  Loop over selected sources, these are all sources in the source table */
/*  with selection by souCode, Sources, timeRange and optionally,         */
/*   prior processing                                                     */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doSources  (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  gchar        Source[17];
  ObitSourceList* doList;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         maxlen, isource, failed=0, good=0;
  gboolean     isBad = FALSE;
  gchar        *Fail="Failed  ", *Done="Done    ";
  gchar        *dataParms[] = {  /* Source selection*/
    "Sources", "souCode", "Qual", "timeRange", "FreqID",
    NULL
  };
  gchar *routine = "doSources";


  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Make sure selector set on inData */
  ObitUVOpen (inData, OBIT_IO_ReadCal, err);
  ObitUVClose (inData, err);
  
  /* Get source list to do */
  doList = ObitUVUtilWhichSources (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Loop over list of sources */
  for (isource = 0; isource<doList->number; isource++) {
    if (!doList->SUlist[isource]) continue; /* removed? */
    maxlen = MIN (16, strlen(doList->SUlist[isource]->SourceName));
    strncpy (Source, doList->SUlist[isource]->SourceName, maxlen);
    Source[maxlen] = 0;

    Obit_log_error(err, OBIT_InfoErr, " ******  Source %s ******", Source);
    ObitTrimTrail(Source);  /* remove trailing blanks */

    /* Save field name */
    dim[0] = 16; dim[1] = 1;
    ObitInfoListAlwaysPut (myInput, "FieldName", OBIT_string, dim, Source);

    /* Process source */
    doChanPoln (Source, myInput, inData, err);
    /* Allow up to 10 failures before first success or up to 10% of large run */
    if (err->error) {
      ObitErrLog(err); /* Show failure messages */
      failed++;
      isBad = TRUE;
      if (((failed>=10) && (good<=0)) || 
	  (((failed>=10)&&(failed>0.1*doList->number))) ||
	  (doList->number<=1)) {  /* Only one? */
	/* This isn't working - Give up */
	Obit_log_error(err, OBIT_Error, "%s: Too many failures, giving up", 
		       routine);
	return;
      }
    } else {
      isBad = FALSE;
      good++;
    } /* OK */


    /* Save processing summary - success or failure? */
    dim[0] = 8; dim[1] = 1;
    if (isBad) 
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Fail);
    else
      ObitInfoListAlwaysPut (myInput, "Status", OBIT_string, dim, Done);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);

    /* ReGet input uvdata */
    if (isource<(doList->number-1)) {
      inData = ObitUnref(inData);
      inData = getInputData (myInput, err);
      if (err->error) Obit_traceback_msg (err, routine, inData->name);
      
      /* Get input parameters from myInput, copy to inData */
      ObitInfoListCopyList (myInput, inData->info, dataParms);
      if (err->error) Obit_traceback_msg (err, routine, inData->name);
      
      /* Make sure selector set on inData */
      ObitUVOpen (inData, OBIT_IO_ReadCal, err);
      ObitUVClose (inData, err);
      
    } /* end reinit uvdata */
  } /* end source loop */

  doList = ObitSourceListUnref(doList);

}  /* end doSources */

/*----------------------------------------------------------------------- */
/*  Loop over frequencies and polarizations for a single source           */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doChanPoln (gchar *Source, ObitInfoList* myInput, ObitUV* inData, 
		 ObitErr* err)
{
  ObitUV       *tempData = NULL;
  ObitImage    *outImage=NULL;
  ofloat xCells, yCells;
  ObitInfoList* saveParmList=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong       ichan=1, bstok=1;
  olong        nx, ny;
  gboolean     doBeam=FALSE;
  gchar        Stokes[5], *chStokes=" IQUVRL";
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "UVRange", "timeRange", "BIF", "EIF", "BChan", "EChan", "Stokes", "subA",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "doBeam",
    NULL
  };
  gchar        *tmpParms[] = {  /* Imaging, weighting parameters */
    "nuGrid", "nvGrid", "MaxBaseline", "MinBaseline", 
    "xCells", "yCells","nx", "ny", "RAShift", "DecShift",
    "nxBeam", "nyBeam", "doCalSelect",
    NULL
  };
  gchar        *saveParms[] = {  /* Imaging, weighting parameters to save*/
    "xCells", "yCells","nx", "ny", "RAShift", "DecShift",
    "nxBeam", "nyBeam",
    NULL
  };
  gchar        *tmpName[] = {  /* Names to use for Image mosaic files */
    "imFileType", "imName", "imClass", "imDisk", "imSeq", "Sources",
    NULL
  };
  gchar *routine = "doChanPoln";

   /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Parameters used here */
   ObitInfoListGetTest(myInput, "doBeam", &type, dim, &doBeam);
 strcpy (Stokes, "I   ");
  ObitInfoListGetTest(myInput, "Stokes",  &type, dim, Stokes);
  xCells = yCells = 1.0;
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &yCells);
  if (yCells <= 0.0) yCells = xCells;
  nx = ny = 10;
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  if (ny <= 0) ny = nx;

  /* Place to save parameters */
  saveParmList = newObitInfoList ();

  
  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create temp file for data */
  tempData = newObitUVScratch (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* set selected Source  */
  dim[0] = strlen(Source); dim[1] = 1;
  ObitInfoListAlwaysPut (inData->info, "Sources", OBIT_string, dim, Source);
  /* Copy input data to output */
  ObitUVCopy (inData, tempData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  /* See if any data */
  if (tempData->myDesc->nvis<=0) {
    Obit_log_error(err, OBIT_InfoWarn,"NO data channel %d", ichan);
    goto endImage;  /* Go on to next */
  }
  ObitInfoListCopyList (myInput, tempData->info, tmpParms);
  ObitInfoListCopyList (inData->info, tempData->info, tmpName);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(tempData->info, "xCells", OBIT_float, dim, &xCells);
  ObitInfoListAlwaysPut(tempData->info, "yCells", OBIT_float, dim, &yCells);
  ObitInfoListAlwaysPut(tempData->info, "nx",     OBIT_long,  dim, &nx);
  ObitInfoListAlwaysPut(tempData->info, "ny",     OBIT_long,  dim, &ny);
  ObitInfoListAlwaysPut(tempData->info, "doBeam", OBIT_bool,  dim, &doBeam);
  
  /* Define output image */
  outImage = setOutputImage (Source, bstok, myInput, tempData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
  /* Save imaging parms for weighting - from defaults in mosaic creation */	
  ObitInfoListCopyList (tempData->info, saveParmList, saveParms);
    
  /* Image */    
  doImage (tempData,  outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, tempData->name);
    
 endImage:

  /* Do history */
  /* Make sure image created */
  Obit_return_if_fail((outImage!=NULL), err, 
		      "%s: No image generated", routine);
  DFTImHistory (Source, chStokes[bstok], myInput, inData, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Cleanup */
  tempData  = ObitUVUnref(tempData);
  if (saveParmList) saveParmList = ObitInfoListUnref(saveParmList);
  
}  /* end doChanPoln */

/*----------------------------------------------------------------------- */
/*  Write History for DFTIm                                               */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      Stoke     Stokes's parameter imaged I, Q, U, V                    */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void DFTImHistory (gchar *Source, gchar Stoke, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage,  ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "BIF", "EIF", "BChan", "EChan",  "UVRange",  "timeRange",
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol",  "PDVer", "xCells", "yCells", "nx", "ny", "RAShift", "DecShift", 
    "doBeam", "nThreads",
    NULL};
  gchar *routine = "DFTImHistory";

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
  } else { /* simply copy history */
    ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Add this programs history */
  ObitHistoryOpen (outHistory, OBIT_IO_ReadWrite, err);
  g_snprintf (hicard, 80, " Start Obit task %s ",pgmName);
  ObitHistoryTimeStamp (outHistory, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Write source and poln */
  if (Stoke==' ') Stoke = 'I';
  g_snprintf (hicard, 80, "%s Source = '%s', Stokes= '%c'", pgmName, Source, Stoke);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
  
} /* end DFTImHistory  */

/**
 * Make arguments for a Threaded ThreadDFTImFunc?
 * \param thread     ObitThread object to be used
 * \param myDesc     Input UV Descriptor
 * \param buffer     Input UV data buffer
 * \param fscale     Frequency scaling of u,v per channel
 * \param X          Array of x values per image cell
 * \param Y          Array of y values per image cell
 * \param ncell      Number of image cells
 * \param doBeam     If true, making a beam
 * \param ThreadArgs [out] Created array of DFTImFuncArg, 
 *                   delete with KillDFTImFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeDFTImFuncArgs (ObitThread *thread, ObitUVDesc *myDesc,
				ofloat *buffer, ofloat *fscale,
				ofloat *X, ofloat *Y, olong ncell,
				gboolean doBeam,
				DFTImFuncArg ***ThreadArgs)

{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(DFTImFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(DFTImFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread    = ObitThreadRef(thread);
    (*ThreadArgs)[i]->myDesc    = ObitUVDescRef(myDesc);
    (*ThreadArgs)[i]->buffer    = buffer;
    (*ThreadArgs)[i]->fscale    = fscale;
    (*ThreadArgs)[i]->X         = X;
    (*ThreadArgs)[i]->Y         = Y;
    (*ThreadArgs)[i]->ncell     = ncell;
    (*ThreadArgs)[i]->doBeam    = doBeam;
    (*ThreadArgs)[i]->work1     = g_malloc0(ncell*sizeof(ofloat));
    (*ThreadArgs)[i]->work2     = g_malloc0(ncell*sizeof(ofloat));
    (*ThreadArgs)[i]->accum     = g_malloc0(ncell*sizeof(odouble));
    (*ThreadArgs)[i]->sumWt     = 0.0;
    (*ThreadArgs)[i]->hiVis     = 0;
    (*ThreadArgs)[i]->loVis     = 1;
    (*ThreadArgs)[i]->ithread   = i;
  }

  return nThreads;
} /*  end MakeInterpImageArgs */

/**
 * Delete arguments for ThreadDFTImFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of DFTImFuncArg
 */
static void KillDFTImFuncArgs (olong nargs, DFTImFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      if (ThreadArgs[i]->myDesc) ObitUVDescUnref(ThreadArgs[i]->myDesc);
      if (ThreadArgs[i]->work1)  g_free(ThreadArgs[i]->work1);
      if (ThreadArgs[i]->work2)  g_free(ThreadArgs[i]->work2);
      if (ThreadArgs[i]->accum)  g_free(ThreadArgs[i]->accum);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillDFTImFuncArgs */

 /* center of DFTIm looping
 * Callable as thread
 * \param arg Pointer to DFTImFuncArg argument with elements:
 * \li myDesc  Data descriptor
 * \li buffer  Data buffer
 * \li fscale  Frequency scaling per channel
 * \li X       Array of X values per cell
 * \li Y       Array of Y values per cell
 * \li work1   work array 1 per cell
 * \li work2   work array 2 per cell
 * \li accum   Image cell accumulation
 * \li ncell   Number of image cells
 * \li doBeam  Making a beam?
 * \li loVis   First vis (1-rel) number
 * \li hiVis   Highest vis (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadDFTIm (gpointer arg)
{
  /* Get arguments from structure */
  DFTImFuncArg *largs = (DFTImFuncArg*)arg;
  ObitUVDesc *myDesc   = largs->myDesc;
  ofloat     *buffer   = largs->buffer;
  ofloat     *fscale   = largs->fscale;
  ofloat     *X   = largs->X;
  ofloat     *Y   = largs->Y;
  ofloat     *work1   = largs->work1;
  ofloat     *work2   = largs->work2;
  odouble    *accum   = largs->accum;
  odouble    sumWt    = largs->sumWt;
  olong      ncell    = largs->ncell;
  gboolean   doBeam   = largs->doBeam;
  olong      loVis    = largs->loVis-1;
  olong      hiVis    = largs->hiVis;

  /* local */
  olong ivis, ichan, nchan, ipos;
  ofloat *u, *v, *vis, *vvis, us, vs, ph, amp, wt;

  if (loVis>hiVis) goto finish;
  /* initialize data pointers */
  nchan = myDesc->ncorr;
  u   = buffer+myDesc->ilocu + loVis*myDesc->lrec;
  v   = buffer+myDesc->ilocv + loVis*myDesc->lrec;
  vis = buffer+myDesc->nrparm + loVis*myDesc->lrec;
  ph  = 0.0;
  amp = 1.0;
  for (ivis=loVis; ivis<hiVis; ivis++) { /* loop over buffer */
    vvis = vis;
    /* Loop over channel */
    for (ichan=0; ichan<nchan; ichan++) {
      if (vvis[2]>0.0) {
	/* Scaled u, v * 2 pi */
	us = (*u)*fscale[ichan]*2*G_PI;
	vs = (*v)*fscale[ichan]*2*G_PI;
	if (!doBeam) {  /* Use data if not making a beam */
	  ph  = -atan2(vvis[1], vvis[0]);                  /* Datum phase */
	  amp = sqrt(vvis[0]*vvis[0] + vvis[1]*vvis[1]);   /* Datum amplitude */
	}
	wt  = vvis[2];                                 /* Datum weight */
	sumWt += wt;
	/* Set phase array */
	for (ipos=0; ipos<ncell; ipos++) work1[ipos] = ph - X[ipos]*us + Y[ipos]*vs;
	ObitSinCosVec(ncell, work1, work2, work1);    /* Cosines sse/avx cos don't work*/
	/* Loop over image cells accumulating */
	for (ipos=0; ipos<ncell; ipos++) accum[ipos] += work1[ipos]*amp*wt;
      } /* end if valid data */
      vvis += myDesc->incf;  /* Update pointer */
    } /* end channel loop */
    /* update data pointers */
    u   += myDesc->lrec;
    v   += myDesc->lrec;
    vis += myDesc->lrec;
  } /* endloop over buffer */

  largs->sumWt = sumWt;  /* Save sum of weights */
  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /*  end ThreadDFTIm */

/*----------------------------------------------------------------------- */
/*  DFT Imaging                                                           */
/*   Input:                                                               */
/*      inData    ObitUV to image                                         */
/*      outImage  Output image                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doImage (ObitUV* inData, ObitImage *outImage, ObitErr* err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOCode retCode=OBIT_IO_SpecErr;
  olong  nx, ny, ix, iy,  plane[5]={1,1,1,1,1};
  olong ipos, i, nxny;
  ofloat xCells, yCells, ycen, xcen;
  ofloat deltax, deltay, x, y, *X, *Y;
  odouble sumWt=0.0;
  ofloat *pixels=NULL;
  gboolean doBeam=FALSE;
  olong nTh, nVis, loVis, hiVis, nVisPerThread, nThreads;
  DFTImFuncArg **threadArgs;
  gchar *routine = "doImage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Parameters used here*/
  /* Image info */
  ObitInfoListGetTest(myInput, "doBeam", &type, dim, &doBeam);
  if (doBeam) Obit_log_error(err, OBIT_InfoErr, "Making Dirty Beam");
  ObitErrLog(err); 
  xCells = yCells = 1.0;
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &xCells);
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &yCells);
  if (yCells <= 0.0) yCells = xCells;
  deltax = (xCells/3600.) * DG2RAD;  /* to radians */
  deltay = (yCells/3600.) * DG2RAD;
  nx = ny = 10;
  ObitInfoListGetTest(myInput, "nx", &type, dim, &nx);
  ObitInfoListGetTest(myInput, "ny", &type, dim, &ny);
  if (ny <= 0) ny = nx;
  nxny = nx*ny;
  xcen = nx/2.; ycen = ny/2.;

  /* Work arrays */
  X     = g_malloc0(nxny*sizeof(ofloat));
  Y     = g_malloc0(nxny*sizeof(ofloat));
  /* Loop over image filling X,Y vectors */
  ipos = 0;
  for (iy=0; iy<ny; iy++) {
    y = (iy - ycen)*deltay;
    for (ix=0; ix<nx; ix++) {
      x = (ix - xcen)*deltax;
      X[ipos] = x;
      Y[ipos] = y;
      ipos++;
    } /* end X loop */
  } /* end Y loop */

  /* Open data */
  retCode=ObitUVOpen (inData, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Need frequency scaling array */
  ObitUVGetFreq (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Initialize Threading */
  nThreads = MakeDFTImFuncArgs (inData->thread, inData->myDesc, 
				inData->buffer, inData->myDesc->fscale, X, Y, nxny, doBeam,
				&threadArgs);
  
  /* Loop through data */
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    retCode = ObitUVRead (inData, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Divide up work */
    nVis = inData->myDesc->numVisBuff;
    nVisPerThread = nVis/nThreads;
    nTh = nThreads;
    if (nVis<100) {nVisPerThread = nVis; nTh = 1;}
    loVis = 1;
    hiVis = nVisPerThread;
    hiVis = MIN (hiVis, nVis);
    
    /* Set up thread arguments */
    for (i=0; i<nTh; i++) {
      if (i==(nTh-1)) hiVis = nVis;  /* Make sure do all */
      threadArgs[i]->loVis   = loVis;
      threadArgs[i]->hiVis   = hiVis;
      if (nTh>1) threadArgs[i]->ithread = i;
      else threadArgs[i]->ithread = -1;
      /* Update which Vis */
      loVis += nVisPerThread;
      hiVis += nVisPerThread;
      hiVis = MIN (hiVis, nVis);
    }
    /* Do operation possibly with threads */
    ObitThreadIterator (inData->thread, nTh, 
			(ObitThreadFunc)ThreadDFTIm,
			(gpointer**)threadArgs);
  
  } /* end loop over data file */
    /* Close UV */
  retCode = ObitUVClose (inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Get pixel array */
  pixels = g_malloc0(nxny*sizeof(ofloat));

  /* Extract Thread values to output */
  sumWt = 0.0;
  for (i=0; i<nThreads; i++) sumWt += threadArgs[i]->sumWt;
  Obit_return_if_fail((sumWt>0.0), err, 
		      "%s: No valid data", routine);
  sumWt = 1.0/sumWt;  /* Normalization */
  for (i=0; i<nThreads; i++) {
    for (ipos=0; ipos<nxny; ipos++) 
      pixels[ipos] += (ofloat)(threadArgs[i]->accum[ipos]*sumWt);
  }

  /* Write image */
  ObitImagePutPlane (outImage, pixels, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Cleanup */
  KillDFTImFuncArgs(nThreads, threadArgs);
  if (pixels) g_free(pixels);
  if (X)      g_free(X);
  if (Y)      g_free(Y);

  } /* end doImage */

