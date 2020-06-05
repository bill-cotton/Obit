/*Need 
  - Beam
*/

/* $Id$  */
/* Obit task to restore and/or flatten an image mosaic                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2020                                               */
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
#include "ObitImageMosaicMF.h"
#include "ObitImageUtil.h"
#include "ObitImageMF.h"
#include "ObitTableCCUtil.h"
#include "ObitDConCleanVisMF.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitFITS.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* RestoreIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void RestoreOut (ObitInfoList* outList, ObitErr *err);
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
/* Get input mosaic of facets */
ObitImageMosaic* getMosaic (ObitInfoList *myInput, ObitErr *err);
/* Create output image */
void setOutputImage (ObitInfoList *myInput, ObitImageMosaic *mosaic, 
		     ObitUV* inData, ObitErr *err);
/* Restore */
void Restore(ObitInfoList *myInput, ObitImageMosaic *mosaic, olong CCver, ObitErr *err);
/* Cross restore */
void XRestore(ObitInfoList *myInput, ObitImageMosaic *mosaic, olong CCver, ObitErr *err);
/* Flatten mosaic */
void Flatten(ObitImageMosaic *mosaic, ObitErr *err);

/** Convolve spectral CCs with a Gaussian. */
ObitFArray* ConvlCC(ObitImage *image, olong CCVer, olong iterm, olong ncomp,
			   ofloat factor, ofloat tmaj, ofloat tmin, ofloat tpa, 
			   ObitErr *err);
/** Cross convolve spectral CCs with a Gaussian. */
void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
	      ObitImage *out, ObitFArray *outGrid, ObitErr *err);

/** Apply Gaussian taper to uv grid. */
void GaussTaper (ObitCArray* uvGrid,  ObitImageDesc *imDesc,
			ofloat gparm[3]);

/* Write history */
void RestoreHistory (ObitInfoList* myInput, ObitUV* inData, 
		     ObitImage* outImage, ObitErr* err);

/* Program globals */
gchar *pgmName = "Restore";       /* Program name */
gchar *infile  = "Restore.in" ;   /* File with program inputs */
gchar *outfile = "Restore.out";   /* File to contain program outputs */
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
/*   Obit task restore/flatten/fit a CLEANed image                        */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImageMosaic *mosaic=NULL;
  ObitUV       *inData = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     doRestore=TRUE, doFlat=TRUE, doFit=TRUE, isMF=FALSE;
  ofloat       antSize=-1.0;
  olong        CCver=0;
  ObitErr      *err= NULL;
 
   /* Startup - parse command line */
  err = newObitErr();
  myInput = RestoreIn (argc, argv, err);
   if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  ObitInfoListGetTest (myInput, "doRestore", &type, dim, &doRestore);
  ObitInfoListGetTest (myInput, "doFlat",    &type, dim, &doFlat);
  ObitInfoListGetTest (myInput, "doFit",     &type, dim, &doFit);
  ObitInfoListGetTest (myInput, "antSize",   &type, dim, &antSize);
  ObitInfoListGetTest (myInput, "CCVer",     &type, dim, &CCver);

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get mosaic */
  mosaic = getMosaic (myInput, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Create output image on mosaic */
  setOutputImage(myInput, mosaic, inData, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Process */
  if (doRestore) {
    Restore  (myInput, mosaic, CCver, err);
    XRestore (myInput, mosaic, CCver, err);
  }
  if (doFlat) Flatten(mosaic, err);
   /* Is this an ImageMF? */
  isMF = ObitImageMFIsA(mosaic->images[0]);
  if (doFit && isMF) ObitImageMFFitSpec ((ObitImageMF*)mosaic->FullField, 
					 antSize, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* History */
  RestoreHistory (myInput, inData, mosaic->FullField,  err);
  /* show any messages and errors */
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  mosaic    = ObitUnref(mosaic);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  err = ObitErrUnref(err);
  
  return ierr;
} /* end of main */

ObitInfoList* RestoreIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "RestoreIn";

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
      
     } else if (strcmp(arg, "-outDType") == 0) { /* Image type AIPS or FITS */
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

    } else if (strcmp(arg, "-out2Seq") == 0) { /* AIPS uv sequence number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "out2Seq", OBIT_oint, dim, &itemp, err);
      
    } else if (strcmp(arg, "-out2Disk") == 0) { /* output uv disk number */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "out2Disk", OBIT_oint, dim, &itemp, err);
      
     } else if (strcmp(arg, "-out2Name") == 0) { /* AIPS uv outName */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "out2Name", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-out2Class") == 0) { /* AIPS uv outClass */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "out2Class", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-out2File") == 0) { /* out2File (FITS uv) */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "out2File", OBIT_string, dim, strTemp);

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
} /* end RestoreIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Restore -input file -output ofile [args]\n");
    fprintf(stderr, "Restore Obit task to SW wideband image/CLEAN data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Restore.in\n");
    fprintf(stderr, "  -output output result file, def Restore.out\n");
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
    fprintf(stderr, "  -outFile output image (FITS Image file\n");  
    fprintf(stderr, "  -outName output image (AIPS file name\n");
    fprintf(stderr, "  -outClass output image (AIPS file class\n");
    fprintf(stderr, "  -outSeq output image (AIPS file sequence\n");
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -out2File output uv FITS Image file\n");  
    fprintf(stderr, "  -out2Name output uv AIPS file name\n");
    fprintf(stderr, "  -out2Class output uv AIPS file class\n");
    fprintf(stderr, "  -out2Seq output uv AIPS file sequence\n");
    fprintf(stderr, "  -out2Disk output uv (AIPS or FITS) disk number (1-rel) \n");
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
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     out2File  Str [?]    output FITS uv file name [def "Image.fits"    */
/*     out2Name  Str [12]   output AIPS uv name  [no def]                 */
/*     out2Class Str [6]    output AIPS uv class  [Restore]                */
/*     out2Seq   Int        output AIPS uv  sequence no  [new]            */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     FOV       Flt (1)    Field of view in deg , NO DEFAULT (0.0)       */
/*     doFlat    Boo (1)    Make full field (flattened) image? def=True   */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  ofloat ftemp;
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
  strTemp = "Restore.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "RestoreName";
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
  strTemp = "RestoreOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "RestoreOut";
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

  /* output FITS UV file name */
  strTemp = "UV.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "out2File", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file name */
  strTemp = "UV";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "out2Name", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV file class */
  strTemp = "IMAGER";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "out2Class", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS UV sequence */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "out2Seq", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output AIPS or FITS UV disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "out2Disk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Field of view in deg def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "FOV", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Make full field (flattened) image?, def = TRUE */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doFlat", OBIT_bool, dim, &btemp, err);
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
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean *booTemp;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Copy doFlat to doFlatten */
  ObitInfoListGetP (myInput, "doFlat", &type, dim, (gpointer)&booTemp);
  ObitInfoListAlwaysPut (myInput, "doFlatten", OBIT_bool, dim, booTemp);

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

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Get input image mosaic of facets                                      */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*      Image Mosaic                                                      */
/*----------------------------------------------------------------------- */
ObitImageMosaic* getMosaic (ObitInfoList *myInput, ObitErr *err)
{
  ObitImageMosaic *mosaic=NULL;
  ObitImage    **image=NULL;
  ObitInfoType type;
  olong        Aseq, disk, cno,i=0, nmaps;
  gchar        *Type, *Type2, *strTemp, inFile[129], inRoot[129];
  gchar        Aname[13], Aclass[7], Aroot[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean     isMF=FALSE;
  gchar        name[101];
  gchar *routine = "getMosaic";

  /* error checks */
  if (err->error) return mosaic;
  g_assert (ObitInfoListIsA(myInput));

  /* How many fields? */
  nmaps = 1;
  ObitInfoListGetTest(myInput, "nfield", &type, dim, &nmaps);

 /* Create local image array */
  image = g_malloc(nmaps*sizeof(ObitImage*));

  /* File type - could be either AIPS or FITS use DataType2 (default DataType) */
  ObitInfoListGetP (myInput, "DataType",  &type, dim, (gpointer)&Type);
  ObitInfoListGetP (myInput, "DataType2", &type, dim, (gpointer)&Type2);
  /* Defauts to DataType */
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
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
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
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	
	/* define object */
	image[i] = newObitImage(name);
	ObitImageSetAIPS(image[i], OBIT_IO_byPlane, disk, cno, AIPSuser,  blc, trc, err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	
	/* Ensure image fully instantiated and OK */
	ObitImageFullInstantiate (image[i], TRUE, err);

	/* Check is this is an ImageMF, and create mosaic of correct type */
	if (mosaic==NULL) {
	  if (!strncmp (image[i]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
	  else                                                      isMF = FALSE;
	  /* Create image mosaic */
	  if (isMF) mosaic = (ObitImageMosaic*)newObitImageMosaicMF ("Mosaic", nmaps);
	  else      mosaic = newObitImageMosaic ("Mosaic", nmaps);
	} /* end create mosaic */
	/* Remake image as ImageMF if needed */
	if (isMF) {
	  image[i] = ObitImageUnref(image[i]);
	  image[i] = (ObitImage*)newObitImageMF(name);
	  ObitImageSetAIPS(image[i], OBIT_IO_byPlane, disk, cno, AIPSuser,  blc, trc, err);
	  /* Ensure image fully instantiated and OK */
	  ObitImageFullInstantiate (image[i], TRUE, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	}
	
	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, i, image[i], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
      } /* end loop over fields */
      
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
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	
	/* Ensure image fully instantiated and OK */
	ObitImageFullInstantiate (image[0], TRUE, err);

	/* Check is this is an ImageMF, and create mosaic of correct type */
	if (mosaic==NULL) {
	  if (!strncmp (image[i]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
	  else                                                      isMF = FALSE;
	  /* Create image mosaic */
	  if (isMF) mosaic = (ObitImageMosaic*)newObitImageMosaicMF ("Mosaic", nmaps);
	  else      mosaic = newObitImageMosaic ("Mosaic", nmaps);
	} /* end create mosaic */
	/* Remake image as ImageMF if needed */
	if (isMF) {
	  image[0] = ObitImageUnref(image[0]);
	  image[0] = (ObitImage*)newObitImageMF(name);
	  ObitImageSetFITS(image[0], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	  /* Ensure image fully instantiated and OK */
	  ObitImageFullInstantiate (image[0], TRUE, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	}
	/* Attach Image */
	ObitImageMosaicSetImage (mosaic, 0, image[0], err);
	if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
     } else { /* Multiple fields */
	
	/* Loop over fields */
	for (i=0; i<nmaps; i++) {
	  /* Set file name */
	  g_snprintf (inFile, 128, "%s%d",inRoot,i);
	  
	  /* define object */
	  g_snprintf (name, 100, "Input image %d",i+1);
	  image[i] = newObitImage(name);
	  ObitImageSetFITS(image[i], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	  
	  /* Ensure image fully instantiated and OK */
	  ObitImageFullInstantiate (image[i], TRUE, err);

	  /* Check is this is an ImageMF, and create mosaic of correct type */
	  if (mosaic==NULL) {
	    if (!strncmp (image[i]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
	    else                                                      isMF = FALSE;
	    /* Create image mosaic */
	    if (isMF) mosaic = (ObitImageMosaic*)newObitImageMosaicMF ("Mosaic", nmaps);
	    else      mosaic = newObitImageMosaic ("Mosaic", nmaps);
	  } /* end create mosaic */
	  /* Remake image as ImageMF if needed */
	  if (isMF) {
	    image[i] = ObitImageUnref(image[i]);
	    image[i] = (ObitImage*)newObitImageMF(name);
	    ObitImageSetFITS(image[0], OBIT_IO_byPlane, disk, inFile, blc, trc, err);
	    /* Ensure image fully instantiated and OK */
	    ObitImageFullInstantiate (image[i], TRUE, err);
	    if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	  }
	  /* Attach Image */
	  ObitImageMosaicSetImage (mosaic, i, image[i], err);
	  if (err->error) Obit_traceback_val (err, routine, "myInput", mosaic);
	} /* end loop over fields */
      }

   } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		     pgmName, Type2);
      return mosaic;
    }
  /* Field of view */
  mosaic->FOV = 0.0;
  ObitInfoListGetTest(myInput, "FOV", &type, dim, &mosaic->FOV);
   
  /* cleanup */
  if (image) {
    for (i=0; i<nmaps; i++) image[i] = ObitImageUnref(image[i]);
    g_free(image);
  }
  return mosaic;
} /* end getMosaic */
  
/*----------------------------------------------------------------------- */
/*  Set output image on mosaic                                            */
/*  Values:                                                               */
/*      myImage   InfoList with inputs                                    */
/*      mosaic    Mosaic on which to attach output Image                  */
/*      err       Obit error/message stack                                */
/*----------------------------------------------------------------------- */
void setOutputImage(ObitInfoList *myInput, ObitImageMosaic *mosaic, 
		    ObitUV* inData, ObitErr *err)
{
  ObitImage    *outImage=NULL;
  ObitInfoType type;
  ObitImageMF  *outMF, *inMF;
  gint32        dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong         j, k, Aseq, disk, cno, nx, ny;
  ofloat       Beam[3]={0.,0.,0.};
  gboolean     exist, isMF = FALSE;
  gchar        *outType, *strTemp=NULL, inFile[129];
  gchar        iAname[16], Aname[16], Aclass[8], *Atype = "MA";
  gchar        tname[101];
  gchar *routine = "setOutputImage";

  if (err->error) return;  /* existing error? */
  
  strcpy(iAname,"Default name");  /* Default for FITS in/AIPS out */
  ObitInfoListGet(myInput, "in2Name", &type, dim, iAname, err);
  iAname[dim[0]] = 0;
 
  /* Is this an ImageMF? */
  if (!strncmp (mosaic->images[0]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
  else                                                               isMF = FALSE;

  /* Override restoring beam */
  ObitInfoListGetTest (myInput, "Beam", &type, dim, Beam);

  ObitInfoListGet (myInput, "DataType", &type, dim, tname, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
  /* output File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&outType);
  /* Defaults to DataType */
  if ((outType==NULL) || (!strncmp(outType,"    ",4)))
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&outType);

  /* Output */
  if (!strncmp (outType, "AIPS", 4)) { /* AIPS output */
    /* Output image */
    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* AIPS name */
    for (k=0; k<12; k++) {Aname[k] = ' ';} Aname[k] = 0;
    ObitInfoListGet(myInput, "outName", &type, dim, Aname, err);
    Aname[dim[0]] = 0;
    /* Default */
    if (!strncmp(Aname,"     ",5)) strcpy (Aname, iAname);

    /* AIPS class */
    for (k=0; k<6; k++) {Aclass[k] = ' ';} Aclass[k] = 0;
    ObitInfoListGet(myInput, "outClass", &type, dim, Aclass, err);
    Aclass[dim[0]] = 0;
    /* Default */
    if (!strncmp(Aclass,"      ",6)) strcpy (Aclass, "Convol");

    /* AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 
    
    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, 
			   &exist, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /* Generate Object name from AIPS name */
    g_snprintf (tname, 100, "%s.%s:%d.%d", Aname, Aclass, Aseq, disk);
    if (isMF) outImage = (ObitImage*)newObitImageMF(tname);
    else      outImage = newObitImage(tname);

    /* reset BLC, TRC */
    for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}
    
    /* define image */
    ObitImageSetAIPS (outImage, OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    

  } else if (!strncmp (outType, "FITS", 4)) {  /* FITS output */
    /* Output image */ 
    /* FITS file name */
    for (j=0; j<128; j++) inFile[j] = 0;
    ObitInfoListGet(myInput, "outFile", &type, dim, inFile, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    inFile[128] = 0;
    ObitTrimTrail(inFile);  /* remove trailing blanks */
    
    /*  FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
    
    /*  Object name from FITS name */
    if (isMF) outImage = (ObitImage*)newObitImageMF(inFile);
    else      outImage = newObitImage(inFile);
    
    /* reset BLC, TRC */
    for (j=0; j<IM_MAXDIM; j++)  {blc[j] = 1; trc[j] = 0;}

    /* define image */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
                   pgmName, strTemp);
  }
  if (err->error) Obit_traceback_msg (err, routine, routine);

  /* set on mosaic */
  mosaic->FullField = outImage;

  /* Define image */
  mosaic->FOV = 0.0;
  ObitInfoListGetTest(myInput, "FOV", &type, dim, &mosaic->FOV);
  mosaic->xCells = 0.0;
  ObitInfoListGetTest(myInput, "xCells", &type, dim, &mosaic->xCells);
  mosaic->xCells /= 3600.; /* to degrees */
  if (mosaic->xCells<=0.0) mosaic->xCells = mosaic->images[0]->myDesc->cdelt[0];
  if (mosaic->images[0]->myDesc->cdelt[0]<0.0) mosaic->xCells = -fabs(mosaic->xCells);
  mosaic->yCells = 0.0;
  ObitInfoListGetTest(myInput, "yCells", &type, dim, &mosaic->yCells);
  mosaic->yCells /= 3600.; /* to degrees */
  if (mosaic->yCells<=0.0) mosaic->yCells = mosaic->images[0]->myDesc->cdelt[1];
  nx = (olong)(5.5 + 2.0*mosaic->FOV/fabs(mosaic->xCells));
  ny = (olong)(5.5 + 2.0*mosaic->FOV/mosaic->yCells);

  /* Copy descriptor from facet 1 */
  mosaic->FullField->myDesc = 
    ObitImageDescCopy(mosaic->images[0]->myDesc, mosaic->FullField->myDesc, err);
  /* Add MF stuff if needed */
  if (isMF) {
    inMF  = (ObitImageMF*)mosaic->images[0];
    ObitImageMFGetSpec(inMF, err);
    outMF = (ObitImageMF*)mosaic->FullField;
    ObitImageMFGetSpec(outMF, err);
    if (err->error) Obit_traceback_msg (err, routine, routine);
  }
  mosaic->FullField->myDesc->inaxes[0] = nx;
  mosaic->FullField->myDesc->inaxes[1] = ny;
  mosaic->FullField->myDesc->cdelt[0]  = mosaic->xCells;
  mosaic->FullField->myDesc->cdelt[1]  = mosaic->yCells;
  mosaic->FullField->myDesc->crpix[0]  = nx/2.+1.;
  mosaic->FullField->myDesc->crpix[1]  = ny/2.+1.;
  if ((Beam[0]>0.0) && (Beam[1]>0.0)){
    /* Override beam - to degrees */
    mosaic->FullField->myDesc->beamMaj = Beam[0]/3600.;
    mosaic->FullField->myDesc->beamMin = Beam[1]/3600.;
    mosaic->FullField->myDesc->beamPA  = Beam[2];
  }
  /* Release buffer to force recreate */
  mosaic->FullField->image  = ObitFArrayUnref(mosaic->FullField->image);
  /* Ensure image fully instantiated and OK */
  ObitImageFullInstantiate (mosaic->FullField, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, routine);
} /* end setOutputImage */

/*----------------------------------------------------------------------- */
/*  Write History for Restore                                             */
/*   Input:                                                               */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void RestoreHistory (ObitInfoList* myInput, 
		     ObitUV* inData, ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "DataType2","in2File", "in2Disk","in2Name","in2Class","in2Seq",
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "CCVer", "EComp", "Beam", "antSize",
    "FOV", "xCells", "yCells", 
    "doRestore", "doFlat", "doFit",
    "nThreads",
    NULL};
  gchar *routine = "RestoreHistory";

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

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end RestoreHistory  */

/** 
 * Restore components removed from the residual image(s)
 * \param myInput Task inputs structure
 * \param mosaic  The mosaic object to restore
 * \param CCver   CC Table version
 * \param err     Obit error stack object.
 */
void Restore(ObitInfoList *myInput, ObitImageMosaic *mosaic, 
	     olong CCver, ObitErr *err)
{
  ObitFArray *convl=NULL;
  ObitImage *image=NULL;
  ofloat Beam[3]={0.,0.,0.},factor=1.0;
  gboolean isMF = FALSE;
  olong iplane, field, num, nOrd, plane[5]={1,1,1,1,1};
  olong ncomp, maxComp, *EComp=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "Restore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Tell user */
  if (err->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Get EComp */
  ObitInfoListGetP(myInput, "EComp",  &type, dim, (gpointer)&EComp);
  if (EComp) maxComp = dim[0];  /* How many? For higher CC tables take all */
  else       maxComp = -1;      /* No list */

  /* Override restoring beam */
  ObitInfoListGetTest (myInput, "Beam", &type, dim, Beam);
  Beam[0] /=3600.; Beam[1] /=3600.; /* To degrees */

  /* Loop over fields */
  for (field = 0; field<mosaic->numberImages; field++) {

    /* How many components this field? 0=> all */
    if (field<maxComp) ncomp = EComp[field];
    else               ncomp = 0;

    /* which Image? */
    image = mosaic->images[field];

    /* Is this an ImageMF? */
    if (!strncmp (mosaic->images[0]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
    else                                                               isMF = FALSE;

    /* Restore Flux then individual channels */
    /* Convolve Gaussians */
    iplane = 0;
    convl = ConvlCC (image, CCver, iplane, ncomp, factor, Beam[0], Beam[1], Beam[2], err);
    if (convl==NULL) continue; /* Components? */
    /* Read image */
    plane[0] = 1;
    ObitImageGetPlane (image, NULL, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
   /* Sum */
    ObitFArrayAdd (image->image, convl, image->image);
    /* Rewrite */
    ObitImagePutPlane (image, NULL, plane, err);
    convl = ObitFArrayUnref(convl);
    if (err->error) Obit_traceback_msg (err, routine, mosaic->name);

    /* Is this an ImageMF? */
    if (!strncmp (mosaic->images[0]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
    else                                                               isMF = FALSE;

    if (isMF) {
      /* Update info */
      ObitImageMFGetSpec((ObitImageMF*)image, err);
      /* individual channels */
      num  = ((ObitImageMF*)image)->nSpec;
      nOrd = ((ObitImageMF*)image)->maxOrder;
      for (iplane=1; iplane<=num; iplane++) {
	/* Convolve Gaussians */
	convl = ConvlCC (image, CCver, iplane, ncomp, 
			 factor, Beam[0], Beam[1], Beam[2], err);
	if (convl==NULL) continue; /* Components? */
	/* Read image */
	plane[0] = iplane+1 + nOrd;
	ObitImageGetPlane (image, NULL, plane, err);
	if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
	/* Sum */
	ObitFArrayAdd (image->image, convl, image->image);
	/* Rewrite */
	ObitImagePutPlane (image, NULL, plane, err);
	convl = ObitFArrayUnref(convl);
	if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
      }
    } /* end isMF */
    /* Reset restoring Beam if needed */
    if ((Beam[0]>0.0) && (Beam[1]>0.0)) {
      ObitImageOpen (image, OBIT_IO_ReadWrite, err);
      image->myDesc->beamMaj = Beam[0];
      image->myDesc->beamMin = Beam[1];
      image->myDesc->beamPA  = Beam[2];
      image->myStatus =  OBIT_Modified; /* Force update */
      ObitImageClose (image, err);
      if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
    }
    /* Free image memory */
    image->image = ObitFArrayUnref(image->image);
  } /* end loop over fields */

} /* end Restore */

/**
 * Restore components removed from one field but also 
 * appearing in another.  Does brute force convolution.
 * Wideband imaging version - does all spectral planes.
 * Spectral orders higher than 0 are flux density weighted averages.
 * \param myInput Task inputs structure
 * \param mosaic  The Mosaic object to restore
 * \param CCver   CC Table version
 * \param err     Obit error stack object.
 */
void XRestore(ObitInfoList *myInput, ObitImageMosaic *mosaic, 
	      olong CCver, ObitErr *err)
{
  ObitImage *image1=NULL, *image2=NULL;
  olong ifield, jfield, iplane, num, nOrd, ncomps, ver, noParms, plane[5]={1,1,1,1,1};
  ofloat BeamTaper1=0.0, BeamTaper2=0.0, factor;
  ofloat gparm[3]={0.0,0.0,0.0}, Beam[3]={0.,0.,0.}, bmaj, bmin, bpa;
  gboolean isAuto, isMF=FALSE;
  ObitFArray *convl=NULL, *accum=NULL;
  olong ncomp, maxComp, *EComp=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType itype, type;
  ObitTableCC *inCC=NULL, *outCC=NULL;
  gchar *routine = "ObitDConCleanVisMFXRestore";

   /* error checks */
  if (err->error) return;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Tell user */
  if (err->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Cross Restoring components");
    ObitErrLog(err);  /* Progress Report */
  }

  /* Get EComp */
  ObitInfoListGetP(myInput, "EComp",  &type, dim, (gpointer)&EComp);
  if (EComp) maxComp = dim[0];  /* How many? For higher CC tables take all */
  else       maxComp = -1;      /* No list */

  /* Override restoring beam */
  ObitInfoListGetTest (myInput, "Beam", &type, dim, Beam);
  Beam[0] /=3600.; Beam[1] /=3600.; /* To degrees */


  /* Double loop over fields */
  for (jfield = 0; jfield<mosaic->numberImages; jfield++) {
    /* output image  */
    image1 = mosaic->images[jfield];
    /* Is this an ImageMF? */
    if (!strncmp (mosaic->images[0]->myDesc->ctype[2], "SPECLNMF", 8)) isMF = TRUE;
    else                                                               isMF = FALSE;

    if (isMF) {
      num  = ((ObitImageMF*)image1)->nSpec;
      nOrd = ((ObitImageMF*)image1)->maxOrder;
    } else {
      num = 1; nOrd = 0;
    }

    /* Get additional beam taper */
    ObitInfoListGetTest(image1->myDesc->info, "BeamTapr", &itype, dim, &BeamTaper1);
    /* Ignore this one if not zero */
    if (BeamTaper1>0.0) continue;

    /* Get accumulation array for image */
    accum = ObitFArrayCreate ("Accum", 2, image1->myDesc->inaxes);
    
    /* Restore Flux then individual channels */
    for (iplane=0; iplane<(num+1); iplane++) {
      /* Read image */
      if (iplane==0) plane[0] = 1;
      else plane[0] = 1+iplane+nOrd;
      ObitImageGetPlane (image1, accum->array, plane, err);

      /* Loop over others */
      for (ifield = 0; ifield<mosaic->numberImages; ifield++) {
	/* Only cross */
	if (ifield==jfield) continue;

	/* which Image? */
	image2 = mosaic->images[ifield];

	/* Overlap? */
	if (!ObitImageDescOverlap(image1->myDesc, image2->myDesc, err)) continue;

	/* How many components this field? 0=> all */
	if (ifield<maxComp) ncomp = EComp[ifield];
	else                ncomp = 0;

	/* Cross convolve Gaussians */
	/* FFT (2D, same grid) or direct convolution */
	isAuto = (mosaic->isAuto[ifield]>0) || (mosaic->isAuto[jfield]>0);
 	if (!isAuto && (!image1->myDesc->do3D && !image2->myDesc->do3D) &&
	    (fabs(image1->myDesc->crval[0]-image2->myDesc->crval[0])<0.01*fabs(image1->myDesc->cdelt[0])) &&
	    (fabs(image1->myDesc->crval[1]-image2->myDesc->crval[1])<0.01*fabs(image1->myDesc->cdelt[1]))) {
	  /* Just to be sure: */
	  if (ObitImageDescAligned(image1->myDesc, image2->myDesc, err)) {
	      /* Can use FFT */
	      ver     = CCver;
	      noParms = 0;
	      inCC    = newObitTableCCValue ("SelectedCC", (ObitData*)image2,
					     &ver, OBIT_IO_ReadOnly, noParms, 
					     err);
	      outCC = ObitTableCCUtilCrossTable (inCC, image2->myDesc, image1, &ncomps, err);
	      if ((ncomps>0) && (outCC!=NULL)) {
		/* Scaling factor  */
		factor = ObitDConCleanGetXRestoreBeam(image2->myDesc, image1->myDesc, 
						      gparm, &bmaj, &bmin, &bpa);
		/* Get additional beam taper - use for convolution  */
		ObitInfoListGetTest(image2->myDesc->info, "BeamTapr", &itype, dim, &BeamTaper2);
		if (BeamTaper2>0.0) {
		  bmaj = bmin = BeamTaper2; bpa = 0.0;
		}
		/* Override restoring beam? */
		if ((Beam[0]>0.0) && (Beam[1]>0.0) && (BeamTaper2<=0.0)) {
		  bmaj = Beam[0]; bmin = Beam[1]; bpa = Beam[2];
		}
		convl = ConvlCC (image1, outCC->tabVer, iplane, ncomp, factor, bmaj, bmin, bpa, err);
		if (convl==NULL) continue;
		if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
		/* Sum */
		ObitFArrayAdd (accum, convl, accum);
		
		convl = ObitFArrayUnref(convl);
	      }
	      inCC = ObitTableCCUnref(inCC);
	      if (outCC!=NULL) {
		ObitImageZapTable (image1, "AIPS CC", outCC->tabVer, err);
		if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
		outCC = ObitTableCCUnref(outCC);  /* Be sure to free memory */
	      }
	    } /* end of if grid aligned */
	} else { /* direct convolution */
	  XConvlCC (image2, CCver, iplane, image1, accum, err);
	  if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
	}
      } /* end inner loop over fields */
      /* Rewrite */
      ObitImagePutPlane (image1, accum->array, plane, err);
      if (err->error) Obit_traceback_msg (err, routine, mosaic->name);
    } /* end loop over planes */
    
    accum = ObitFArrayUnref(accum);
  } /* end outer loop over fields */
} /* end XRestore */

/**
 * Flatten multiple facets if needed
 * Does Flatten if FullField member of mosaic member is defined.
 * \param mosaic  The mosaic object to flatten
 * \param err     Obit error stack object.
 */
void Flatten(ObitImageMosaic *mosaic, ObitErr *err)
{
  ObitImageMosaicClassInfo* mosaicClass; 
  gchar *routine = "Flatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitImageMosaicIsA(mosaic));

  if ((mosaic->FullField!=NULL) && (mosaic->numberImages>1)) {
    /* Tell user */
    if (err->prtLv>1) {
      Obit_log_error(err, OBIT_InfoErr,"Flattening images");
      ObitErrLog(err);  /* Progress Report */
    }
    mosaicClass = (ObitImageMosaicClassInfo*)mosaic->ClassInfo;
    mosaicClass->ObitImageMosaicFlatten (mosaic, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, mosaic->name);

} /* end Flatten */

/** 
 * Convolve a set of Clean components with a beam returning an image array
 * Gaussian will be obtained from the CC Table unless overridden by tmaj etc.
 * 
 * \param image  Image with CC table and defines size of image grid and
 *               with beam to convolve with.
 * \param CCVer  CC table number
 * \param iterm  Select spectral term, 0=flux, higher is a spectral channel..
 * \param nncomp Max CC number, 0=>all
 * \param factor Scaling factor
 * \param tmaj   if > 0 then the major axis size (deg) of convolving Gaussian
 * \param tmin   if > 0 then the minor axis size (deg) of convolving Gaussian
 * \param tpa    Position angle (deg) of convolving Gaussian if tmaj>=0
 * \param err    Obit error stack object.
 * \return An array with the Clean components convolved, NULL if no components
 */
ObitFArray* ConvlCC(ObitImage *image, olong CCVer, olong iterm, olong nncomp, 
			   ofloat factor, ofloat tmaj, ofloat tmin, ofloat tpa, 
			   ObitErr *err)
{
  ObitIOCode retCode;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitImageDesc *imDesc = NULL;
  olong first, last, ncomp, ver, ndim, naxis[2], ddim[2];
  gchar *tabType = "AIPS CC";
  ofloat gparm[3], bmaj, bmin, bpa;
  ObitFArray *grid = NULL;
  /* ObitFArray *tempFArray=NULL; DEBUG */
  ObitCArray *uvGrid = NULL;
  ObitFFT *forFFT = NULL, *revFFT = NULL;
  gchar *routine = "ConvlCC";

  /* error checks */
  if (err->error) return grid;

  /* Open Image */
  image->extBuffer = TRUE;   /* No need for buffer */
  retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);
  
  /* Restoring beam */
  imDesc = image->myDesc;
  bmaj   = imDesc->beamMaj;
  bmin   = imDesc->beamMin;
  bpa    = imDesc->beamPA;
    
  /* Close */
  retCode = ObitImageClose (image, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);
  image->extBuffer = FALSE;   /* May need buffer later */

  /* Get CC table */
  ver = CCVer;
  tempTable = newObitImageTable (image, OBIT_IO_ReadWrite, tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) Obit_traceback_val (err, routine, image->name, grid);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid);

  /* Open and close to get header */
  ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, grid);

  /* Any components? */
  if (CCTable->myDesc->nrow<=0) 
    {CCTable = ObitTableCCUnref(CCTable); return grid;}

  /* Grid components */
  first = 0;  /* all */
  last = 0;
  if (nncomp>0) last = nncomp;
  /* Need routine to use product of flux with another parameter */
  /* Return gridded flux*term */
  /* Spectral or normal */
  if ((CCTable->noParms>4) && (iterm>0)) { /* Spectral */
    retCode = ObitTableCCUtilGridSpect (CCTable, 1, iterm, &first, &last, FALSE, 
					factor, 0.0, 1.0e20, imDesc, &grid, gparm, &ncomp, 
					err);
  } else { /* normal */
    retCode = ObitTableCCUtilGrid (CCTable, 1, &first, &last, FALSE, 
				   factor, 0.0, 1.0e20, imDesc, &grid, gparm, &ncomp, 
				   err);
  }
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, image->name, grid);
  
  /* Free CC table */
  CCTable = ObitTableCCUnref(CCTable);
  
  /* FFT to image plane */
  ObitFArray2DCenter (grid); /* Swaparoonie to FFT order */
  /* Make Output of FFT if needed */
  ndim = 2;
  naxis[0] = 1+grid->naxis[1]/2; naxis[1] = grid->naxis[0]; 
  if (uvGrid) uvGrid = ObitCArrayRealloc(uvGrid, ndim, naxis);
  else uvGrid = ObitCArrayCreate ("FFT output", ndim, naxis);
  /* Create Forward FFT or reuse if OK */
  ddim[0] = grid->naxis[1]; ddim[1] = grid->naxis[0];
  if ((!forFFT) || (ddim[0]!=forFFT->dim[0]) || (ddim[1]!=forFFT->dim[1])) {
    forFFT = ObitFFTUnref(forFFT);
    forFFT = newObitFFT("FFT:FTImage", OBIT_FFT_Forward, 
			OBIT_FFT_HalfComplex, 2, ddim);
  }
  /* FFT */
  ObitFFTR2C (forFFT, grid, uvGrid);
  /* Put the center at the center */
  ObitCArray2DCenter (uvGrid);
  
  /* Taper for Gaussian - Use restoring beam or Gaussians from table if any 
     add rotation of image */
  if (gparm[0]<0.0) { /* restoring beam */
    gparm[0] = bmaj;
    gparm[1] = bmin;
    gparm[2] = bpa + image->myDesc->crota[image->myDesc->jlocd];
  } else { /* Gaussians from table */
    gparm[0] = gparm[0];
    gparm[1] = gparm[1];
    gparm[2] = gparm[2] + image->myDesc->crota[image->myDesc->jlocd];
  }
  /* Check for override in call */
  if (tmaj>0.0) {
    gparm[0] = tmaj;
    gparm[1] = tmin;
    gparm[2] = tpa + image->myDesc->crota[image->myDesc->jlocd];
  }
  GaussTaper (uvGrid, imDesc, gparm);
  
  /* FFT back to image */
  ObitCArray2DCenter (uvGrid); /* Swaparoonie to FFT order */
  /* Create reverse FFT or reuse if OK */
  ddim[0] = grid->naxis[1]; ddim[1] = grid->naxis[0];
  if ((!revFFT) || (ddim[0]!=revFFT->dim[0]) || (ddim[1]!=revFFT->dim[1])) {
    revFFT = ObitFFTUnref(revFFT);
    revFFT = newObitFFT("FFT:FTuv", OBIT_FFT_Reverse, 
			OBIT_FFT_HalfComplex, 2, ddim);
  }
  /* FFT */
  ObitFFTC2R (revFFT, uvGrid, grid);
  /* Put the center at the center */
  ObitFArray2DCenter (grid);
  
  /* Cleanup */
  forFFT = ObitFFTUnref(forFFT);
  revFFT = ObitFFTUnref(revFFT);
  uvGrid = ObitCArrayUnref(uvGrid);

  return grid;
} /* end ConvlCC */

/** 
 * Cross convolve a set of Clean components with a beam returning an image array
 * \param in      Image with CC table and defines size of image grid and
 *                with beam to convolve with.
 * \param CCVer   CC table number
 * \param iterm   Select spectral term, 0=flux, higher are spectral channels
 * \param out     Image defining output grid
 * \param outGrid Grid onto which to accumulate
 * \param err     Obit error stack object.
 */
void XConvlCC(ObitImage *in, olong CCVer, olong iterm, 
	      ObitImage *out, ObitFArray *outGrid, ObitErr *err)
{
  ObitImageDesc *imDesc1=NULL, *imDesc2=NULL;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitFArray *list = NULL, *tmpArray = NULL;
  olong j, ver, ncomp, ndim, naxis[2];
  ofloat gparm[3], gauss[3], bmaj, bmin, bpa, sr, cr, cellx, celly;
  ofloat scale, BeamTaper1=0.0, BeamTaper2=0.0;
  gchar *tabType = "AIPS CC";
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType itype;
  gchar *routine = "ObitDConCleanVisMF:XConvlCC";

  /* error checks */
  if (err->error) return;
  
  imDesc1 = in->myDesc;
  imDesc2 = out->myDesc;
  
  /* Any overlap? */
  if (!ObitImageDescOverlap(imDesc1, imDesc2, err)) return;
  
  /* Get additional beam taper for output */
  ObitInfoListGetTest(imDesc2->info, "BeamTapr", &itype, dim, &BeamTaper2);
  /* Ignore this one if not zero */
  if (BeamTaper2>0.0) return;

  /* Get additional beam taper for input */
  ObitInfoListGetTest(imDesc1->info, "BeamTapr", &itype, dim, &BeamTaper1);

  /* Get CC table */
  ver = CCVer;
  tempTable = newObitImageTable (in, OBIT_IO_ReadWrite, tabType, &ver, err);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Get list from jfield */
  if (iterm>0)
    list = ObitTableCCUtilCrossListSpec (CCTable, imDesc1, imDesc2, 
					 gparm, &ncomp, iterm, err);
  else
    list = ObitTableCCUtilCrossList (CCTable, imDesc1, imDesc2, 
				     gparm, &ncomp, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Free CC table */
  CCTable = ObitTableCCUnref(CCTable);
  
  /* Anything to do? */
  if (ncomp<=0) {ObitFArrayUnref(list); return;}

  /* Setup for ifield if needed */
  if (!tmpArray) {  /* Create FArray (zero filled) */
    ndim = 2; naxis[0] = imDesc2->inaxes[0]; naxis[1] = imDesc2->inaxes[1];
    tmpArray = ObitFArrayCreate ("Image for CCs", ndim, naxis);
  }
  
  /* get restoring beam and scaling */
  scale = ObitDConCleanGetXRestoreBeam(imDesc1, imDesc2, gparm, &bmaj, &bmin, &bpa);

  /* Scale list flux if needed */
  if (scale!=1.0) {
    if (iterm>0) {  /* Spectral flux density */
      for (j=0; j<list->naxis[1]; j++)
	list->array[2+j*list->naxis[0]] = list->array[3+j*list->naxis[0]] * scale;
    } else {   /* Normal flux density */
      for (j=0; j<list->naxis[1]; j++)
	list->array[2+j*list->naxis[0]] *= scale;
    }
  }

  /* Actually convolve with imaging taper if given */
  if (BeamTaper1>0.0) {
    bmaj = BeamTaper1;
    bmin = BeamTaper1;
    bpa  = 0.0;
  } else {
    bmaj = imDesc1->beamMaj;
    bmin = imDesc1->beamMin;
    bpa  = imDesc1->beamPA;
  }
  
  cellx = imDesc1->cdelt[0];
  celly = imDesc1->cdelt[1];
  cr = cos ((bpa + imDesc1->crota[imDesc1->jlocd])*DG2RAD);
  sr = sin ((bpa + imDesc1->crota[imDesc1->jlocd])*DG2RAD);
  gauss[0] = ((cr*cr)/(bmin*bmin) + (sr*sr)/(bmaj*bmaj)) *
    cellx*cellx*4.0*log(2.0);
  gauss[1] =  ((sr*sr)/(bmin*bmin) + (cr*cr)/(bmaj*bmaj)) *
    celly*celly*4.0*log(2.0);
  gauss[2] = (1.0/(bmin*bmin) - 1.0/(bmaj*bmaj)) *
    sr*cr*fabs(celly*celly)*8.0*log(2.0);

  /* Convolve list to tmpArray */
  ObitFArrayConvGaus (tmpArray, list, ncomp, gauss);
  
  /* Accumulate */
  ObitFArrayAdd (outGrid, tmpArray, outGrid);
    
  /* Cleanup */
  list = ObitFArrayUnref(list);
  tmpArray = ObitFArrayUnref(tmpArray);
  
} /* end XConvlCC */

/**
 * Apply Gaussian taper to a half Plane Complex grid 
 * assumed in the form from an ObitFFT.
 * NOTE: the uv grid is different in Obit (FFTW) and AIPS.
 * \param uvGrid Grid to be tapered
 * \param imDesc Image descriptor for image of which uvGrid is the FFT.
 * \param Gaussian in units of degrees, bmaj, bmin, bpa
 */
void GaussTaper (ObitCArray* uvGrid, ObitImageDesc *imDesc,
		 ofloat gparm[3])
{
  ofloat dU, dV, UU, VV, texp;
  ofloat konst, xmaj, xmin, cpa, spa, b1, b2, b3, bb2, bb3;
  ofloat taper, norm, *grid, tx, ty;
  olong i, j, nx, ny, naxis[2];

  /* Image info - descriptor should still be valid */
  nx = imDesc->inaxes[imDesc->jlocr];
  ny = imDesc->inaxes[imDesc->jlocd];
  
  /* Normalization factor */
  norm = ((ofloat)nx) * ((ofloat)ny);
  tx = MAX (1.0/sqrt(1.1331), gparm[0]/fabs(imDesc->cdelt[imDesc->jlocr]));
  ty = MAX (1.0/sqrt(1.1331), gparm[1]/fabs(imDesc->cdelt[imDesc->jlocd]));

  norm = 1.1331 * tx * ty / norm;

  /* UV cell spacing */
  dU = RAD2DG /  (nx * fabs(imDesc->cdelt[imDesc->jlocr]));
  dV = RAD2DG /  (ny * fabs(imDesc->cdelt[imDesc->jlocd]));
  
  konst = DG2RAD * G_PI * sqrt (0.5) / 1.17741022;
  xmaj = gparm[0] * konst;
  xmin = gparm[1] * konst;
  cpa = cos (DG2RAD * (90.0 + gparm[2])); /* is this right? */
  spa = sin (DG2RAD * (90.0 + gparm[2]));
  b1 = -(((cpa*xmaj)*(cpa*xmaj)) + ((spa*xmin)*(spa*xmin)));
  b2 = -(((spa*xmaj)*(spa*xmaj)) + ((cpa*xmin)*(cpa*xmin)));
  b3 = - 2.0 * spa * cpa * (xmaj*xmaj - xmin*xmin);
  
  /* pointer to complex grid */
  naxis[0] = naxis[1] = 0;
  grid = ObitCArrayIndex(uvGrid, naxis);
  
  /* loop over uv array */  
  for (i=0; i<ny; i++) {
    VV = dV * (i-nx/2);
    UU = 0.0;
    bb2 = b2 * VV * VV;
    bb3 = b3 * VV;
    /* Loop down row computing, applying taper */
    for (j=0; j<1+nx/2; j++) {
      texp = b1 * UU * UU + bb2 + bb3 * UU;
      if (texp>-14.0) taper = norm * exp (texp);
      else  taper = 0.0;
      UU = UU + dU;
      grid[2*j]   *= taper;
      grid[2*j+1] *= taper;
    }
    grid += 2*uvGrid->naxis[0];
  }

} /* end GaussTaper */

