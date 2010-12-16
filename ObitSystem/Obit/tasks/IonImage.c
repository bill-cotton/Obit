/* $Id$  */
/* Obit task to image/CLEAN a uv data set with field-based cal        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2010                                          */
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
#include "ObitSkyModelVMIon.h"
#include "ObitImageMosaic.h"
#include "ObitImageUtil.h"
#include "ObitUVImager.h"
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitDConCleanVis.h"
#include "ObitUVSelfCal.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitDisplay.h"
#include "ObitIonCal.h"
#include "ObitUVImagerIon.h"
#include "ObitTableSNUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitUVPeelUtil.h"
#include "ObitTablePSUtil.h"
#include "ObitUVUtil.h"
#include "ObitFITS.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* IonImageIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void IonImageOut (ObitInfoList* outList, ObitErr *err);
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
ObitUV* setOutputUV (gchar *Source, ObitInfoList *myInput, ObitUV* inData, 
		     ObitErr *err);
/* Set output info on uvdata, create output image */
void setOutputData (gchar *Source, olong iStoke, ObitInfoList *myInput, 
		    ObitUV* inData, ObitImage **outImage, ObitErr *err);

/* Loop over sources */
void doSources (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Loop over Channels/Poln */
void doChanPoln (gchar *Source, ObitInfoList* myInput, ObitUV* inData, 
		 ObitErr* err);

/* Field based calibration */
void doFieldCal (gchar *Source, ObitInfoList* myInput, ObitUV* inData, 
		 ObitErr* err);

/* Field based Image */
void doFieldImage (gchar *Stokes, ObitInfoList* myInput, ObitUV* inData, 
		   ObitDConCleanVis *myClean, olong *nfield, olong **ncomp, 
		   ObitErr* err);

/* Subtract Stokes I model from data */
void subIPolModel (ObitUV* outData,  ObitSkyModel *skyModel, olong *selFGver, 
		   olong nfield, olong *ncomp, ObitErr* err);

/* Write history */
void IonImageHistory (gchar *Source, gchar Stok, ObitInfoList* myInput, 
		      ObitUV* inData, ObitImage* outImage, ObitUV* outData,
		      ObitErr* err);

/* Baseline dependent time averaging */
void BLAvg (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
	    ObitErr* err);

/* Program globals */
gchar *pgmName = "IonImage";       /* Program name */
gchar *infile  = "IonImage.in" ;   /* File with program inputs */
gchar *outfile = "IonImage.out";   /* File to contain program outputs */
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
/*   Obit task to image a uv data set with field-based calibration        */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = IonImageIn (argc, argv, err);
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

ObitInfoList* IonImageIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "IonImageIn";

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

     } else if (strcmp(arg, "-dispURL") == 0) { /* Display server URL */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "dispURL", OBIT_string, dim, strTemp);
      
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
} /* end IonImageIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: IonImage -input file -output ofile [args]\n");
    fprintf(stderr, "IonImage Obit task to image with field-based calibration\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def IonImage.in\n");
    fprintf(stderr, "  -output output result file, def IonImage.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType AIPS or FITS type for input image\n");
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
    fprintf(stderr, "  -dispURL display server URL \n");
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
/*     out2Class Str [6]    output AIPS uv class  [IonImage]                */
/*     out2Seq   Int        output AIPS uv  sequence no  [new]            */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     FOV       Flt (1)    Field of view in deg , NO DEFAULT (0.0)       */
/*     UVRange   Flt (2)    Range n uv plane in klambda, def=all          */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     Robust    Flt (1)    Briggs robust factor (AIPS version), def=0.0  */
/*     UVTaper   Flt (2)    Taper in uv plane in klambda in u, v, def=all */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=TRUE       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     doBand    Int (1)    If >0.5 apply bandpass cal.                   */
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     BPVer     Int (1)    Bandpass table version, 0=highest, def=0      */
/*     doPol     Boo (1)    Apply polarization calibration?, def=False    */
/*     doFull    Boo (1)    Make full field (flattened) image? def=True   */
/*     Catalog   Str (48)   Outlier catalog name, def 'NVSSVZ.FIT'        */
/*     OutlierDist Flt (1)  Maximum distance to add outlyers (deg), def=0 */
/*     OutlierFlux Flt (1)  Min. estimated outlier flux den. (Jy), def=0  */
/*     OutlierS  I Flt (1)  Spectral index to est. flux den., def=-0.7    */
/*     OutlierSize Int (1)  Size of outlier field, def=50                 */
/*     CLEANBox  Int[4,?]   Clean box, def=all                            */
/*     autoWindow Boo(1)    If true set windows automatically, def=FALSE  */
/*     Gain      Flt (1)    Clean gain, def=0.1                           */
/*     minFlux   Flt (1)    Clean minimum flux density, def=0             */
/*     Niter     Int (1)    Maximum # of CLEAN comp., def=No CLEAN        */
/*     minPatch  Int (1)    Clean Min. BEAM half-width, def=200           */
/*     Beam      Flt (3)    Clean beam maj, min, PA (", ", deg) def=fit   */
/*     CCFilter  Flt (2)    CC filter, [min. sum, radius in pix.], def=no */
/*     maxPixel  Int (1)    Max. pixels in inner cycle, def=50000         */
/*     ionVer    Int (1)    NI table version [1]                          */
/*     nZern     Int (1)    Number Zernike terms[5]                       */
/*     solInt    Flt (1)    Solution interval [1 min]                     */
/*     MaxRMS    Flt (1)    Target RMS residual [30]                      */
/*     UpdateInt Flt (1)    Interval (min) between cal. updates in apply  */
/*     Alpha     Flt (1)    default spectral index (0)                    */
/*     dispURL   Str(48)    Display derver URL                            */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp, iarray[4];
  ofloat ftemp, farray[3];
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
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input FITS file name */
  strTemp = "IonImage.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "IonImageName";
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
  strTemp = "IonImageOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "IonImageOut";
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

  /* Field of view in deg def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "FOV", OBIT_float, dim, &ftemp, err);
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

  /* Briggs robust factor (AIPS version), def=0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Robust", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Taper in uv plane in klambda in u, v, def=all */
  dim[0] = 2;dim[1] = 1;
  farray[0] = 0.0; farray[1] = 0.0;
  ObitInfoListPut (out, "UVTaper", OBIT_float, dim, farray, err);
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
  
  /* Make full field (flattened) image?, def = TRUE */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doFull", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Outlier catalog name, def = NVSSVZ.FIT" */
  strTemp = "NVSSVZ.FIT";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Catalog", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum distance to add outlyers (deg),, def=0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "OutlierDist", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Minimum estimated outlier flux density (Jy), def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "OutlierFlux", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Spectral index to estimate flux density, def = -0.7 */
  dim[0] = 1;dim[1] = 1;
  ftemp = -0.7; 
  ObitInfoListPut (out, "OutlierSI", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Size of outlyer field, def = 50 */
  dim[0] = 1;dim[1] = 1;
  itemp = 50; 
  ObitInfoListPut (out, "OutlierSize", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Clean box, def all */
  dim[0] = 4; dim[1] = 1;
  iarray[0] = iarray[1] = iarray[2] = iarray[3] = 0;
  ObitInfoListPut (out, "CLEANBox", OBIT_oint, dim, &iarray, err);

  /* autoWindow?, def= False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "autoWindow", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Clean gain, def = 0.1 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.1; 
  ObitInfoListPut (out, "Gain", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Clean minimum flux density, def = 0.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "minFlux", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum # of CLEAN comp., def = 0 (no clean) */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "Niter", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Clean Min. BEAM half-width, def = 200 */
  dim[0] = 1;dim[1] = 1;
  itemp = 200; 
  ObitInfoListPut (out, "minPatch", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Clean beam maj, min, PA (asec, asec, deg), def = 0 (fit) */
  dim[0] = 3;dim[1] = 1;
  farray[0] = 0.0; farray[1] = 0.0; farray[2] = 0.0;
  ObitInfoListPut (out, "Beam", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* CC filter, [min. sum flux, radius in pixels], def = no filter */
  dim[0] = 2;dim[1] = 1;
  farray[0] = 0.0; farray[1] = 0.0;
  ObitInfoListPut (out, "CCFilter", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  
  /* Maximum pixels searched in inner cycle, def = 50000 */
  dim[0] = 1;dim[1] = 1;
  itemp = 50000; 
  ObitInfoListPut (out, "maxPixel", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Subarray */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "subA", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* ionVer */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "ionVer", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* nZern */
  dim[0] = 1;dim[1] = 1;
  itemp = 5; 
  ObitInfoListPut (out, "nZern", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* UpdateInt  */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.5; 
  ObitInfoListPut (out, "UpdateInt", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* solInt */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1.0; 
  ObitInfoListPut (out, "solInt", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* MaxRMS */
  dim[0] = 1;dim[1] = 1;
  ftemp = 30.0; 
  ObitInfoListPut (out, "MaxRMS", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* default Spectral index */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.0; 
  ObitInfoListPut (out, "Alpha", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Display URL, def = "None" */
  strTemp = "None";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "dispRL", OBIT_string, dim, strTemp, err);
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
  gboolean *booTemp, btemp;
  olong itemp;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

   /* Default NField is 0 */
  itemp = 0;  type = OBIT_long; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListGet(myInput, "NField", &type, dim, &itemp, err);
  if (itemp<0) itemp = 0;
  ObitInfoListAlwaysPut (myInput, "NField", type, dim, &itemp);

  /* If   NField=1 doFull = FALSE */
  if (itemp==1) {
    btemp = FALSE; type = OBIT_bool; dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (myInput, "doFull", type, dim, &btemp);
  }

 /* Copy doFull to doFlatten */
  ObitInfoListGetP (myInput, "doFull", &type, dim, (gpointer)&booTemp);
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
  olong         Aseq, disk, cno, nvis, nThreads;
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
    
    /* define object */
    nvis = 1000;
    nThreads = 1;
    ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
    nvis *= nThreads;
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
    nvis = 1000;
    nThreads = 1;
    ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
    nvis *= nThreads;
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
/*  Create output uv data                                                 */
/*  Sets AIPS like Name,CLASS,seq info even for FITS files                */
/*  One output image per requested Stokes                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    Input ObitUV from which to clone output                 */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/* Returns the output uv data                                             */
/*----------------------------------------------------------------------- */
ObitUV* setOutputUV (gchar *Source, ObitInfoList *myInput, ObitUV* inData, 
 		     ObitErr *err)
{
  ObitUV    *outUV = NULL;
  ObitInfoType type;
  ObitIOType IOType;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, out2File[129], *out2Name, *out2F=NULL;
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129], *fullPath, *today=NULL;
  gchar     *routine = "setOutputUV";

  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output UV Object */
  g_snprintf (tname, 100, "output UV data");
  outUV = newObitUV(tname);
    
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
    /* Generate output name from Source, out2Name */
    ObitInfoListGetP (myInput, "out2Name", &type, dim, (gpointer)&out2Name);
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) g_snprintf (tname, 100, "%s", out2Name);
    else g_snprintf (tname, 100, "%s%s", Source, out2Name);
    /* If no name use input name */
    if ((tname[0]==' ') || (tname[0]==0)) {
      ObitInfoListGetP (myInput, "inName", &type, dim, (gpointer)&strTemp);
      g_snprintf (tname, 128, "%s", strTemp);
    }
      
    IOType = OBIT_IO_AIPS;  /* Save file type */
    /* input AIPS disk - default is outDisk */
    ObitInfoListGet(myInput, "out2Disk", &type, dim, &disk, err);
    if (disk<=0)
       ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* output AIPS sequence */
    ObitInfoListGet(myInput, "out2Seq", &type, dim, &Aseq, err);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    strncpy (Aname, tname, 13); Aname[12] = 0;
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "out2Class", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 6);
    } else { /* Didn't find */
      strncpy (Aclass, "IonImage", 6);
    }
    /* Default for blank */
    if (!strncmp (Aclass, "    ", 4)) strncpy (Aclass, "IonImage", 7);
    Aclass[6] = 0;

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "out2Seq", OBIT_oint, dim, &Aseq);
    } 

    /* Allocate catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
    /* define object */
    nvis = 1000;
    ObitInfoListGetTest(inData->info, "nVisPIO", &type, dim, &nvis);
    ObitUVSetAIPS (outUV, nvis, disk, cno, AIPSuser, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS UV data %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, out2Name */
    ObitInfoListGetP (myInput, "out2File", &type, dim, (gpointer)&out2F);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = out2F[i]; tname[i] = 0;
    /* If blank use ".uvtab" */
    if ((tname[0]==' ') || (tname[0]==0)) g_snprintf (tname, 128, ".uvtab");
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (out2File, 128, "!%s", tname);
    else g_snprintf (out2File, 128, "!%s%s", Source, tname);
    ObitTrimTrail(out2File);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "out2Disk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to outDisk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* Delete any previous version */
    fullPath = ObitFITSFilename (disk, out2File, err);
    /* Does filename exist? */
    if (ObitFileExist (fullPath, err) || err->error) {
      ObitFileZapFile (fullPath, err);
      g_free(fullPath);
    }

    /* define object */
    nvis = 1000;
    ObitInfoListGetTest(inData->info, "nVisPIO", &type, dim, &nvis);
    ObitUVSetFITS (outUV, nvis, disk, out2File, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output FITS UV data %s on disk %d", out2File, disk);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  /* Clone output uv from input UV *
  outUV->myDesc = ObitUVDescCopy(inData->myDesc, outUV->myDesc, err);
  outUV->myDesc->nvis = 0;  */
  /*ObitUVClone (inData, outUV, err);*/
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, 9);
  if (today) g_free(today);
 
  /* Full instantiation
  ObitUVFullInstantiate (outUV, FALSE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", outUV); */
  /*  outUV = ObitUVCopy (inData, outUV, err); DEBUG */

  return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Set output info on uvdata, create output image                        */
/*  Sets AIPS like Name,CLASS,seq info even for FITS files                */
/*  One output image per requested Stokes                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      iStoke    Stokes number (1-rel), I, Q, U, V                       */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to accept parameters defining output             */
/*   Output:                                                              */
/*      outImage  Output image depending on Stokes request                */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void setOutputData (gchar *Source, olong iStoke, ObitInfoList *myInput, 
		    ObitUV* inData, ObitImage **outImage, ObitErr *err)
{
  ObitInfoType type;
  ObitIOType IOType;
  olong      i, n, Aseq, disk, cno, lType;
  gchar     *Type, *strTemp, outFile[129], *outName, *outF;
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gboolean  exist;
  gchar     tname[129], *chStokes="IQUV";
  gchar     *routine = "setOutputData";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image %cPol",chStokes[iStoke-1]);
  *outImage = newObitImage(tname);
    
 /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  lType = dim[0];
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
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "xClass", 7);
    }
    Aclass[0] = chStokes[iStoke-1];

    /* if ASeq==0 create new, high+1 */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, "myInput");
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "outSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirAlloc(disk, AIPSuser, Aname, Aclass, Atype, Aseq, &exist, err);
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    
    /* define object */
    ObitImageSetAIPS ((*outImage), OBIT_IO_byPlane, disk, cno, AIPSuser, 
		      blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    Obit_log_error(err, OBIT_InfoErr, 
		   "Making output AIPS image %s %s %d on disk %d cno %d",
		   Aname, Aclass, Aseq, disk, cno);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&outF);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = outF[i]; tname[i] = 0;
    /* If blank use ".fits" */
    if ((tname[0]==' ') || (tname[0]==0)) g_snprintf (tname, 128, ".fits");
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (outFile, 128, "%c%s", chStokes[iStoke-1], tname);
    else g_snprintf (outFile, 128, "%s%c%s", Source, chStokes[iStoke-1], tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* define object */
    ObitImageSetFITS ((*outImage), OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
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
    return;
  }
  
  /* Copy Field info to inData InfoList */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(inData->info, "imSeq",  OBIT_long, dim, &Aseq);
  ObitInfoListAlwaysPut(inData->info, "imDisk", OBIT_long, dim, &disk);
  ObitInfoListAlwaysPut(inData->info, "imFileType", OBIT_long, dim, &IOType);
  dim[0] = 12;
  ObitInfoListAlwaysPut(inData->info, "imName", OBIT_string, dim, Aname);
  dim[0] = 6;
  ObitInfoListAlwaysPut(inData->info, "imClass", OBIT_string, dim, Aclass);

  return;
} /* end setOutputData */

/*----------------------------------------------------------------------- */
/*  Loop over sources in source list in myInput                           */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doSources  (ObitInfoList* myInput, ObitUV* inData, ObitErr* err)
{
  gchar        Source[17], *strTemp;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         i, isource, nsource, failed=0, good=0;
  gboolean     isBad = FALSE;
  gchar        *Fail="Failed  ", *Done="Done    ";
  gchar *routine = "doSources";


  /* For now only loop over list of sources */
  ObitInfoListGetP(myInput, "Sources",  &type, dim, (gpointer)&strTemp);
  nsource = dim[1];

  for (isource = 0; isource<nsource; isource++) {
    for (i=0; i<16; i++) Source[i] = ' '; Source[i] = 0;
    for (i=0; i<dim[0]; i++)  Source[i] = strTemp[i]; Source[i] = 0;
    /* Only accept blank name in the first Source */
    strTemp += dim[0];
    if (!strncmp (Source, "    ", 4) && (isource>0)) continue;
    Obit_log_error(err, OBIT_InfoErr, " ******  Source %s ******", Source);
    ObitTrimTrail(Source);  /* remove trailing blanks */

    doChanPoln (Source, myInput, inData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    /* Allow up to 10 failures before first success or up to 10% of large run */
    if (err->error) {
      ObitErrLog(err); /* Show failure messages */
      failed++;
      isBad = TRUE;
      if (((failed>=10) && (good<=0)) || 
	  (((failed>=10)&&(failed>0.1*nsource)))) {
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
    ObitTablePSSummary (inData, myInput, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  } /* end source loop */

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
  ObitUVImagerIon   *imager=NULL;
  ObitSkyModelVMIon *skyModel=NULL;
  ObitDConCleanVis  *myClean=NULL;
  ObitUV       *outData = NULL;
  ObitImage    *outField=NULL;
  ObitImage    *outImage[4]={NULL,NULL,NULL,NULL};
  ObitInfoList* saveParmList=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         ochan, ichan, nchan, chInc, chAvg, BChan, EChan, RChan, 
    inBChan, bchan, echan, istok, nstok, bstok, estok, nfield, *ncomp=NULL;
  gboolean     first, doFlat, doFCal, btemp, autoWindow, Tr=TRUE, do3D;
  olong        inver, outver, selFGver, plane[5] = {0,1,1,1,1};
  gchar        Stokes[5], *chStokes=" IQUV", *CCType = "AIPS CC";
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "UVRange", "timeRange", "UVTape",
    "BIF", "EIF", "BChan", "EChan", "subA", "FreqID", "souCode", "Qual",
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", "doPol",
    "Mode",
    NULL
  };
  gchar        *tmpParms[] = {  /* Imaging, weighting parameters */
    "doFull", "do3D", "FOV", "PBCor", "antSize", 
    "Catalog", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
    "Robust", "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "WtPower",
    "MaxBaseline", "MinBaseline", "rotate", "Beam",
    "NField", "xCells", "yCells", "nx", "ny", "RAShift", "DecShift",
    "nxBeam", "nyBeam", "dispURL", "Alpha", "doCalSelect", 
    NULL
  };
  gchar        *saveParms[] = {  /* Imaging, weighting parameters to save*/
    "NField", "xCells", "yCells","nx", "ny", "RAShift", "DecShift",
    "nxBeam", "nyBeam",
    NULL
  };
  gchar        *tmpName[] = {  /* Names to use for Image mosaic files */
    "imFileType", "imName", "imClass", "imDisk", "imSeq", "Sources",
    NULL
  };
  gchar        *CLEANParms[] = {  /* Clean parameters */
    "CLEANBox", "autoWindow", "Gain", "minFlux", "Niter", "minPatch", "Beam", 
    "Mode", "CCFilter", "maxPixel", "dispURL", "autoWindow",
    NULL
  };
  gchar        *SMParms[] = {  /* Sky model parameters */
    "ionVer", "UpdateInt", "prtLv", 
    NULL
  };
  gchar *NIlist[] = {"AIPS NI", NULL};
  gchar *routine = "doChanPoln";

   /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Total number of channels */
  nchan = inData->myDesc->inaxes[inData->myDesc->jlocf];

  /* Parameters used here */
  chInc = 1;
  ObitInfoListGetTest(myInput, "chInc",  &type, dim, &chInc);
  chAvg = 0;
  ObitInfoListGetTest(myInput, "chAvg",  &type, dim, &chAvg);
  /* Average everything = Continuum ? */
  if (chAvg<=0) {chAvg = nchan; chInc = nchan;}
  if (chInc<=0) chInc = chAvg;
  chInc = MAX (1, chInc);
  inBChan = 1;
  ObitInfoListGetTest(myInput, "BChan",  &type, dim, &inBChan);
  if (inBChan<=0) inBChan = 1;
  BChan = 1;
  EChan = nchan - inBChan+1;
  ObitInfoListGetTest(myInput, "EChan",  &type, dim, &EChan);
  if (EChan<=0) EChan = nchan;
  EChan = MIN (EChan, nchan) - inBChan+1;
  /* Update inputs */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(myInput, "EChan",  OBIT_long, dim, &EChan);
  RChan = 0;
  ObitInfoListGetTest(myInput, "RChan",  &type, dim, &RChan);
  RChan -= inBChan-1;
  RChan = MIN (RChan, EChan);
  RChan = MAX (BChan, RChan);
  strcpy (Stokes, "I   ");
  ObitInfoListGetTest(myInput, "Stokes",  &type, dim, Stokes);
  doFlat = TRUE;
  ObitInfoListGetTest(myInput, "doFlatten", &type, dim, &doFlat);
  autoWindow = FALSE;
  ObitInfoListGetTest(myInput, "autoWindow", &type, dim, &autoWindow);
  doFCal = TRUE;
  ObitInfoListGetTest(myInput, "doFCal", &type, dim, &doFCal);
  do3D = TRUE;
  ObitInfoListGetTest(myInput, "do3D", &type, dim, &do3D);
  
  /* Place to save parameters */
  saveParmList = newObitInfoList ();
  
  /* Number of stokes parameter, I, Q, U, V */
  nstok = 0;
  bstok = 1;
  if ((Stokes[0]=='I') || (Stokes[0]=='F') || (Stokes[0]==' ')) nstok = 1;
  if ((nstok==1) && (Stokes[1]=='Q')) nstok = 2;
  if ((nstok==2) && (Stokes[2]=='U')) nstok = 3;
  if ((nstok==3) && (Stokes[3]=='V')) nstok = 4;
  estok = bstok + nstok - 1;
  /* Single poln cases */
  if (Stokes[0]=='Q') {bstok=2; estok=2; nstok = 1;}
  if (Stokes[0]=='U') {bstok=3; estok=3; nstok = 1;}
  if (Stokes[0]=='V') {bstok=4; estok=4; nstok = 1;}
  /* Make sure Stokes OK */
  Obit_return_if_fail((estok>=bstok), err, 
			"%s: Problem with Stokes %s", routine, Stokes);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create output file for data */
  outData = setOutputUV (Source, myInput, inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy/select data to output file */
  dim[0] = 4;
  sprintf (Stokes, "    ");
  ObitInfoListAlwaysPut (inData->info, "Stokes", OBIT_string, dim, Stokes);
  /* set selected Source  */
  dim[0] = 16; dim[1] = 1;
  ObitInfoListAlwaysPut (inData->info, "Sources", OBIT_string, dim, Source);
  
  /* Copy or average input data to output */
  BLAvg (myInput, inData, outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Copy NI table to output if not calibrating */
  if (!doFCal) ObitUVCopyTables (inData, outData, NULL, NIlist, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Copy control info */
  ObitInfoListCopyList (myInput, outData->info, tmpParms);
  ObitInfoListCopyList (inData->info, outData->info, tmpName);

  /* Do field based calibration if requested */
  if (doFCal) {
    /* Do calibration */
    doFieldCal (Source, myInput, outData, err);
    /* Zap any iNdeX table */
    inver = -1;
    ObitUVZapTable (outData, "AIPS NX", inver, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  }

  /* Loop over channels */
  first = TRUE;
  ochan = (RChan - BChan) / chInc;
  for (ichan = RChan; ichan<=EChan; ichan+=chInc) {
    ochan++; /* output channel number */
    
    /* set selected channels */
    bchan = ichan; 
    echan = bchan + chAvg - 1;
    echan = MIN (echan, nchan);
    dim[0] = 1;
    ObitInfoListAlwaysPut (outData->info, "BChan", OBIT_long, dim, &bchan);
    ObitInfoListAlwaysPut (outData->info, "EChan", OBIT_long, dim, &echan);
    
    Obit_log_error(err, OBIT_InfoErr, " **** Start Channel %d - %d", bchan, echan);
    
    /* Loop over poln */
    for (istok=bstok; istok<=estok; istok++) {
      
      /* set selected Stokes  */
      dim[0] = 4;
      sprintf (Stokes, "%c   ", chStokes[istok]);
      ObitInfoListAlwaysPut (outData->info, "Stokes", OBIT_string, dim, Stokes);
      dim[0] = 1;
      ObitInfoListAlwaysPut (outData->info, "doCalSelect", OBIT_bool, dim, &Tr);
      
      /* Tell about it */
      Obit_log_error(err, OBIT_InfoErr, " ** Stokes %s", Stokes);

      /* Output image stuff */ 
      if (ichan==RChan) {
	/* Define output image(s) - potentially one per poln */
	setOutputData (Source, istok, myInput, outData, &outImage[istok-1], err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);
      } /* end first plane of Stokes */
      
      /* initialization */
      if (first ) {
	first = FALSE;
	
	/* Create Imager */
	imager = ObitUVImagerIonUnref(imager);
	imager = ObitUVImagerIonCreate("Ion imager", outData, err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);

	/* Create SkyModel */
	skyModel = ObitSkyModelVMIonUnref(skyModel);
	skyModel = ObitSkyModelVMIonCreate("Ion SkyModel", imager->mosaic);
	/* Get input parameters from myInput, copy to sykModel */
	ObitInfoListCopyList (myInput, skyModel->info, SMParms);

	/* Make CleanVis */
	myClean = ObitDConCleanVisUnref(myClean);
	myClean = ObitDConCleanVisCreate2("Clean Object", outData, 
					  (ObitUVImager*)imager, 
					  (ObitSkyModel*)skyModel, err);
	if (err->error) Obit_traceback_msg (err, routine, outData->name);
	
	/* Get input parameters from myInput, copy to myClean */
	ObitInfoListCopyList (myInput, myClean->info, CLEANParms);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);	

  	/* Save imaging parms for weighting - from defaults in mosaic creation */	
	ObitInfoListCopyList (outData->info, saveParmList, saveParms);
      }
      
      /* (Re)Set windows for Stokes I */
      if (istok==bstok) ObitDConCleanVisDefWindow((ObitDConClean*)myClean, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      
      /* Save imaging parms for weighting */	
      ObitInfoListCopyList (saveParmList, outData->info, saveParms);

      /* More Output image stuff */ 
      if (ichan==BChan) {
	/* Create output image(s) */
	if (doFlat) 
	  outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
	else
	  outField = ObitImageMosaicGetImage (myClean->mosaic, 0, err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);
		
	ObitImageUtilMakeCube (outField->myDesc, inData->myIO->myDesc, 
			       outImage[istok-1]->myDesc, 
			       Stokes, BChan, EChan, chInc, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	ObitImageFullInstantiate (outImage[istok-1], FALSE, err);
	outField = ObitImageUnref(outField);
	if (err->error) Obit_traceback_msg (err, routine, outImage[istok-1]->name);
	/* end of create output */
      } else if ((RChan>BChan) && (ichan==RChan)) { /* Restarting - should already exist */
	ObitImageFullInstantiate (outImage[istok-1], TRUE, err);
	if (err->error) Obit_traceback_msg (err, routine, outImage[istok-1]->name);
      }
      
      /* Automatic windowing  */
      btemp = autoWindow;
      /*** if (istok>bstok) btemp = Fl; */
      ObitInfoListAlwaysPut (myClean->info, "autoWindow", OBIT_bool, dim, &btemp);

      /* Set Stokes on SkyModel */
      ObitInfoListAlwaysPut (myClean->skyModel->info, "Stokes", OBIT_string, dim, Stokes);

      /* Do actual processing */
      doFieldImage (Stokes, myInput, outData, myClean, &nfield, &ncomp, err);
      if (err->error) Obit_traceback_msg (err, routine, outData->name);

      /* Subtract sky model from outData for I  */
      selFGver = -1;
      if ((istok==1) && (nstok==1) && (ichan==RChan)) {
	subIPolModel (outData, myClean->skyModel, &selFGver, nfield, ncomp, err);
      }
      if (err->error) Obit_traceback_msg (err, routine, outData->name);

      /* If 2D imaging concatenate CC tables */
      if ((myClean->nfield>1) && (!myClean->mosaic->images[0]->myDesc->do3D) && doFlat)
	ObitImageMosaicCopyCC (myClean->mosaic, err);
      
      /* Copy result to output */
      plane[0] = ochan;
      if (doFlat) 
	outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
      else { /* Copy the first image */
	outField = ObitImageMosaicGetImage (myClean->mosaic, 0, err);
	/* and copy the CC table */
	inver = 1;
	outver = plane[0];
	ObitDataCopyTable ((ObitData*)outField, (ObitData*)outImage[istok-1],
			   CCType, &inver, &outver, err);
      }
      ObitImageUtilInsertPlane (outField, outImage[istok-1], plane, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);

      /* For 2D imaging copy CC Table */
      if (!do3D && doFlat) {
 	inver   = 1;
	outver  = plane[0];
	ObitDataCopyTable ((ObitData*)outField, (ObitData*)outImage[istok-1],
			   CCType, &inver, &outver, err);
	if (err->error) Obit_traceback_msg (err, routine, outField->name);
      }
      outField = ObitImageUnref(outField);

      if (ncomp) g_free(ncomp); ncomp = NULL;
    } /* End Stokes loop */

    
  } /* end loop over channels */

  /* Do history */
  for (istok=bstok; istok<=estok; istok++) {
    /* Make sure image created */
    Obit_return_if_fail((outImage[istok-1]!=NULL), err, 
			"%s: No image generated", routine);

    if (istok==bstok)
      IonImageHistory (Source, chStokes[istok], myInput, inData, outImage[istok-1], 
		       outData, err);
    else
      IonImageHistory (Source, chStokes[istok], myInput, inData, outImage[istok-1], 
		       NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    outImage[istok-1]  = ObitUnref(outImage[istok-1]);
 } /* end History/Unref loop */
  
  
  /* Cleanup */
  if (myClean) {
    /* Leave facet images if not myClean->mosaic->doFull and myClean->mosaic->numberImages >1 */
    if (!(!myClean->mosaic->doFull && (myClean->mosaic->numberImages>1)))
      ObitImageMosaicZapImage (myClean->mosaic, -1, err); /* Delete mosaic members */
    if (doFlat) {  /* Delete flattened as well if not output */
      outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
      if (outField) outField = ObitImageZap(outField, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    }
    myClean  = ObitDConCleanVisUnref(myClean);
  }
  outData  = ObitUVUnref(outData);
  imager   = ObitUVImagerIonUnref(imager);
  skyModel = ObitSkyModelVMIonUnref(skyModel);

}  /* end doChanPoln */

/*----------------------------------------------------------------------- */
/*  Field based calibration                                               */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doFieldCal (gchar *Source, ObitInfoList* myInput, ObitUV* inUV, 
		 ObitErr* err)
{
  ObitIonCal *ionCal = NULL;
  gboolean autoWindow, do3D=TRUE, saveDoFull=FALSE, Tr=TRUE, Fl=FALSE;
  ObitInfoType type;
  ObitIOType IOType;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar        Aname[13], Aclass[7], *Type;
  olong         NField, itemp, Niter, Aseq=1, disk, ver;
  ofloat Gain, minFlux, seeing;
  gchar        *FCParms[] = {  /* Field based calibration parameters */
    "Catalog",  "OutlierDist",  "OutlierFlux",  "OutlierSI", "OutlierSize",  
    "solInt", "nZern", "FitDist", "MaxDist", "MinPeak", 
    "MaxWt", "MaxQual", "MaxRMS", "MinRat", "FCStrong", "prtLv", 
    NULL
  };
  gchar *routine = "doFieldCal";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inUV));
  
  /* Create IonCal */
  ionCal = newObitIonCal ("Ionospheric calibration");
  ObitIonCalSetData (ionCal, inUV);   /* Attach uv data */
  
  /* Copy control info */
  ObitInfoListCopyList (myInput, ionCal->info, FCParms);

  /* Save input do3D and set to TRUE on uv data */
  ObitInfoListGetTest(myInput, "do3D",  &type, dim, &do3D);
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "do3D", OBIT_bool, dim, &Tr);

  /* Special field based calibration parameters  used here and renamed */
  Niter = 200; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListGetTest(myInput, "FCNiter",  &type, dim, &Niter);
  ObitInfoListAlwaysPut(inUV->info, "Niter", OBIT_long, dim, &Niter);
  Gain = 0.1;
  ObitInfoListGetTest(myInput, "FCGain",  &type, dim, &Gain);
  ObitInfoListAlwaysPut(inUV->info, "Gain", OBIT_float, dim, &Gain);
  minFlux = 1.0;
  ObitInfoListGetTest(myInput, "FCminFlux",  &type, dim, &minFlux);
  ObitInfoListAlwaysPut(inUV->info, "minFlux",  OBIT_float, dim, &minFlux);
  autoWindow = FALSE;
  ObitInfoListGetTest(myInput, "autoWindow",  &type, dim, &autoWindow);
  ObitInfoListAlwaysPut(inUV->info, "autoWindow", OBIT_bool, dim, &autoWindow);
  NField = 0; itemp = 0;
  ObitInfoListGetTest(myInput, "NField",  &type, dim, &NField);
  ObitInfoListAlwaysPut(inUV->info, "NField", OBIT_long, dim, &itemp);
 

  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) IOType = OBIT_IO_AIPS;
  else IOType = OBIT_IO_FITS;
  disk = 1;
  ObitInfoListGetTest(myInput, "outDisk", &type, dim, &disk);

  /* Scratch file name info on inUV InfoList */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(inUV->info, "imSeq",  OBIT_long, dim, &Aseq);
  ObitInfoListAlwaysPut(inUV->info, "imDisk", OBIT_long, dim, &disk);
  ObitInfoListAlwaysPut(inUV->info, "imFileType", OBIT_long, dim, &IOType);
  dim[0] = 12;
  g_snprintf (Aname, 12, "TMP%s", Source);
  ObitInfoListAlwaysPut(inUV->info, "imName", OBIT_string, dim, Aname);
  dim[0] = 6;
  g_snprintf (Aclass, 6, "IC0000");
  ObitInfoListAlwaysPut(inUV->info, "imClass", OBIT_string, dim, Aclass);

  /* Save/reset any doFull */
  ObitInfoListGetTest(myInput, "doFull",  &type, dim, &saveDoFull);
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doFull", OBIT_bool, dim, &Fl);

  /* Calibrate */
  ObitIonCaldoCal (ionCal, err); 
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

 /* Reset doFull */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doFull", OBIT_bool, dim, &saveDoFull);

  /* Say which NI table to use; 0=> highest*/
  dim[0] = dim[1] = 1;
  ver = 0;
  ObitInfoListAlwaysPut(myInput, "ionVer", OBIT_long, dim, &ver);

  /* Save Seeing */
  dim[0] = dim[1] = 1;
  seeing = -1.0;
  ObitInfoListGetTest(ionCal->info, "seeing",  &type, dim, &seeing);
  ObitInfoListAlwaysPut(myInput, "seeing", type, dim, &seeing);

  /* Restore NField */
  ObitInfoListAlwaysPut(inUV->info, "NField", type, dim, &NField);

  /* Restore do3D on uv data */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "do3D", OBIT_bool, dim, &do3D);

  /* Cleanup*/
  ionCal = ObitIonCalUnref(ionCal);
} /* end doFieldCal */

/*----------------------------------------------------------------------- */
/*  Field based Imaging                                                   */
/*   Input:                                                               */
/*      Stokes    Input Stokes type                                       */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to image                                         */
/*      myClean   CLEAN object                                            */
/*   Output:                                                              */
/*      nfield Number of fields in image mosaic                           */
/*      ncomp  Array of number of components in data (i.e. after any peel)*/
/*             per field,  this array should be deallocated after use     */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doFieldImage (gchar *Stokes, ObitInfoList* myInput, ObitUV* inUV, 
		   ObitDConCleanVis *myClean, olong *nfield, olong **ncomp, 
		   ObitErr* err)
{
  ObitInfoType type;
  olong        prtLv;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat       autoCen, PeelFlux, ftemp;
  gboolean     Fl = FALSE, Tr = TRUE, doRestore, doFlatten, reimage;
  gboolean     doneRecenter=FALSE;
  gchar        *FCParms[] = {  /* Imaging parameters */
    "ionVer", "prtLv", 
    NULL
  };
  gchar *routine = "doFieldImage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inUV));

  /* Don't restore and flatten before done */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(myClean->info, "doRestore", OBIT_bool, dim, &Fl);
  ObitInfoListAlwaysPut(myClean->info, "doFlatten", OBIT_bool, dim, &Fl);
  /* Explicitly do weighting  */
  ObitInfoListAlwaysPut(myClean->info, "doWeight", OBIT_bool, dim, &Tr);

  prtLv = 1;
  ObitInfoListGetTest(myInput, "prtLv",  &type, dim, &prtLv);
  ObitInfoListAlwaysPut(myClean->info, "prtLv", type, dim, &prtLv);
 
  /* Copy control info */
  ObitInfoListCopyList (myInput, inUV->info, FCParms);

  /* Peeling trip level */
  PeelFlux = 1.0e20;
  ObitInfoListGetTest(myInput, "PeelFlux", &type, dim, &PeelFlux);
 
  /* Auto center trip level */
  autoCen = 1.0e20;
  ObitInfoListGetTest(myInput, "autoCen", &type, dim, &autoCen);
  /* Recentering trip level in CLEAN */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1.1 * MIN (autoCen, PeelFlux); /* Fudge a bit due to shallow CLEAN */
  ObitInfoListAlwaysPut (myClean->info, "autoCen", OBIT_float, dim, &ftemp);

  /* Image/Clean */
  ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);

  /* Only recenter once */
  ftemp = 1.0e20;
  ObitInfoListAlwaysPut (myClean->info, "autoCen", OBIT_float, dim, &ftemp);

  /* Need to recenter bright sources? */
  if (((myClean->peakFlux>autoCen) || (myClean->peakFlux> PeelFlux)) && !doneRecenter) {
    /* Compress CC files */
    ObitSkyModelCompressCC (myClean->skyModel, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(myClean->mosaic->info, "restartFlux", 
			  OBIT_float, dim, &autoCen);
    reimage = ObitDConCleanVisReimage (myClean, inUV, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    
    /* Always reImage/Clean if you get here */
    /* Don't need to remake beams  */
    dim[0] = 1;dim[1] = 1;
    ObitInfoListAlwaysPut(myClean->info, "doBeam", OBIT_bool, dim, &Fl);
    
    if (reimage) {      
      Obit_log_error(err, OBIT_InfoErr, 
		     "Redoing image/deconvolution to center strong source on pixel");
    } else {
      Obit_log_error(err, OBIT_InfoErr, 
		     "AutoCenter false alarm - continue CLEANing");
      /* Reuse all prior components - note 0 here means none */
      ftemp  = 0.01*autoCen;
      dim[0] = 1;dim[1] = 1;
      ObitInfoListAlwaysPut (myClean->info, "reuseFlux", OBIT_float, dim, &ftemp);
    }
    ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    
    autoCen = 1.0e20;  /* only once */
    doneRecenter = TRUE;
  } /* End auto center */

  /* Loop peeling sources */
  ObitUVPeelUtilLoop (myInput, inUV, myClean, nfield, ncomp, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);

  /* Restore if requested */
  doRestore = TRUE;
  ObitInfoListGetTest(myInput, "doRestore", &type, dim, &doRestore);
  if (doRestore) { 
    ObitDConCleanRestore((ObitDConClean*)myClean, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    /* Cross restore? */
    if (myClean->nfield>1)
      ObitDConCleanXRestore((ObitDConClean*)myClean, err);
  }

  /* Flatten if requested */
  doFlatten = TRUE;
  ObitInfoListGetTest(myInput, "doFlatten", &type, dim, &doFlatten);
  if (doFlatten) {
    ObitDConCleanFlatten((ObitDConClean*)myClean, err);

    /* Display? */
    if (myClean->display)
      ObitDisplayShow (myClean->display, (Obit*)myClean->mosaic->FullField, NULL, 
		       1, err);
  } else {
    /* Display mosaic? */
    if (myClean->display)
      ObitDisplayShow (myClean->display, (Obit*)myClean->mosaic, NULL, 
		       1, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);

  /* Cleanup */
} /* end doFieldImage */

/*----------------------------------------------------------------------- */
/* Subtract IPol skyModel from inData with possible application of        */
/* calibration.  The subtraction cannot apply the calibration so          */
/* copy/calibrate to a scratch file.                                      */
/* Calibration is turned off on inData                                    */
/*   Input:                                                               */
/*      inData    ObitUV to copy history from                             */
/*      skyModel  Skymodel to subtract                                    */
/*      selFGver  Continuum channel selection FG flag, -1 if none         */
/*      nfield    Number of fields in image mosaic                        */
/*      ncomp     Array of number of components in data to use -1=>none   */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void subIPolModel (ObitUV* outData,  ObitSkyModel *skyModel, olong *selFGver, 
		   olong nfield, olong *ncomp, ObitErr* err)
{
  ObitUV *scrUV = NULL;
  oint otemp;
  olong i, jtemp, nfld, *itemp=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar IStokes[5];
  gchar *include[] = {"AIPS FG", "AIPS NI", NULL};
  gchar *routine = "subIPolModel";

  /* Message */
  Obit_log_error(err, OBIT_InfoErr, "Subtracting IPol model from output uv data");
  ObitErrLog(err); 

  /* unset selection flagging */
  dim[0] = dim[1] = dim[2] = 1;
  otemp = -1;
  ObitInfoListAlwaysPut (outData->info, "flagVer", OBIT_oint, dim, &otemp);
 
    /* Reset Sky Model to use ncomp components */
  nfld = skyModel->mosaic->numberImages;
  itemp = ObitMemAlloc(nfld*sizeof(olong));  /* temp. array */
  dim[0] = nfld;
  for (i=0; i<nfld; i++) {
    if (i<nfield) {
      if (ncomp[i]>=0) {
	itemp[i] = 1;
      } else { /* None */
	itemp[i] = 2;
      }
    } else {  /* In case but should never happen */
      itemp[i] = 0;
    }
  }
  ObitInfoListAlwaysPut(skyModel->info, "BComp", OBIT_long, dim, itemp);
  for (i=0; i<nfld; i++) {
    if (i<nfield) {
      if (ncomp[i]>=0) {
	itemp[i] = ncomp[i];
      } else { /* None */
	itemp[i] = 1;
      }
    } else {  /* In case but should never happen */
      itemp[i] = 0;
    }
  }
  ObitInfoListAlwaysPut(skyModel->info, "EComp", OBIT_long, dim, itemp);
  itemp = ObitMemFree(itemp);  /* Deallocate */
    
  /* No translation in Stokes */ 
  dim[0] = 4;
  sprintf (IStokes, "    ");
  ObitInfoListAlwaysPut (outData->info, "Stokes", OBIT_string, dim, IStokes);
  /* No calibration */
  dim[0] = 1; jtemp = -1;
  ObitInfoListAlwaysPut (outData->info, "doCalib", OBIT_long, dim, &jtemp);
  
  /* Open and close uvdata to set descriptor for scratch file */
  ObitUVOpen (outData, OBIT_IO_ReadCal, err);
  ObitUVClose (outData, err);

  /* Copy to scratch with calibration */
  scrUV = newObitUVScratch (outData, err);
  scrUV = ObitUVCopy (outData, scrUV, err);
  if (err->error) goto cleanup;

  /* Any selection flagging table gets lost in the shuffle unless we copy it 
     also need NI table on scrUV */
  ObitUVCopyTables (outData, scrUV, NULL, include, err);
  if (err->error) goto cleanup;
  
  /* No more calibration */
  dim[0] = 1; jtemp = -1;
  ObitInfoListAlwaysPut (scrUV->info, "doCalib", OBIT_long, dim, &jtemp);
  ObitInfoListAlwaysPut (scrUV->info, "doBand",  OBIT_long, dim, &jtemp);

    /* Subtract */
  ObitSkyModelSubUV (skyModel, scrUV, outData, err);
  if (err->error) goto cleanup;

  /* Make sure something happened */
  Obit_return_if_fail((outData->myDesc->nvis>1), err, 
			"%s: NO Data written to output", routine);

  /* Copy any selection flagging/NI table back */
  ObitUVCopyTables (scrUV, outData, NULL, include, err);
  if (err->error) goto cleanup;
  
  /* reset selection flagging */
  dim[0] = dim[1] = dim[2] = 1;
  otemp = (oint)(*selFGver);
  ObitInfoListAlwaysPut (outData->info, "flagVer", OBIT_oint, dim, &otemp);

  /* No more calibration */
  dim[0] = 1; jtemp = -1;
  ObitInfoListAlwaysPut (outData->info, "doCalib", OBIT_long, dim, &jtemp);
  ObitInfoListAlwaysPut (outData->info, "doBand",  OBIT_long, dim, &jtemp);

  /* Cleanup */
 cleanup:
  scrUV = ObitUVUnref(scrUV);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  
} /* end subIPolModel */

/*----------------------------------------------------------------------- */
/*  Write History for IonImage                                            */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      Stoke     Stokes's parameter imaged I, Q, U, V                    */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outImage  ObitImage to write history to                           */
/*      outData   If non NULL, output UV data to write history to         */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void IonImageHistory (gchar *Source, gchar Stoke, ObitInfoList* myInput, 
		    ObitUV* inData, ObitImage* outImage, ObitUV* outData, 
		    ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "souCode", "Qual", 
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "UVRange",  "timeRange",  "Robust",  "UVTaper",  "WtBox", "WtFunc", 
    "BIF", "EIF", "BChan", "EChan",  "chInc", "chAvg", "BLFact", "BLFOV",  "BLchAvg",
    "doCalSelect",  "doCalib",  "gainUse", "doBand ",  "BPVer",  "flagVer", 
    "doPol",  "doFull", "do3D", "Catalog",  "OutlierDist",  "OutlierFlux",  "OutlierSI",
    "OutlierSize",  "CLEANBox",  "Gain",  "minFlux",  "Niter",  "minPatch",
    "FOV", "xCells", "yCells", "nx", "ny", "RAShift", "DecShift", "doRestore",
    "Beam",  "CCFilter",  "maxPixel", "autoWindow", "subA",
    "doFCal", "ionVer", "UpdateInt", "solInt", "nZern", "FitDist", "MaxDist", "MinPeak", 
    "MaxWt", "MaxQual", "MaxRMS", "MinRat", "FCStrong", "FCNiter", "FCGain", "FCminFlux", 
    "seeing", "autoCen","PBCor", "antSize", "Alpha",
    "PeelFlux", "PeelRefAnt", "PeelSNRMin", "PeelSolInt", "PeelNiter", "PeelMinFlux",
    "PeelAvgPol", "PeelAvgIF", "PeelType", "PeelMode",
    "nThreads",
    NULL};
  gchar *routine = "IonImageHistory";

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
 
  /* Do UV history -f outData given  */
  if (outData==NULL) return;
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

  /* Write source and poln */
  if (Stoke==' ') Stoke = 'I';
  g_snprintf (hicard, 80, "%s Source = '%s', Stokes= '%c'", pgmName, Source, Stoke);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end IonImageHistory  */

#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif
/*----------------------------------------------------------------------- */
/*  Copy or baseline dependent time average data                          */
/*  If BLFact>1.00 then use baseline dependent time averaging, else       */
/*  just a straight copy from inData to outData                           */
/*  Uses minimum of solPint or solAInt as the maximum time average        */
/*  or if this is zero then 1 min.                                        */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList use:                       */
/*       "BLFact"  OBIT_float  (1,1,1) Maximum time smearing factor       */
/*       "BLFOV"   OBIT_float  (1,1,1) Field of view (radius, deg)        */
/*                                     Default FOV or 0.5*lambda/25.0 m   */
/*       "BLchAvg" OBIT_bool   (1,1,1) Also average channels? [def FALSE] */
/*       "solPInt" OBIT_float  (1,1,1) Phase self-cal soln. interval (min)*/
/*       "solAInt" OBIT_float  (1,1,1) Amp self-cal soln. interval (min)  */
/*      inData    ObitUV to copy data from                                */
/*      outData   Output UV data to write                                 */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void BLAvg (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
	    ObitErr* err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NumChAvg=1;
  odouble Freq;
  gboolean BLchAvg=FALSE;
  ofloat BLFact=0.0, FOV=0.0, solPInt=0.0, solAInt=0.0, maxInt;
  gchar *blank = "    ";
  gchar *routine = "BLAvg";

  /* What to do? */
  ObitInfoListGetTest(myInput, "BLFact", &type, dim, &BLFact);
  if (BLFact>1.00) { /* Average */
    /* Set parameters */
    ObitInfoListGetTest(myInput, "BLchAvg",   &type, dim, &BLchAvg);
    ObitInfoListGetTest(myInput, "BLFOV",   &type, dim, &FOV);
    if (FOV<=0.0) ObitInfoListGetTest(myInput, "FOV",   &type, dim, &FOV);
    ObitInfoListGetTest(myInput, "solInt",     &type, dim, &solPInt);
    solPInt /= 4.0;
    ObitInfoListGetTest(myInput, "PeelSolInt", &type, dim, &solAInt);
    if (solAInt<=0.0) solAInt = solPInt;
    if (solPInt<=0.0) solPInt = solAInt;
    maxInt = MIN (solPInt, solAInt);
    if (maxInt<=0.0) maxInt = 1.0;  /* Default 1 min */
    /* Default FOV 0.5 lambda/diameter */
    if (FOV<=0.0) {
      Freq = inData->myDesc->crval[inData->myDesc->jlocf];
      FOV = 0.5 * Freq * VELIGHT / 25.0;;
    }

    /* Average channels? */
    if (BLchAvg) {
      NumChAvg = ObitUVUtilNchAvg(inData, BLFact, FOV, err);
      if (err->error) Obit_traceback_msg (err, routine, inData->name);
      NumChAvg = MAX (1, NumChAvg);
      Obit_log_error(err, OBIT_InfoErr, 
		     "Averaging %d channels", NumChAvg);
    }
    
    dim[0] = dim[1] = dim[2] = dim[3] = 1;
    ObitInfoListAlwaysPut (inData->info, "FOV",      OBIT_float, dim, &FOV);
    ObitInfoListAlwaysPut (inData->info, "maxInt",   OBIT_float, dim, &maxInt);
    ObitInfoListAlwaysPut (inData->info, "maxFact",  OBIT_float, dim, &BLFact);
    ObitInfoListAlwaysPut (inData->info, "NumChAvg", OBIT_long,  dim, &NumChAvg);
    outData = ObitUVUtilBlAvgTF(inData, FALSE, outData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* (Re)set stuff 
       No source selection */
    dim[0] = 4; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (inData->info, "Sources", OBIT_string,   dim, blank);
    /* Reset FOV to previous value */
    ObitInfoListGetTest(myInput, "FOV",   &type, dim, &FOV);
    dim[0] = 1; dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut (inData->info, "FOV", OBIT_float,  dim, &FOV);
    
  } else { /* Straight copy */
    ObitUVCopy (inData, outData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  }
} /* end BLAvg */
