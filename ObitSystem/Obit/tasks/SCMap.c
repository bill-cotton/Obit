/* $Id$  */
/* Obit task to image/CLEAN/selfcalibrate a uv data set             */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2014                                          */
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
#include "ObitUVUtil.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* SCMapIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void SCMapOut (ObitInfoList* outList, ObitErr *err);
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
void setOutputData (gchar *Source, ObitInfoList *myInput, 
		    ObitUV* inData, ObitImage **outImage, ObitErr *err);

/* Loop over sources */
void doSources (ObitInfoList* myInput, ObitUV* inData, ObitErr* err);

/* Loop over Channels/Poln */
void doChanPoln (gchar *Source, ObitInfoList* myInput, ObitUV* inData, 
		 ObitErr* err);

/* Image/self cal loop */
void doImage (ObitInfoList* myInput, ObitUV* inData, 
	      ObitDConCleanVis *myClean, gboolean autoWindow, ObitErr* err);
/* Write history */
void SCMapHistory (gchar *Source, ObitInfoList* myInput, 
		    ObitUV* inData, ObitData* outData, ObitErr* err);

/* Baseline dependent time averaging */
void BLAvg (ObitInfoList* myInput, ObitUV* inData, ObitUV* outData, 
	    ObitErr* err);

/* Program globals */
gchar *pgmName = "SCMap";       /* Program name */
gchar *infile  = "SCMap.in" ;   /* File with program inputs */
gchar *outfile = "SCMap.out";   /* File to contain program outputs */
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
/*   Obit task to image a uv data set                                     */
/*----------------------------------------------------------------------- */
{
  oint         ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitUV       *inData = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = SCMapIn (argc, argv, err);
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
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto done;

  /* Get input uvdata */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto done;

  /* Process */
  doSources (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto done;

  /* Set up output */
 done:
  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  
  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* SCMapIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "SCMapIn";

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
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

  return list;
} /* end SCMapIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: SCMap -input file -output ofile [args]\n");
    fprintf(stderr, "SCMap Obit task to image/CLEAN data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def SCMap.in\n");
    fprintf(stderr, "  -output output result file, def SCMap.out\n");
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
/*     outClass  Str [6]    output AIPS image class  [IClean]             */
/*     outSeq    Int        output AIPS image sequence no  [new]          */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     out2File  Str [?]    output FITS uv file name [def "Image.fits"    */
/*     out2Name  Str [12]   output AIPS uv name  [no def]                 */
/*     out2Class Str [6]    output AIPS uv class  [SCmap]                 */
/*     out2Seq   Int        output AIPS uv  sequence no  [new]            */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     Sources   Str (16,1) Sources selected, blank = all                 */
/*     Stokes    Str (4)    Stokes parameter to image, def=I              */
/*     FOV       Flt (1)    Field of view in deg , NO DEFAULT (0.0)       */
/*     UVRange   Flt (2)    Range n uv plane in klambda, def=all          */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     Robust    Flt (1)    Briggs robust factor (AIPS version), def=0.0  */
/*     UVTaper   Flt (2)    Taper in uv plane in klambda in u, v, def=all */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
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
/*     Reuse     Flt (1)    Restart multiple of field 1 RMS, def = 10.0   */
/*     minPatch  Int (1)    Clean Min. BEAM half-width, def=200           */
/*     Beam      Flt (3)    Clean beam maj, min, PA (", ", deg) def=fit   */
/*     CCFilter  Flt (2)    CC filter, [min. sum, radius in pix.], def=no */
/*     maxPixel  Int (1)    Max. pixels in inner cycle, def=50000         */
/*     maxSCLoop Int (1)    Maximum number of selfcal loops, def=0        */
/*     subA      Int (1)    Subarray, def=1                               */
/*     minFluxPSC Flt(1)    min peak flux for phase selfcal               */
/*     minFluxASC Flt(1)    min peak flux for A&P selfcal                 */
/*     Alpha      Flt(1)    default spectral index (0)                    */
/*     dispURL    Str(48)   Display derver URL                            */
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
  strTemp = "SCMap.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "SCMapName";
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
  strTemp = "SCMapOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Output AIPS Image file name */
  strTemp = "SCMapOut";
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

  /* Restart multiple of field 1 RMS, def = 10.0 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 10.0; 
  ObitInfoListPut (out, "Reuse", OBIT_float, dim, &ftemp, err);
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

  /* Modeling method: 'DFT','GRID','    ', def = '    ' (chose fastest) */
  strTemp = "    ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Cmethod ", OBIT_string, dim, strTemp, err);
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

 /* min peak flux for phase selfcal  */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1.0e20; 
  ObitInfoListPut (out, "minFluxPSC", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

 /* minFluxPSC min peak flux for A&P selfcal  */
  dim[0] = 1;dim[1] = 1;
  ftemp = 1.0e20; 
  ObitInfoListPut (out, "minFluxASC", OBIT_float, dim, &ftemp, err);
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
  gchar *strTemp;
  gboolean *booTemp, btemp;
  olong itemp;
  ObitSkyModelMode modelMode;
  gchar *routine = "digestInputs";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Convert test Cmethod to enum  Mode */
  ObitInfoListGetP (myInput, "Cmethod", &type, dim, (gpointer)&strTemp);
  if (!strncmp (strTemp, "GRID", 4)) modelMode = OBIT_SkyModel_Grid;
  else if (!strncmp (strTemp, "DFT", 3)) modelMode = OBIT_SkyModel_DFT;
  else modelMode = OBIT_SkyModel_Fastest;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (myInput, "Mode", OBIT_long, dim, &modelMode);

  /* noScrat - no scratch files for AIPS disks */
  ObitAIPSSetnoScrat(myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, "task Input");

  /* Default NField is 1 */
  itemp = 1;  type = OBIT_long; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListGet(myInput, "NField", &type, dim, &itemp, err);

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
    
    /* define object */
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
/*  Create output uv data                                                 */
/*  Sets AIPS like Name,CLASS,seq info even for FITS files                */
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
  gchar     *Type, *strTemp, out2File[129], *out2Name, *out2F;
  gchar     Aname[13], Aclass[7], *Atype = "UV";
  olong      nvis;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean  exist;
  gchar     tname[129];
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
    strncpy (Aname, tname, 13);
    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "out2Class", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "SCMap ", 7);
    }
    /* Default for blank */
    if (!strncmp (Aclass, "    ", 4)) strncpy (Aclass, "SCMap ", 7);

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
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, out2Name */
    ObitInfoListGetP (myInput, "out2File", &type, dim, (gpointer)&out2F);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = out2F[i]; tname[i] = 0;
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (out2File, 128, "%s", tname);
    else g_snprintf (out2File, 128, "%s%s", Source, tname);
    ObitTrimTrail(out2File);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "out2Disk", &type, dim, &disk, err);
    if (disk<=0) /* defaults to outDisk */
      ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* define object */
    nvis = 1000;
    ObitInfoListGetTest(inData->info, "nVisPIO", &type, dim, &nvis);
    ObitUVSetFITS (outUV, nvis, disk, out2File, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outUV);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
		   pgmName, Type);
    return outUV;
  }
  
  return outUV;
} /* end setOutputUV */

/*----------------------------------------------------------------------- */
/*  Set output info on uvdata, create output image                        */
/*  Sets AIPS like Name,CLASS,seq info even for FITS files                */
/*  One output image per requested Stokes                                 */
/*   Input:                                                               */
/*      Source    Source name                                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to accept parameters defining output             */
/*   Output:                                                              */
/*      outImage  Output image depending on Stokes request                */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void setOutputData (gchar *Source, ObitInfoList *myInput, 
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
  gchar     tname[129];
  gchar     *routine = "setOutputData";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Create basic output Image Object */
  g_snprintf (tname, 100, "output Image IPol");
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
      strncpy (Aclass, "xClean", 7);
    }
    /* Default for blank */
    if (!strncmp (Aclass, "    ", 4)) strncpy (Aclass, "xClean", 7);
    Aclass[0] = 'I';

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
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS output */
    /* Generate output name from Source, outName */
    ObitInfoListGetP (myInput, "outFile", &type, dim, (gpointer)&outF);
    n = MIN (128, dim[0]);
    for (i=0; i<n; i++) tname[i] = outF[i]; tname[i] = 0;
    /* Something in source name? */
    if ((Source[0]==' ') || (Source[0]==0)) 
      g_snprintf (outFile, 128, "I%s", tname);
    else g_snprintf (outFile, 128, "%sI%s", Source, tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
	   
    IOType = OBIT_IO_FITS;  /* Save file type */

    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    
    /* define object */
    ObitImageSetFITS ((*outImage), OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_msg (err, routine, "myInput");
    
    /* Make up AIPS-like name, class...  */
    if (strncmp (Source, "    ", 4))
      strncpy (Aname, Source, 13);
    else
      strncpy (Aname, outFile, 13);
    strncpy (Aclass, "IMap", 7);
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
  dim[0] = strlen (Aname);
  ObitInfoListAlwaysPut(inData->info, "imName", OBIT_string, dim, Aname);
  dim[0] = strlen (Aclass);
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
  olong         i, isource, nsource;
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
  ObitDConCleanVis *myClean=NULL;
  ObitUV       *outData = NULL;
  ObitImage    *outField=NULL;
  ObitImage    *outImage=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        chInc, BChan, EChan, nchan, oldSN, win[4];
  gboolean     doFlat, btemp, autoWindow, Tr=TRUE, do3D;
  olong        inver, outver, *CLEANBox, plane[5] = {0,1,1,1,1};
  gchar        Stokes[5], *CCType = "AIPS CC";
  ofloat       modelFlux=0.0;
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "UVRange", "timeRange", "UVTape",
    "BIF", "EIF", "subA", "FreqID", "souCode", "Qual", 
    "doCalSelect", "doCalib", "gainUse", "doBand", "BPVer", "flagVer", 
    "doPol", "PDVer", "Mode",
    NULL
  };
  gchar        *tmpParms[] = {  /* Imaging, weighting parameters */
    "doFull", "do3D", "FOV", "PBCor", "antSize", 
    "Catalog", "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
    "Robust", "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "WtPower",
    "MaxBaseline", "MinBaseline", "rotate", "Beam", "minFlux",
    "NField", "xCells", "yCells","nx", "ny", "RAShift", "DecShift",
    "nxBeam", "nyBeam", "Alpha", "doCalSelect", 
    NULL
  };
  gchar        *tmpName[] = {  /* Names to use for Image mosaic files */
    "imFileType", "imName", "imClass", "imDisk", "imSeq", "Sources",
    NULL
  };
  gchar        *CLEANParms[] = {  /* Clean parameters */
    "CLEANBox", "Gain", "minFlux", "Niter", "minPatch", "Beam", 
    "Mode", "CCFilter", "maxPixel", "dispURL",
    NULL
  };
  gchar *routine = "doChanPoln";

   /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));

  /* Parameters used here */
  /* Average everything selected = Continuum  */
  nchan = inData->myDesc->inaxes[inData->myDesc->jlocf];
  BChan = 1;
  ObitInfoListGetTest(myInput, "BChan",  &type, dim, &BChan);
  BChan = MAX (1, BChan);
  EChan = nchan;
  ObitInfoListGetTest(myInput, "EChan",  &type, dim, &EChan);
  EChan = MIN (EChan, nchan);
  if (EChan<=0) EChan = nchan;
  /* Update inputs */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(myInput, "EChan",  OBIT_long, dim, &EChan);
  chInc = EChan - BChan + 1;
  chInc = MAX (1, chInc);
  doFlat = TRUE;
  ObitInfoListGetTest(myInput, "doFlatten", &type, dim, &doFlat);
  autoWindow = FALSE;
  ObitInfoListGetTest(myInput, "autoWindow", &type, dim, &autoWindow);
  do3D = TRUE;
  ObitInfoListGetTest(myInput, "do3D", &type, dim, &do3D);

  /* Restarting soln? */
  oldSN = -1;
  ObitInfoListGetTest(myInput, "oldSN", &type, dim, &oldSN);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Create output file for data */
  outData = setOutputUV (Source, myInput, inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* set data selection */
  dim[0] = 1;
  ObitInfoListAlwaysPut (inData->info, "BChan", OBIT_long, dim, &BChan);
  ObitInfoListAlwaysPut (inData->info, "EChan", OBIT_long, dim, &EChan);
  
  /* Calibrate/edit/copy data as correlator data to output file */
  dim[0] = 4;
  sprintf (Stokes, "    ");
  ObitInfoListAlwaysPut (inData->info, "Stokes", OBIT_string, dim, Stokes);
  /* set selected Source  */
  dim[0] = 16; dim[1] = 1;
  ObitInfoListAlwaysPut (inData->info, "Sources", OBIT_string, dim, Source);
  /* Copy data unless restarting */
  if (oldSN<0) {
    BLAvg (myInput, inData, outData, err);
    /* Old outData = ObitUVCopy (inData, outData, err);*/
  } else {
    /* Ensure inData fully instantiated and OK */
    ObitUVFullInstantiate (outData, TRUE, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  ObitInfoListCopyList (myInput, outData->info, tmpParms);
  ObitInfoListCopyList (inData->info, outData->info, tmpName);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
  /* Image/calibrate Ipol  */
  sprintf (Stokes, "I   ");
  dim[0] = 4;
  ObitInfoListAlwaysPut (outData->info, "Stokes", OBIT_string, dim, Stokes);
  dim[0] = 1;
  ObitInfoListAlwaysPut (outData->info, "doCalSelect", OBIT_bool, dim, &Tr);
  
  /* Define output image */
  setOutputData (Source, myInput, outData, &outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
  /* Make CleanVis */
  myClean = ObitDConCleanVisCreate("Clean Object", outData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Get input parameters from myInput, copy to myClean */
  ObitInfoListCopyList (myInput, myClean->info, CLEANParms);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);	
    
  /* Set CLEAN windows  */
  ObitDConCleanVisDefWindow((ObitDConClean*)myClean, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);

  /* Check if need to add a small window at the origin? */    
  ObitInfoListGetTest(myInput, "modelFlux",  &type, dim, &modelFlux);
  ObitInfoListGetP(myInput, "CLEANBox",  &type, dim, (gpointer)&CLEANBox);
  if ((modelFlux>0.0) && CLEANBox && (CLEANBox[0]==0)) {
    /* `yes add 5 pixel circle at center - change default window */
    win[0] = 5;
    win[1] = myClean->mosaic->nx[0]/2 + 1;
    win[2] = myClean->mosaic->ny[0]/2 + 1;
    win[3] = 0;
    ObitDConCleanWindowUpdate (myClean->window, 1, 1, OBIT_DConCleanWindow_round,
			       win, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  }

  if (myClean->mosaic->numberImages<=1) doFlat = FALSE;  /* Don't flatten 1 */
  /* Create output image(s) */
  if (doFlat) 
    outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
  else
    outField = ObitImageMosaicGetImage (myClean->mosaic, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  ObitImageUtilMakeCube (outField->myDesc, inData->myIO->myDesc, 
			 outImage->myDesc,  Stokes, BChan, EChan, chInc, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  ObitImageFullInstantiate (outImage, FALSE, err);
  outField = ObitImageUnref(outField);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
  /* Set Stokes on SkyModel */
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (myClean->skyModel->info, "Stokes", OBIT_string, dim, Stokes);
  
  /* Automatic windowing?  */
  btemp = autoWindow;
  
  /* Do actual processing */
  doImage (myInput, outData, myClean, btemp, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Copy result to output */
  plane[0] = 1;
  if (doFlat) 
    outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
  else { /* Copy the first image */
    outField = ObitImageMosaicGetImage (myClean->mosaic, 0, err);
    /* and copy the CC table */
    inver = 1;
    outver = plane[0];
    ObitDataCopyTable ((ObitData*)outField, (ObitData*)outImage,
		       CCType, &inver, &outver, err);
  }
  ObitImageUtilInsertPlane (outField, outImage, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  
  /* For 2D imaging with flatten copy CC Table */
  if (!do3D && doFlat) {
    inver   = 1;
    outver  = plane[0];
    ObitDataCopyTable ((ObitData*)outField, (ObitData*)outImage,
		       CCType, &inver, &outver, err);
    if (err->error) Obit_traceback_msg (err, routine, outField->name);
  }
  outField = ObitImageUnref(outField);
  
  /* Make sure image created */
  Obit_return_if_fail((outImage!=NULL), err, "%s: No image generated", routine);

  /* Do history */
  SCMapHistory (Source, myInput, inData, (ObitData*)outImage, err);
  SCMapHistory (Source, myInput, inData, (ObitData*)outData, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  outImage  = ObitUnref(outImage);
  
  
  /* Cleanup */
  if (myClean) {
    /* Leave facet images if not myClean->mosaic->doFull and myClean->mosaic->numberImages >1 */
    if (!(!myClean->mosaic->doFull && (myClean->mosaic->numberImages>1)))
      ObitImageMosaicZapImage (myClean->mosaic, -1, err); /* Delete mosaic members */
    if (doFlat) {  /* Delete flattened as well if not output */
      outField = ObitImageMosaicGetFullImage (myClean->mosaic, err);
      outField = ObitImageZap(outField, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    }
    myClean  = ObitDConCleanVisUnref(myClean);
  }
  outData  = ObitUVUnref(outData);

}  /* end doChanPoln */

/*----------------------------------------------------------------------- */
/*  Imaging/Deconvolution self calibration loop                           */
/*   Input:                                                               */
/*      myInput    Input parameters on InfoList                           */
/*      inUV       ObitUV to image                                        */
/*      myClean    CLEAN object                                           */
/*      autoWindow If True, do autoWindow at end of each self cal cycle   */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doImage (ObitInfoList* myInput, ObitUV* inUV, 
	      ObitDConCleanVis *myClean, gboolean autoWindow, ObitErr* err)
{
  ObitUVSelfCal *selfCal = NULL;
  ObitUV       *scrUV = NULL;
  ObitImage    *outImage=NULL;
  ObitInfoType type;
  olong         maxPSCLoop, maxASCLoop, SCLoop, oldSN, itemp, jtemp;
  ofloat       minFluxPSC, minFluxASC, modelFlux, reuse, ftemp, autoCen;
  ofloat       solInt, minFlux=0.0;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean     Fl = FALSE, Tr = TRUE, init=TRUE, doRestore, doFlatten, doSC;
  gboolean     btemp, noNeg, noSCNeed, reimage, didSC=FALSE, imgOK=FALSE, converged = FALSE;
  gchar        Stokes[5], soltyp[5], solmod[5], stemp[5];
  gchar        *SCParms[] = {  /* Self parameters */
    "maxSCLoop", "minFluxPSC", "minFluxASC", "refAnt", 
    "WtUV", "avgPol", "avgIF", "noNeg", "doMGM", "minSNR", 
    "minNo", "doSmoo", "prtLv", "modelFlux", "modelPos", "modelParm",
    "dispURL", 
    NULL
  };
  const ObitDConCleanClassInfo *myCleanClass = (ObitDConCleanClassInfo*)myClean->ClassInfo;
  gchar *routine = "doImage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inUV));

  /* Parameters used here*/
  ObitInfoListGet(myInput, "maxPSCLoop",  &type, dim, &maxPSCLoop,  err);
  ObitInfoListGet(myInput, "maxASCLoop",  &type, dim, &maxASCLoop,  err);
  ObitInfoListGet(myInput, "minFluxPSC", &type, dim, &minFluxPSC, err);
  ObitInfoListGet(myInput, "minFluxASC", &type, dim, &minFluxASC, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  reuse = 10.0;
  ObitInfoListGetTest(myInput, "Reuse",&type, dim, &reuse);
  autoCen = 1.0e20;
  ObitInfoListGetTest(myInput, "autoCen", &type, dim, &autoCen);
  oldSN = -1;
  ObitInfoListGetTest(myInput, "oldSN", &type, dim, &oldSN);
  minFlux = 0.0;
  ObitInfoListGetTest(myInput, "minFlux", &type, dim, &minFlux);
  noNeg = TRUE;
  ObitInfoListGetTest(myInput, "noNeg", &type, dim, &noNeg);

  /* Get Stokes being imaged */
  strncpy (Stokes, "F   ", 4); 
  ObitInfoListGetTest (inUV->info, "Stokes", &type, dim, Stokes);

  /* Only do self cal for Stokes I (or F) */
  if ((Stokes[0]!='I') && (Stokes[0]!='F') && ((Stokes[0]!=' '))) {
    maxPSCLoop  = 0;
    maxASCLoop  = 0;
    minFluxPSC = 1.0e20;
    minFluxASC = 1.0e20;
  }
  
  /* Don't restore and flatten before done */
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(myClean->info, "doRestore", OBIT_bool, dim, &Fl);
  ObitInfoListAlwaysPut(myClean->info, "doFlatten", OBIT_bool, dim, &Fl);
  /* Explicitly do weighting - needed to apply SC solutions */
  ObitInfoListAlwaysPut(myClean->info, "doWeight", OBIT_bool, dim, &Tr);

  /* No reuse of components the first time */
  ftemp = -1.0;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (myClean->info, "reuseFlux", OBIT_float, dim, &ftemp);
  /* Recentering trip level in CLEAN */
  ObitInfoListAlwaysPut (myClean->info, "autoCen", OBIT_float, dim, &autoCen);

  /* Create selfCal */
  doSC = ((maxPSCLoop>0) || (maxASCLoop>0));
  if (doSC>0) {
    selfCal = ObitUVSelfCalCreate ("SelfCal", myClean->skyModel);
    /* Copy control info */
    ObitInfoListCopyList (myInput, selfCal->info, SCParms);
    
    /* Set vis vs baseline histogram */
    ObitUVSelfCalFluxHist(selfCal, inUV, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);

    /* Need initial model calibration? oldSN overrides */
    if (oldSN<0) {
      modelFlux = 0.0;
      ObitInfoListGetTest(myInput, "modelFlux", &type, dim, &modelFlux);
      if (modelFlux>0.0) {
	solInt = 0.0; dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListGetTest(myInput, "solPInt", &type, dim, &solInt);
	ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float, dim, &solInt);
	soltyp[0] = soltyp[1] = soltyp[2] = soltyp[3] = ' '; soltyp[4] = 0;
	ObitInfoListGetTest(myInput, "solPType", &type, dim, soltyp);
	dim[0] = 4;
	ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string, dim, soltyp);
	solmod[0] = solmod[1] = solmod[2] = solmod[3] = ' '; solmod[4] = 0;
	ObitInfoListGetTest(myInput, "solPMode", &type, dim, solmod);
	ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
	ObitUVSelfCalModel (selfCal, inUV, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      }
    }

    /* Restarting with existing SN table? */
    if (oldSN>=0) {
      dim[0] = 1; dim[1] = 1;
      ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
      itemp = 2;
      ObitInfoListAlwaysPut(inUV->info, "doCalib", OBIT_long, dim, &itemp);
      itemp = oldSN;
      ObitInfoListAlwaysPut(inUV->info, "gainUse", OBIT_long, dim, &itemp);
      Obit_log_error(err, OBIT_InfoErr, 
		     "Restarting selfcal with old SN %d (0=highest)", oldSN);
    }
  }
  
  /******** Phase Self cal loop  or no self cal ********/
  if ((maxPSCLoop>0) || ((maxPSCLoop<=0)&&(maxASCLoop<=0))) {
    init = TRUE;
    converged = FALSE;
    didSC = FALSE;  /* Until proven otherwise */
    for (SCLoop=0; SCLoop<=maxPSCLoop; SCLoop++) {
      
      /* Set specific selfcal parameters */   
      if (doSC) {
	solInt = 0.0; dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListGetTest(myInput, "solPInt", &type, dim, &solInt);
	ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float, dim, &solInt);
	soltyp[0] = soltyp[1] = soltyp[2] = soltyp[3] = ' '; soltyp[4] = 0;
	ObitInfoListGetTest(myInput, "solPType", &type, dim, soltyp);
	dim[0] = 4;
	ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string, dim, soltyp);
	solmod[0] = solmod[1] = solmod[2] = solmod[3] = ' '; solmod[4] = 0;
	ObitInfoListGetTest(myInput, "solPMode", &type, dim, solmod);
	ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
      }
    
      /* Image/Clean */
      ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      imgOK = TRUE; 
     
      /* Only recenter once */
      ftemp = 1.0e20;
      dim[0] = 1;
      ObitInfoListAlwaysPut (myClean->info, "autoCen", OBIT_float, dim, &ftemp);
      
      /* Need to recenter bright sources? */
      if (myClean->peakFlux>autoCen) {
	/* Compress CC files */
	ObitSkyModelCompressCC (myClean->skyModel, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	
	dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut(myClean->mosaic->info, "restartFlux", 
			      OBIT_float, dim, &autoCen);
	reimage = ObitDConCleanVisReimage (myClean, inUV, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	
	/* reImage/Clean */
	if (reimage) {      
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Redoing image/deconvolution to center strong source on pixel");
	  ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
	  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	}
	
	autoCen = 1.0e20;  /* only once */
      }
      
      /* Convergence test */
      if (doSC) {
	if (converged || (SCLoop>=maxPSCLoop)) break;
	Obit_log_error(err, OBIT_InfoErr, " ******  Phase Self Calibration number %d", SCLoop+1);
	
	/* Pass peak found in deconvolution to self cal */
	dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut(selfCal->info, "peakFlux", OBIT_float, dim, &myClean->peakFlux);
	
	/* Reset minFlux disturbed by Clean */
	dim[0] = dim[1] = dim[2] = 1;
	ftemp = 0.0;
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "noNeg", OBIT_bool, dim, &noNeg);
	
	/* Do self cal */
	didSC = TRUE;
	imgOK = FALSE;  /* Need new image */
	converged = ObitUVSelfCalSelfCal (selfCal, inUV, init, &noSCNeed, 
					  myClean->window, err);
	if (err->error) Obit_traceback_msg (err, routine, selfCal->name);
	if (noSCNeed) didSC = FALSE;
	if (converged || noSCNeed)  break;
	init = FALSE;
	
	/* reset flux limit for next Clean to 1 sigma */
	dim[0] = 1;dim[1] = 1;
	ObitInfoListAlwaysPut (myClean->info, "minFlux", OBIT_float, dim, &selfCal->RMSFld1);
	btemp = FALSE;
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "noNeg", OBIT_bool, dim, &btemp);
	
	/* Possibly reuse some of CLEAN model to start next time */
	if (reuse>0.0) {
	  ftemp = reuse*selfCal->RMSFld1;
	  dim[0] = 1;dim[1] = 1;
	  ObitInfoListAlwaysPut (myClean->info, "reuseFlux", OBIT_float, dim, &ftemp);
	} /* end set reuse level */
	
	/* Auto Windowing? */
	if (autoWindow) {
	  myCleanClass->ObitDConCleanAutoWindow ((ObitDConClean*)myClean, 
						 myClean->currentFields, NULL, err);
	  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	}
	
      } /* end if self cal */
      if (noSCNeed) break;
    
    }  /* End Self cal loop */
  } /* End Phase self cal */

  /* If doing two sets of calibration apply previous to data */
  if ((didSC) && ((maxASCLoop>0) && (myClean->peakFlux>minFluxASC))) {
    /* Message */
    Obit_log_error(err, OBIT_InfoErr, "Applying previous calibration to output uv data");
    ObitErrLog(err); 
    /* Apply calibration */
    dim[0] = 1; jtemp = 2;
    ObitInfoListAlwaysPut (inUV->info, "doCalib", OBIT_long, dim, &jtemp);
    jtemp = 0;
    ObitInfoListAlwaysPut (inUV->info, "gainUse",  OBIT_long, dim, &jtemp);
    /* No translation in Stokes */ 
    dim[0] = 4;
    stemp[0] = stemp[1] = stemp[2] = stemp[3] = ' ';  stemp[4] = 0;
    ObitInfoListAlwaysPut (inUV->info, "Stokes", OBIT_string, dim, stemp);
    
    /* Copy to scratch with calibration */
    scrUV = newObitUVScratch (inUV, err);
    scrUV = ObitUVCopy (inUV, scrUV, err);
    /* And then back */
    inUV = ObitUVCopy (scrUV, inUV, err);
    scrUV = ObitUVUnref(scrUV);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
    
    /* No more calibration for now */
    dim[0] = 1; jtemp = -1;
    ObitInfoListAlwaysPut (inUV->info, "doCalib", OBIT_long, dim, &jtemp);

    /* Set Stokes Desired */
    dim[0] = 4;
    ObitInfoListAlwaysPut (inUV->info, "Stokes", OBIT_string, dim, Stokes);
    
  } /* end apply prior calibration to data */

  /******** Amp&PPhase Self cal loop  ********/
  if (maxASCLoop>0) {
    init = TRUE;
    converged = FALSE;
    for (SCLoop=0; SCLoop<=maxASCLoop; SCLoop++) {
      
      /* Set specific selfcal parameters */   
      if (doSC) {
	solInt = 0.0; dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListGetTest(myInput, "solAInt", &type, dim, &solInt);
	ObitInfoListAlwaysPut(selfCal->info, "solInt", OBIT_float, dim, &solInt);
	soltyp[0] = soltyp[1] = soltyp[2] = soltyp[3] = ' '; soltyp[4] = 0;
	ObitInfoListGetTest(myInput, "solAType", &type, dim, soltyp);
	dim[0] = 4;
	ObitInfoListAlwaysPut(selfCal->info, "solType", OBIT_string, dim, soltyp);
	solmod[0] = solmod[1] = solmod[2] = solmod[3] = ' '; solmod[4] = 0;
	ObitInfoListGetTest(myInput, "solAMode", &type, dim, solmod);
	ObitInfoListAlwaysPut(selfCal->info, "solMode", OBIT_string, dim, solmod);
      }
    
      /* Image/Clean */
      if (!imgOK) 
	ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      imgOK = FALSE;
    
      /* Only recenter once */
      ftemp = 1.0e20;
      dim[0] = 1;
      ObitInfoListAlwaysPut (myClean->info, "autoCen", OBIT_float, dim, &ftemp);
      
      /* Need to recenter bright sources? */
      if (myClean->peakFlux>autoCen) {
	/* Compress CC files */
	ObitSkyModelCompressCC (myClean->skyModel, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	
	dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut(myClean->mosaic->info, "restartFlux", 
			      OBIT_float, dim, &autoCen);
	reimage = ObitDConCleanVisReimage (myClean, inUV, err);
	if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	
	/* reImage/Clean */
	if (reimage) {      
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Redoing image/deconvolution to center strong source on pixel");
	  ObitDConCleanVisDeconvolve ((ObitDCon*)myClean, err);
	  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	}
	
	autoCen = 1.0e20;  /* only once */
      }
      
      /* Convergence test */
      if (doSC) {
	if (converged || (SCLoop>=maxASCLoop)) break;
	Obit_log_error(err, OBIT_InfoErr, " ******  Amp&Phase Self Calibration number %d", SCLoop+1);
	
	/* Pass peak found in deconvolution to self cal */
	dim[0] = dim[1] = dim[2] = 1;
	ObitInfoListAlwaysPut(selfCal->info, "peakFlux", OBIT_float, dim, &myClean->peakFlux);
	
	/* Reset minFlux disturbed by Clean */
	dim[0] = dim[1] = dim[2] = 1;
	ftemp = 0.0;
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "noNeg", OBIT_bool, dim, &noNeg);
	
	/* Do self cal */
	converged = ObitUVSelfCalSelfCal (selfCal, inUV, init, &noSCNeed, 
					  myClean->window, err);
	if (err->error) Obit_traceback_msg (err, routine, selfCal->name);
	if (converged || noSCNeed)  break;
	init = FALSE;
	
	/* reset flux limit for next Clean to 1 sigma */
	dim[0] = 1;dim[1] = 1;
	ObitInfoListAlwaysPut (myClean->info, "minFlux", OBIT_float, dim, &selfCal->RMSFld1);
	btemp = FALSE;
	ObitInfoListAlwaysPut(selfCal->skyModel->info, "noNeg", OBIT_bool, dim, &btemp);
	
	/* Possibly reuse some of CLEAN model to start next time */
	if (reuse>0.0) {
	  ftemp = reuse*selfCal->RMSFld1;
	  if (SCLoop==0) ftemp = 1.0e20;  /* Not on first loop */
	  dim[0] = 1;dim[1] = 1;
	  ObitInfoListAlwaysPut (myClean->info, "reuseFlux", OBIT_float, dim, &ftemp);
	} /* end set reuse level */
	
	/* Auto Windowing? */
	if (autoWindow) {
	  myCleanClass->ObitDConCleanAutoWindow ((ObitDConClean*)myClean, 
						 myClean->currentFields, NULL, err);
	  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
	}
	
      } /* end if self cal */
      if (noSCNeed) break;
    
    }  /* End Self cal loop **/
  } /* End Amp&Phase self cal */

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
  if (doFlatten && (myClean->mosaic->numberImages>1)) {
    ObitDConCleanFlatten((ObitDConClean*)myClean, err);
    /* If 2D imaging concatenate CC tables */
    if (!myClean->mosaic->images[0]->myDesc->do3D) 
      ObitImageMosaicCopyCC (myClean->mosaic, inUV, err);
    
    outImage = ObitImageMosaicGetFullImage (myClean->mosaic, err);
  } else  outImage = ObitImageMosaicGetImage (myClean->mosaic, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  
  /* Display? */
  if (selfCal && selfCal->display && outImage)
    ObitDisplayShow (selfCal->display, (Obit*)outImage, NULL, 1, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
   
  /* Cleanup */
  selfCal  = ObitUVSelfCalUnref(selfCal);
  outImage = ObitImageUnref(outImage);

} /* end doImage */

/*----------------------------------------------------------------------- */
/*  Write History for SCMap                                               */
/*   Input:                                                               */
/*      Source    Name of source being imaged                             */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitUV to copy history from                             */
/*      outData   Output Image or UV data to write history to             */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void SCMapHistory (gchar *Source, ObitInfoList* myInput, 
		    ObitUV* inData, ObitData* outData, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "FreqID", "souCode", "Qual", 
    "outFile",  "outDisk", "outName", "outClass", "outSeq",
    "BChan", "EChan", "BIF", "EIF", "BLFact", "BLFOV", "BLchAvg",
    "FOV",  "UVRange",  "timeRange",  "Robust",  "UVTaper",  
    "doCalSelect",  "doCalib",  "gainUse",  "doBand ",  "BPVer",  "flagVer", 
    "doPol",  "PDVer", "doFull", "do3D", "Catalog", "CatDisk", 
    "OutlierDist",  "OutlierFlux",  "OutlierSI",
    "OutlierSize",  "CLEANBox",  "Gain",  "minFlux",  "Niter",  "minPatch",
    "Reuse", "autoCen", "Beam",  "Cmethod",  "CCFilter",  "maxPixel", 
    "autoWindow", "subA", "Alpha", "WtUV", "avgPol", "avgIF", "noNeg", "minSNR",
    "maxPSCLoop", "minFluxPSC", "solPInt", "solPType", "solPMode", 
    "maxASCLoop", "minFluxASC", "solAInt", "solAType", "solAMode", 
    "oldSN", "modelFlux", "modelPos", "modelParm", 
    "nThreads",
    NULL};
  gchar *routine = "SCMapHistory";

  /* error checks */
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitUVIsA(inData));
  g_assert (ObitDataIsA(outData));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory (outData, OBIT_IO_WriteOnly, err);

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

  g_snprintf (hicard, 80, "%s Source = '%s'", pgmName, Source);
  ObitHistoryWriteRec (outHistory, -1, hicard, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Copy selected values from myInput */
  ObitHistoryCopyInfoList (outHistory, pgmName, hiEntries, myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);
  ObitHistoryClose (outHistory, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  inHistory  = ObitHistoryUnref(inHistory);  /* cleanup */
  outHistory = ObitHistoryUnref(outHistory);
 
} /* end SCMapHistory  */

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
  gchar *routine = "BLAvg";

  /* What to do? */
  ObitInfoListGetTest(myInput, "BLFact", &type, dim, &BLFact);
  if (BLFact>1.00) { /* Average */
    /* Set parameters */
    ObitInfoListGetTest(myInput, "BLFOV",   &type, dim, &FOV);
    if (FOV<=0.0) ObitInfoListGetTest(myInput, "FOV",   &type, dim, &FOV);
    ObitInfoListGetTest(myInput, "BLchAvg",   &type, dim, &BLchAvg);
    ObitInfoListGetTest(myInput, "solPInt", &type, dim, &solPInt);
    ObitInfoListGetTest(myInput, "solAInt", &type, dim, &solAInt);
    if (solAInt<=0.0) solAInt = solPInt;
    if (solPInt<=0.0) solPInt = solAInt;
    maxInt = MIN (solPInt, solAInt);
    if (maxInt<=0.0) maxInt = 1.0;  /* Default 1 min */
    /* Default FOV 0.5 lambda/diameter */
    if (FOV<=0.0) {
      Freq = inData->myDesc->crval[inData->myDesc->jlocf];
      FOV = RAD2DG * 0.5 * VELIGHT / (Freq * 25.0);
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

  } else { /* Straight copy */
    ObitUVCopy (inData, outData, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
  }
} /* end BLAvg */
