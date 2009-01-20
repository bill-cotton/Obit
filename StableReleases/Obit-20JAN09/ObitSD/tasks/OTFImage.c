/* $Id$  */
/* Obit  Radio Single dish On The Fly imaging software                */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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
#include "ObitSystem.h"
#include "ObitMem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"
#include "ObitOTF.h"
#include "ObitOTFCal.h"
#include "ObitOTFUtil.h"
#include "ObitIOOTFFITS.h"
#include "ObitDConCleanOTFRec.h"
#include "ObitHistory.h"
#include "ObitData.h"
#include "ObitDisplay.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitThread.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* OTFImageIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void OTFImageOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Digest inputs */
void digestInputs(ObitInfoList *myInput, ObitErr *err);
/* Get input data */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err);
/* Define output image */
ObitImage* setOutputImage (ObitInfoList *myInput, ObitErr *err);
/* Define output weight image */
ObitImage* setOutputWeight (ObitInfoList *myInput, ObitErr *err);
/* Loop over Channels */
void doChan (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err);
/* Write history */
void OTFImageHistory (ObitInfoList* myInput, ObitOTF* inData, 
		      ObitImage* outImage, ObitErr* err);


/* Program globals */
gchar *pgmName = "OTFImage";       /* Program name */
gchar *infile  = "OTFImage.in" ;   /* File with program inputs */
gchar *outfile = "OTFImage.out";   /* File to contain program outputs */
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
  ObitOTF      *inData = NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = OTFImageIn (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) return 1;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Digest input */
  digestInputs(myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input data */
  inData = getInputData (myInput, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Process */
  doChan (myInput, inData, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* show any messages and errors */
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inData    = ObitUnref(inData);
  
  /*  DEBUG  End DEBUG
      ObitMemPrint (stdout); */
     

  /* Shutdown Obit */
 exit: 
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput  = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* OTFImageIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "OTFImageIn";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return list;

  /* Make default inputs/output InfoList */
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
      
     } else if (strcmp(arg, "-inFile") == 0) { /*inFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "inFile", OBIT_string, dim, strTemp);
      
     } else if (strcmp(arg, "-outFile") == 0) { /*outFile */
      strTemp = argv[++ax];
      dim[0] = strlen (strTemp);
      ObitInfoListAlwaysPut (list, "outFile", OBIT_string, dim, strTemp);

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
} /* end OTFImageIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: OTFImage -input file -output ofile [args]\n");
    fprintf(stderr, "OTFImage Obit task to image OTF data\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def OTFImage.in\n");
    fprintf(stderr, "  -output output result file, def OTFImage.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -inFile input FITS OTF file\n");
    fprintf(stderr, "  -inDisk input image disk number (1-rel) \n");
    fprintf(stderr, "  -outFile output image FITS Image file\n");  
    fprintf(stderr, "  -outDisk output image ((AIPS or FITS) disk number (1-rel) \n");
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
/*     inDisk    Int        input  FITS OTF disk no  [def 1]              */
/*     outDisk   Int        output FITS image disk no  [def 1]            */
/*     outFile   Str [?]    output FITS image file name [def "Image.fits" */
/*     out2Disk  Int        output AIPS or FITS uv disk no  [def 1]       */
/*     Targets   Str (16,1) Targets selected, blank = all                 */
/*     Scans     Int (2)    Scans selected, all                           */
/*     timeRange Flt (2)    Timerange in days , def=all                   */
/*     doCalSelect Boo (1)  Apply calibration/selection?  def=False       */
/*     doCalib   Int (1)    >0 => apply calibration, 2=> cal. wt, def=-1  */
/*     gainUse   Int (1)    Gain table (CL/SN) table to apply, 0=> highest*/
/*     flagVer   Int (1)    Flagging table version, def=0                 */
/*     CLEANBox  Int[4,?]   Clean box, def=all                            */
/*     autoWindow Boo(1)    If true set windows automatically, def=FALSE  */
/*     Gain      Flt (1)    Clean gain, def=0.1                           */
/*     minFlux   Flt (1)    Clean minimum flux density, def=0             */
/*     Niter     Int (1)    Maximum # of CLEAN comp., def=No CLEAN        */
/*     Patch     Int (1)    Clean Min. BEAM half-width, def=100           */
/*     BeamSize  Flt (1)    Clean beam (asec)                             */
/*     deMode    Boo (1)    subtract mode when forming image, def=FALSE   */
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
  strTemp = "OTFData.uvtab";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* FITS disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "inDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS Image file name */
  strTemp = "OTFImageOut.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* output FITS Image disk number */
  dim[0] = 1;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "outDisk", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* TARGETS selected, blank = all */
  strTemp = "                ";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "Targets", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
    
  /* SCANS */
  iarray[0] = 0; iarray[1] = 0;
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Scans", OBIT_string, dim, iarray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Timerange in days */
  dim[0] = 2;dim[1] = 1;
  farray[0] = -1.0e20; farray[1] = 1.0e20;
  ObitInfoListPut (out, "timeRange", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Fraction of peak per major cycle, def=0.75 */
  dim[0] = 1;dim[1] = 1;
  ftemp = 0.75; 
  ObitInfoListPut (out, "fracPeak", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Apply calibration/selection?, def=True */
  dim[0] = 1; dim[1] = 1;
  btemp = TRUE;
  ObitInfoListPut (out, "doCalSelect", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Subtract mode when making image?, def=False */
  dim[0] = 1; dim[1] = 1;
  btemp = FALSE;
  ObitInfoListPut (out, "deMode", OBIT_bool, dim, &btemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  >0 => apply gain calibration, 2=> cal. wt, def=no cal. */
  dim[0] = 1;dim[1] = 1;
  itemp = -1; 
  ObitInfoListPut (out, "doCalib", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /*  Gain table (Cal/Soln) table to apply, 0=> highest, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "gainUse", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Flagging table version, def=0 */
  dim[0] = 1;dim[1] = 1;
  itemp = 0; 
  ObitInfoListPut (out, "flagVer", OBIT_oint, dim, &itemp, err);
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
  ObitInfoListPut (out, "Patch", OBIT_oint, dim, &itemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Clean beam maj, min, PA (asec, asec, deg), def = 0 (fit) */
  dim[0] = 1;dim[1] = 1;
  farray[0] = 0.0; 
  ObitInfoListPut (out, "BeamSize", OBIT_float, dim, farray, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Maximum pixels searched in inner cycle, def = 50000 */
  dim[0] = 1;dim[1] = 1;
  itemp = 50000; 
  ObitInfoListPut (out, "maxPixel", OBIT_oint, dim, &itemp, err);
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
  oint         doCalib;
  gchar        *strTemp, *tname, tmpFile[129];
  ofloat       ftemp;
  gboolean doCalSelect;
  /*gchar *routine = "digestInputs";*/

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));

  /* Make sure doCalSelect set properly */
  doCalSelect = FALSE;
  ObitInfoListGetTest(myInput, "doCalSelect",  &type, dim, &doCalSelect);
  doCalib = -1;
  ObitInfoListGetTest(myInput, "doCalib",  &type, dim, &doCalib);
  doCalSelect = doCalSelect || (doCalib>0);
  ObitInfoListAlwaysPut (myInput, "doCalSelect", OBIT_bool, dim, &doCalSelect);

  /* Set temporary FITS image name as outName */
  /* output FITS file name */
  if (ObitInfoListGetP(myInput, "outFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (tmpFile, strTemp, 128);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    tname = g_strconcat ("tmp", tmpFile, NULL);
    strncpy (tmpFile, tname, 128);
    g_free(tname);
  } else { 
    strncpy (tmpFile, "tmpImage.fits", 128);
  }
  dim[0] = MIN (128, strlen(tmpFile));
  ObitInfoListAlwaysPut (myInput, "outName", OBIT_string, dim, tmpFile);

  /* Convert "RACenter", "DecCenter" to "RA", "Dec" */
  ftemp = 0.0; dim[0] = dim[1] = 1;
  ObitInfoListGetTest(myInput, "RACenter",  &type, dim, &ftemp);
  ObitInfoListAlwaysPut(myInput, "RA",  OBIT_float, dim, &ftemp);
  ftemp = 0.0; dim[0] = dim[1] = 1;
  ObitInfoListGetTest(myInput, "DecCenter",  &type, dim, &ftemp);
  ObitInfoListAlwaysPut(myInput, "Dec",  OBIT_float, dim, &ftemp);

  /* Make sure xCells negative */
  ftemp = 0.0;
  ObitInfoListGetTest(myInput, "xCells",  &type, dim, &ftemp);
  ftemp = -fabs(ftemp);
  ObitInfoListAlwaysPut(myInput, "xCells",  OBIT_float, dim, &ftemp);

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
/*       ObitOTF with input data                                          */
/*----------------------------------------------------------------------- */
ObitOTF* getInputData (ObitInfoList *myInput, ObitErr *err)
{
  ObitOTF       *inData = NULL;
  ObitInfoType type;
  olong         disk, nrec, nThreads;
  gchar        *strTemp, inFile[129];
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "getInputData";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inData;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input OTF data Object */
  inData = newObitOTF("input OTF data");
  
  /* input FITS file name */
  if (ObitInfoListGetP(myInput, "inFile", &type, dim, (gpointer)&strTemp)) {
    strncpy (inFile, strTemp, 128);
  } else { 
    strncpy (inFile, "No_Filename_Given", 128);
  }
  
  /* input FITS disk */
  ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);

  /* define object */
  nrec = 5000;
  nThreads = 1;
  ObitInfoListGetTest(myInput, "nThreads", &type, dim, &nThreads);
  nrec *= nThreads;
  ObitOTFSetFITS (inData, nrec, disk, inFile,  err); 
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);
  
  /* Ensure inData fully instantiated and OK */
  ObitOTFFullInstantiate (inData, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inData);

  return inData;
} /* end getInputData */

/*----------------------------------------------------------------------- */
/*  Loop over frequencies                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to image                                        */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void doChan (ObitInfoList* myInput, ObitOTF* inData, ObitErr* err)
{
  ObitDConCleanOTFRec *myClean=NULL;
  ObitImage    *outField=NULL;
  ObitImage    *outImage=NULL, *outWeight=NULL;
  ObitTableOTFTarget* targetTable=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         ochan, ichan, nchan, chInc, chAvg, BChan, EChan, RChan, 
    bchan, echan, qual;
  gboolean     btemp, autoWindow;
  ofloat       BeamSize, RACenter, DecCenter;
  odouble      RA, Dec;
  olong        ver, inver, outver, plane[5] = {0,1,1,1,1};
  gchar        *targets, target[24], *CCType = "AIPS CC";
  gchar        *dataParms[] = {  /* Parameters to calibrate/select data */
    "Targets", "Scans", "Feeds", "timeRange", "keepCal",
    "doCalSelect", "doCalib", "gainUse", "flagVer", 
    NULL
  };
  gchar        *imgparms[] = {  /* Imaging parameters */
    "RA", "Dec", "xCells", "yCells", 
    "nx", "ny", "minWt", "deMode", "deBias", "Proj", "ConvType", "ConvParm",
    "outName", "outDisk", "doScale", "doFilter",
    NULL
  };
  gchar        *CLEANParms[] = {  /* Clean parameters */
    "CLEANBox", "autoWindow", "Gain", "minFlux", "Niter", "Patch", 
    "fracPeak", "noResid", "doRestore", "dispURL",
    NULL
  };
  gchar *routine = "doChan";

   /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitOTFIsA(inData));

  /* Default position */
  RACenter = 0.0;
  ObitInfoListGetTest(myInput, "RA",  &type, dim, &RACenter);
  DecCenter = 0.0;
  ObitInfoListGetTest(myInput, "Dec",  &type, dim, &DecCenter);
  if ((RACenter==0.0) && (DecCenter==0.0)) {
    if (ObitInfoListGetP(myInput, "Targets", &type, dim, (gpointer*)&targets)) {
      strncpy (target, targets, MIN(24,dim[0])); target[dim[0]]=0;
      if (target[0]!=' ') { /* Check if given */
	ver = 1;
	targetTable = 
	  newObitTableOTFTargetValue ("TargetTable", (ObitData*)inData, &ver, OBIT_IO_ReadOnly, err);
	qual = 0;
	ObitTableOTFTargetGetByName(targetTable, target, qual, &RA, &Dec, err);
	if (err->error) Obit_traceback_msg (err, routine, inData->name);
	targetTable = ObitTableOTFTargetUnref(targetTable);
	RACenter  = (ofloat)RA;
	DecCenter = (ofloat)Dec;
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "RA",  OBIT_float, dim, &RACenter);
	ObitInfoListAlwaysPut(myInput, "Dec",  OBIT_float, dim, &DecCenter);
      }
    }
  }

  /* Define output image */
  outImage = setOutputImage (myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Define output weight image */
  outWeight = setOutputWeight (myInput, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);

  /* Total number of channels */
  nchan = inData->myDesc->inaxes[inData->myDesc->jlocf];

  /* Parameters used here */
  chInc = 1;
  ObitInfoListGetTest(myInput, "chInc",  &type, dim, &chInc);
  chAvg = 0;
  ObitInfoListGetTest(myInput, "chAvg",  &type, dim, &chAvg);
  /* Average everything = Continuum ? */
  if (chAvg<=0) {chAvg = nchan; chInc = nchan;}
  BChan = 1;
  ObitInfoListGetTest(myInput, "BChan",  &type, dim, &BChan);
  BChan = MAX (1, BChan);
  EChan = nchan;
  ObitInfoListGetTest(myInput, "EChan",  &type, dim, &EChan);
  EChan = MIN (EChan, nchan);
  if (EChan<=0) EChan = nchan;
  RChan = 0;
  ObitInfoListGetTest(myInput, "RChan",  &type, dim, &RChan);
  RChan = MIN (RChan, EChan);
  RChan = MAX (BChan, RChan);
  autoWindow = FALSE;
  ObitInfoListGetTest(myInput, "autoWindow", &type, dim, &autoWindow);

  /* Get input parameters from myInput, copy to inData */
  ObitInfoListCopyList (myInput, inData->info, dataParms);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
 
  /* Set imaging parameters */
  ObitInfoListCopyList (myInput, inData->info, imgparms);
    
  /* Make CleanOTFRec */
  myClean = ObitDConCleanOTFRecCreate("Clean Object", inData, err);
  if (err->error) Obit_traceback_msg (err, routine, inData->name);
  
  /* Get input Clean parameters from myInput, copy to myClean */
  ObitInfoListCopyList (myInput, myClean->info, CLEANParms);
  ObitInfoListGetTest(myInput, "BeamSize", &type, dim, &BeamSize);
  BeamSize /= 3600.0;  /* to Deg */
  ObitInfoListAlwaysPut(myClean->info, "BeamSize", type, dim, &BeamSize);
  
  /* Set Clean windows - no it's incomplete
  ObitDConCleanDefWindow((ObitDConClean*)myClean, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name); */
    
  /* Loop over channels */
  ochan = RChan - BChan;
  for (ichan = RChan; ichan<=EChan; ichan+=chInc) {
    ochan++; /* output channel number */
    
    if (BChan<=EChan)
      Obit_log_error(err, OBIT_InfoErr, " **** Start Channel %d", ichan);
    
    /* set selected channels */
    bchan = ichan; 
    echan = bchan + chAvg - 1;
    echan = MIN (echan, nchan);
    dim[0] = 1;
    ObitInfoListAlwaysPut (inData->info, "BChan", OBIT_long, dim, &bchan);
    ObitInfoListAlwaysPut (inData->info, "EChan", OBIT_long, dim, &echan);
    
    /* Automatic windowing  */
    btemp = autoWindow;
    ObitInfoListAlwaysPut (myClean->info, "autoWindow", OBIT_bool, dim, &btemp);
    
    /* Create image, Do CLEAN */
    ObitDConCleanOTFRecDeconvolve ((ObitDCon*)myClean, err);
    if (err->error) Obit_traceback_msg (err, routine, inData->name);
    
    /* Output image stuff */ 
    if (ichan==BChan) {
      /* Create output image(s) */
      outField = myClean->clean;

      /* Image */
      ObitOTFUtilMakeCube (outField->myDesc, inData->myIO->myDesc, 
			   outImage->myDesc, 
			   "I ", BChan, EChan, chInc, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      ObitImageFullInstantiate (outImage, FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);

      /* Weight image */
      ObitOTFUtilMakeCube (outField->myDesc, inData->myIO->myDesc, 
			   outWeight->myDesc, 
			   "I ", BChan, EChan, chInc, err);
      if (err->error) Obit_traceback_msg (err, routine, myClean->name);
      ObitImageFullInstantiate (outWeight, FALSE, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
      /* end of create output */
    } else if ((RChan>BChan) && (ichan==RChan)) { 
      /* Restarting - should already exist */
      ObitImageFullInstantiate (outImage, TRUE, err);
      ObitImageFullInstantiate (outWeight, TRUE, err);
      if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    }
    
    /* Copy result to output */
    plane[0] = ochan;
    outField = myClean->clean;
    /* and copy the CC table */
    inver = 1;
    outver = plane[0];
    ObitDataCopyTable ((ObitData*)outField, (ObitData*)outImage,
		       CCType, &inver, &outver, err);

    /* Image */
    ObitImageUtilInsertPlane (outField, outImage, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    
    /* Weight image */
    ObitImageUtilInsertPlane (myClean->weight, outWeight, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    
  } /* end loop over channels */
  
  /* Do history */
  /* Make sure image created */
  Obit_return_if_fail((outImage!=NULL), err, 
		      "%s: No image generated", routine);
  
  OTFImageHistory (myInput, inData, outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  outImage  = ObitUnref(outImage);
  
  /* Make sure image created */
  Obit_return_if_fail((outWeight!=NULL), err, 
		      "%s: No weight image generated", routine);
  
  OTFImageHistory (myInput, inData, outWeight, err);
  if (err->error) Obit_traceback_msg (err, routine, myClean->name);
  outWeight  = ObitUnref(outWeight);
  
  /* Cleanup */
  if (myClean) {
    /* delete Image Mosaic */
    ObitImageMosaicZapImage (myClean->mosaic, -1, err); /* Delete mosaic members */
    if (err->error) Obit_traceback_msg (err, routine, myClean->name);
    myClean  = ObitDConCleanOTFRecUnref(myClean);
  }

}  /* end doChan */

/*----------------------------------------------------------------------- */
/*  Write History for OTFImage                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inData    ObitOTF to copy history from                            */
/*      outImage  ObitImage to write history to                           */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void OTFImageHistory (ObitInfoList* myInput, ObitOTF* inData, 
		      ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "inFile",  "inDisk", 
    "outDType", "outName", "outClass", "out2Class", "outSeq", "outFile",  "outDisk", 
    "Targets", "Scans", "Feeds", "timeRange",  "keepCal", "minWt", "deMode", "deBias",
    "doScale", "doFilter", "doCalSelect",  "doCalib",  "gainUse", "flagVer", 
    "CLEANBox",  "Gain",  "minFlux",  "Niter",  "Patch",
    "BeamSize",  "fracPeak", "noResid", "doRestore", "fracPeak", 
    "autoWindow", "nThreads",
    NULL};
  gchar *routine = "OTFImageHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitOTFIsA(inData));
  g_assert (ObitImageIsA(outImage));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inData, OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

  /* FITS - copy header */
  ObitHistoryCopyHeader (inHistory, outHistory, err);
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
 
} /* end OTFImageHistory  */

/*----------------------------------------------------------------------- */
/*  Create output image                                                   */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Returns image object                                                 */
/*----------------------------------------------------------------------- */
ObitImage* setOutputImage (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  gchar     outFile[129], tmpFile[129], *tname, *outName, *strTemp, *Type=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  olong      Aseq, cno;
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong      i, disk;
  gboolean  exist;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputImage";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  
  /* Create basic output Image Object */
  outImage = newObitImage("Output Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */ 
    /* Output AIPS name */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    if ((outName!=NULL) && (strlen(outName)>0) && (outName[0]!=' ') && (outName[1]!=' ')) {
      strncpy (Aname, outName, 13); 
    } else { /* No output name given */
      strncpy (Aname, "noname", 13); 
    }
    Aname[12] = 0;

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "outClass", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Map  ", 7);
    }
    if ((strlen(Aclass)<=0) || (Aclass[0]==' ') || (Aclass[1]==' '))
      strncpy (Aclass, "Map  ", 7);
    Aclass[6] = 0;

    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

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

    /* Get name */
    ObitInfoListGet (myInput, "outFile", &type, dim, tmpFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    
    /* add .fits */
    tname = g_strconcat (tmpFile, ".fits", NULL);
    strncpy (outFile, tname, 128);
    g_free(tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    
    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } /* End FITS output image */
 
 return outImage;
} /* end setOutputImage  */

/*----------------------------------------------------------------------- */
/*  Create output weight image                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*   Returns image object                                                 */
/*----------------------------------------------------------------------- */
ObitImage* setOutputWeight (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage *outImage=NULL;
  ObitInfoType type;
  gchar     outFile[129], tmpFile[129], *tname, *outName, *strTemp, *Type=NULL;
  gint32    dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar     Aname[13], Aclass[7], *Atype = "MA";
  olong      Aseq, cno;
  olong      blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong      trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong      i, disk;
  gboolean  exist;
  gchar     *FITS = "FITS";
  gchar     *routine = "setOutputWeight";

  /* error checks */
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));
  
  /* Create basic output Image Object */
  outImage = newObitImage("Output Image");
  
  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "outDType", &type, dim, (gpointer)&Type);
  if ((Type==NULL) || (!strncmp(Type,"    ",4))) Type = FITS;
  
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS output */ 
    /* Output AIPS name */
    ObitInfoListGetP (myInput, "outName", &type, dim, (gpointer)&outName);
    for (i=0; i<12; i++) Aname[i] = ' '; Aname[i] = 0;
    if ((outName!=NULL) && (strlen(outName)>0) && (outName[0]!=' ') && (outName[1]!=' ')) {
      strncpy (Aname, outName, 13); 
    } else { /* No output name given */
      strncpy (Aname, "noname", 13); 
    }
    Aname[12] = 0;

    /* output AIPS class */
    if (ObitInfoListGetP(myInput, "out2Class", &type, dim, (gpointer)&strTemp)) {
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "Wt  ", 7);
    }
    if ((strlen(Aclass)<=0) || (Aclass[0]==' ') || (Aclass[1]==' '))
      strncpy (Aclass, "Wt  ", 7);
    Aclass[6] = 0;

    /* AIPS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "outSeq", &type, dim, &Aseq, err);

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

    /* Get name */
    ObitInfoListGet (myInput, "outFile", &type, dim, tmpFile, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    ObitTrimTrail(tmpFile);  /* remove trailing blanks */
    
    /* add .fits */
    tname = g_strconcat (tmpFile, ".wt", NULL);
    strncpy (outFile, tname, 128);
    g_free(tname);
    ObitTrimTrail(outFile);  /* remove trailing blanks */
    
    /* output FITS disk */
    ObitInfoListGet(myInput, "outDisk", &type, dim, &disk, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
    
    /* define object */
    ObitImageSetFITS (outImage, OBIT_IO_byPlane, disk, outFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  } /* End FITS output weight image */
 
 return outImage;
} /* end setOutputWeight  */
