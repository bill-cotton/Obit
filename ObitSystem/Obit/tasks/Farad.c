/* $Id$  */
/* Rotation measure analysis  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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

#include "ObitCArray.h"
#include "ObitFaraSyn.h"
#include "ObitImage.h"
#include "ObitThread.h"
#include "ObitSinCos.h"
#include "ObitComplex.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitFArray.h"
#include "ObitUtil.h"
#include "ObitPlot.h"
void sincosf(float x, float *sin, float *cos); /* grumble */

/* internal prototypes */
/* Get inputs */
ObitInfoList* FaradIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void FaradOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get input image */
ObitImage* getInputImage (ObitInfoList *myInput, gchar Stok, ObitErr *err);
/* Define output image */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitErr *err);
/* Do Faraday analysis */
void doFarad (ObitInfoList *myInput, ObitImage* inQImage, ObitImage* inUImage,
	      ObitImage* outImage, ObitErr *err);
/* Write history */
void FaradHistory (ObitInfoList* myInput, ObitImage* inQImage, 
		   ObitImage* outImage, ObitErr* err);
/* Program globals */
ObitSystem   *mySystem= NULL;
gchar *pgmName = "Farad";       /* Program name */
gchar *infile  = "Farad.inp";   /* File with program inputs */
gchar *outfile = "Farad.out";   /* File to contain program outputs */
olong  pgmNumber;      /* Program number (like POPS no.) */
olong  AIPSuser;       /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;        /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;        /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
/* Working */
ofloat maxRMSyn, minRMSyn, delRMSyn; /* RM selection */
olong  nInChan=0;              /* Number of input channels */
ofloat *QRMS=NULL, *URMS=NULL; /* Q/U Channel rmses (nInChan)  */
ofloat qMedian, uMedian, qSigma, uSigma;  /* Statistics */
ofloat *specCor=NULL;          /* chennel spectral correction */
odouble* lamb2=NULL;           /* Lambda^2 per channel  */
double   refLamb2;             /* reference Lambda^2 */
olong  nOutChan=0;             /* Number of output channels (nOutChan) */
/*---------------Private structures----------------*/

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Plot images                                                          */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitImage    *inQImage=NULL, *inUImage=NULL, *outImage= NULL;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = FaradIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get input images */
  inQImage = getInputImage (myInput, 'Q', err);
  inUImage = getInputImage (myInput, 'U', err);
  outImage = getOutputImage (myInput, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Evaluate Faraday function */
  doFarad (myInput, inQImage, inUImage, outImage, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Do history */
  FaradHistory ( myInput, inQImage, outImage, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}

  /* cleanup */
  myInput    = ObitInfoListUnref(myInput);    /* delete input list */
  inQImage   = ObitUnref(inQImage);
  inUImage   = ObitUnref(inUImage);
  outImage   = ObitUnref(outImage);
 
  /* Shutdown  */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* FaradIn (int argc, char **argv, ObitErr *err)
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
  oint    itemp, i, j, k, iarr[IM_MAXDIM];
  ObitInfoList* list;
  gchar *routine = "FaradIn";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

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
      
     } else if (strcmp(arg, "-BLC") == 0) { /* BLC */
      dim[0] = 1;
      /* read until something starts with "-" of hit end */
      i = 0;
      while ( ((ax+1)<argc) && (argv[ax+1][0]!='-')) {
	iarr[i++] = strtol(argv[++ax], NULL, 0);
      }
      dim[0] = i;
      ObitInfoListAlwaysPut (list, "BLC", OBIT_oint, dim, &iarr);
      
     } else if (strcmp(arg, "-TRC") == 0) { /* TRC */
      dim[0] = 1;
      /* read until something starts with "-" of hit end */
      i = 0;
      while ( ((ax+1)<argc) && (argv[ax+1][0]!='-')) {
	iarr[i++] = strtol(argv[++ax], NULL, 0);
      }
      dim[0] = i;
      ObitInfoListAlwaysPut (list, "TRC", OBIT_oint, dim, &iarr);
      
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

  /* Initialize Threading */
  ObitThreadInit (list);

  return list;
} /* end FaradIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Farad -input file -output ofile [args]\n");
    fprintf(stderr, "Farad Obit task to do Faraday analysis \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Farad.in\n");
    fprintf(stderr, "  -output output result file, def Farad.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inQFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inQName input AIPS file name\n");
    fprintf(stderr, "  -inQClass input AIPS file class\n");
    fprintf(stderr, "  -inQSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inQDisk input (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BLC bottom-left (1-rel) pixel of image def. [1,1,..] \n");
    fprintf(stderr, "  -TRC top-right (1-rel) pixel of image def. [nx,ny,..] \n");
    /*/exit(1);  bail out */
  }/* end Usage */

/*----------------------------------------------------------------------- */
/*  Create default input ObitInfoList                                     */
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
/*     inFile    Str [?]    input FITS image file name [def "Image.fits"] */
/*     inName    Str [12]   input AIPS image name  [no def]               */
/*     inClass   Str [6]    input AIPS image class  [no def]              */
/*     inSeq     Int        input AIPS image sequence no  [no def]        */
/*     inDisk    Int        input AIPS or FITS image disk no  [def 1]     */
/*     blc       Int  [7]   bottom-left (1-rel) corner[def {1,1,...)]     */
/*     trc       Int  [7]   top-right (1-rel) corner[def {0,0,...)]       */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultInputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *strTemp;
  oint   itemp;
  olong   blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong   trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultInputs";

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

  /* output FITS file names */
  strTemp = "FaradImage.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "outFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "Name";
  dim[0] = 12; dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file class */
  strTemp = "QPol";
  dim[0] = 6; dim[1] = 1;
  ObitInfoListPut (out, "inClass", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  strTemp = "UPol";
  ObitInfoListPut (out, "in2Class", OBIT_string, dim, strTemp, err);
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

  /* BLC*, TRC* */
  dim[0] = IM_MAXDIM;dim[1] = 1;
  itemp = 1; 
  ObitInfoListPut (out, "BLC", OBIT_oint, dim, blc, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);
  ObitInfoListPut (out, "TRC", OBIT_oint, dim, trc, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultInputs */

/*----------------------------------------------------------------------- */
/*  Create default output ObitInfoList                                    */
/*   Return                                                               */
/*       ObitInfoList  with default values                                */
/*  Values:                                                               */
/*     Center    Center Pixel  [0., 0.]                                   */
/*     Width     Width in pixels  [0., 0.]                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  /* gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
     ofloat ftemp[2] = {0.0, 0.0};*/
  ObitInfoList *out = newObitInfoList();
  /*gchar *routine = "defaultOutputs";*/

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Get input image, sets BLC, TRC selection                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      Stokes    'Q' or 'U'                                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with first input image                                 */
/*----------------------------------------------------------------------- */
ObitImage* getInputImage (ObitInfoList *myInput, gchar Stok, ObitErr *err)
{
  ObitImage    *inImage = NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        i, noffset;
  gchar *routine = "getInputImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inImage;
  g_assert (ObitInfoListIsA(myInput));

  /* By type */
  if (Stok=='Q')      inImage = ObitImageFromFileInfo ("inQ", myInput, err);
  else if (Stok=='U') inImage = ObitImageFromFileInfo ("inU", myInput, err);
  else  Obit_log_error(err, OBIT_Error, "%s: Unknown Stokes type %c", 
                   routine, Stok);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);

  /* Ensure inImage fully instantiated and OK */
  ObitImageFullInstantiate (inImage, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

   /* Set defaults BLC, TRC */
  for (i=0; i<IM_MAXDIM; i++) {
    if (blc[i]<=0) blc[i] = 1;
    blc[i] = MAX (1,  blc[i]);
    if (trc[i]<=0) trc[i] = inImage->myDesc->inaxes[i];
    trc[i] = MIN (trc[i], inImage->myDesc->inaxes[i]);
  }

  /* What 0-rel plane does spectral data start on? */
  if (!strncmp(inImage->myDesc->ctype[inImage->myDesc->jlocf], "SPECLNMF", 8)) {
    noffset = 2;  /* What plane does spectral data start on */
    ObitInfoListGetTest (inImage->myDesc->info, "NTERM", &type, dim, &noffset);
  } else {   /* Normal spectral cube */
    noffset = 0;
  }
  /* start with first spectral plane or later */
  blc[2] = MAX(blc[2],noffset+1);

  /* Make sure at least two planes */
  Obit_retval_if_fail(((trc[2]-blc[2])>=2), err, inImage,
		      "%s: MUST have at least two planes", routine);

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  /* Why??? for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;*/
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  return inImage;
} /* end getInputImage */

/*----------------------------------------------------------------------- */
/*  Define output image cube                                              */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage for output image                                       */
/*----------------------------------------------------------------------- */
ObitImage* getOutputImage (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *outImage = NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean F=FALSE;
  gchar *routine = "getOutputImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return outImage;
  g_assert (ObitInfoListIsA(myInput));

  /* May not exist */
  ObitInfoListAlwaysPut(myInput, "outExist", OBIT_bool, dim, &F);

  outImage = ObitImageFromFileInfo ("out", myInput, err);

  if (err->error) Obit_traceback_val (err, routine, "myInput", outImage);
  return outImage;
} /* end getOutputImage */

/*----------------------------------------------------------------------- */
/*  Evaluate RM function                                                  */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inQImage  Input Q Image                                           */
/*      inUImage  Input U Image                                           */
/*   Output:                                                              */
/*      outImage   Output cube                                            */
/*      err        Obit Error stack                                       */
/*----------------------------------------------------------------------- */
void doFarad (ObitInfoList *myInput, ObitImage* inQImage, ObitImage* inUImage, 
	      ObitImage* outImage, ObitErr *err)
{
  olong i, iplane, nOut, plane[5]={1,1,1,1,1}, noffset=0;
  unsigned short jplane;
  olong lambplane, naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong   blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong   trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ObitFArray **inQFArrays=NULL, **inUFArrays=NULL, *workAmp=NULL;
  ObitImageDesc *outDesc=NULL;
  ObitFaraSyn *myFaraSyn=NULL;
  gboolean doRMSyn=TRUE, doError=TRUE;
  ofloat minQUSNR=0.1, minFrac=0.5, maxChi2=1000.0;
  odouble freq, refFreq;
  olong nlamb2, nx, ny;
  ofloat alpha=0.0, *plnWt=NULL, fact, *work=NULL, fblank = ObitMagicF();
  gchar *today=NULL, keyword[12];
  gchar *SPECRM   = "SPECRM  ", *MaxRMSyn   = "MaxRMSyn";
  gchar *routine = "doFarad";
  /* DEBUG 
  olong pos[2] = {100,100};
  ofloat *cptr;*/
  if (err->error) return; /* previous error? */
  g_assert(ObitImageIsA(inQImage));
  g_assert(ObitImageIsA(inUImage));
  g_assert(ObitImageIsA(outImage));

  /*  Max RM for RM syn function */
  InfoReal.flt = 0.0; type = OBIT_float;
  maxRMSyn  = 0.0;
  ObitInfoListGetTest(myInput, "maxRMSyn", &type, dim, &InfoReal);
  if (type==OBIT_float)       maxRMSyn = InfoReal.flt;
  else if (type==OBIT_double) maxRMSyn = (ofloat)InfoReal.dbl;

  /*  Min RM for RM syn function */
  InfoReal.flt = 1.0; type = OBIT_float;
  minRMSyn = 1.0;
  ObitInfoListGetTest(myInput, "minRMSyn", &type, dim, &InfoReal);
  if (type==OBIT_float)       minRMSyn = InfoReal.flt;
  else if (type==OBIT_double) minRMSyn = (ofloat)InfoReal.dbl;

  /* Delta RM for RM syn function */
  InfoReal.flt = 1.0; type = OBIT_float;
  delRMSyn = 1.0;
  ObitInfoListGetTest(myInput, "delRMSyn", &type, dim, &InfoReal);
  if (type==OBIT_float)       delRMSyn = InfoReal.flt;
  else if (type==OBIT_double) delRMSyn = (ofloat)InfoReal.dbl;

  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* Doing only direct search */
  ObitInfoListGetTest(myInput, "doRMSyn", &type, dim, &doRMSyn);

  /* Error analysis? */
  ObitInfoListGetTest(myInput, "doError",   &type, dim, &doError);

  /* other parameters */
  ObitInfoListGetTest(myInput, "minQUSNR",  &type, dim, &minQUSNR);
  ObitInfoListGetTest(myInput, "minFrac",   &type, dim, &minFrac);
  ObitInfoListGetTest(myInput, "maxChi2",   &type, dim, &maxChi2);

  /* Open input images to get info */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (inQImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inQImage->extBuffer = FALSE;
  retCode = ObitImageOpen (inQImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inQImage->name);
  ObitInfoListAlwaysPut (inUImage->info, "IOBy", OBIT_long, dim, &IOBy);
  inUImage->extBuffer = FALSE;
  retCode = ObitImageOpen (inUImage, OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUImage->name);

  /* Check compatability */
  Obit_return_if_fail(((inQImage->myDesc->inaxes[0]==inUImage->myDesc->inaxes[0]) && 
		       (inQImage->myDesc->inaxes[1]==inUImage->myDesc->inaxes[1]) &&
		       (inQImage->myDesc->inaxes[2]==inUImage->myDesc->inaxes[2])), err,
		      "%s: Input images incompatable", routine);

  refFreq = inQImage->myDesc->crval[inQImage->myDesc->jlocf]; /* reference frequency */
  refLamb2= (VELIGHT/refFreq)*(VELIGHT/refFreq);  /* Reference lambda^2  */
  refLamb2= 1.0e-6;  /* Reference lambda^2 - no zero divide */
  
  /* Save refLamb2 for History */
  Obit_log_error(err, OBIT_InfoErr, "refLambda2=%f ",refLamb2);
  dim[0] = dim[1] = dim[2] = 1;;
  ObitInfoListAlwaysPut (myInput, "refLambda2", OBIT_double, dim, &refLamb2);

  /* What plane does spectral data start on? */
  if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
    noffset = 2;  /* What plane does spectral data start on */
    ObitInfoListGetTest (inQImage->myDesc->info, "NTERM", &type, dim, &noffset);
  } else {   /* Normal spectral cube */
    noffset = 0;
  }

  /* Determine number of frequency planes and initialize in */
  nlamb2  = inQImage->myDesc->inaxes[inQImage->myDesc->jlocf]-noffset;
  nInChan = nlamb2;  /* to global */
  lamb2   = g_malloc0(nlamb2*sizeof(odouble));
  inQFArrays = g_malloc0(nlamb2*sizeof(ObitFArray*));
  inUFArrays = g_malloc0(nlamb2*sizeof(ObitFArray*));
  QRMS       = g_malloc0(nlamb2*sizeof(ofloat));
  URMS       = g_malloc0(nlamb2*sizeof(ofloat));
  plnWt      = g_malloc0(nlamb2*sizeof(ofloat));
 
  /* How many output planes? nterm always 2 */
  if (doError) nOut = 2+2*2;
  else         nOut = 2+2;
  /* always 3 for doRMSyn */
  if (doRMSyn) nOut = 3;
  nOutChan = nOut;  /* to global */

  /* Image size */
  nx = inQImage->myDesc->inaxes[0];
  ny = inQImage->myDesc->inaxes[1];
  /* Storage arrays */
  naxis[0] = (olong)nx;  naxis[1] = (olong)ny; 
  workAmp  = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nlamb2; i++) inQFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nlamb2; i++) inUFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inQImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

    /* Change third axis to type "SPECRM  " and leave the reference frequency
     as the "CRVAL" */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
  outImage->myDesc->crval[outImage->myDesc->jlocf]  =  VELIGHT/sqrt(refLamb2);
  outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
  outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  1.0;
  if (doRMSyn) strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], MaxRMSyn, IMLEN_KEYWORD);
  else         strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], SPECRM ,  IMLEN_KEYWORD);
  outImage->myDesc->bitpix = -32;  /* Float it */

  /* Creation date today */
  today = ObitToday();
  strncpy (outImage->myDesc->date, today, IMLEN_VALUE);
  if (today) g_free(today);

  /* save descriptor */
  outDesc = ObitImageDescCopy (outImage->myDesc, outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

  /* Open output */
  IOBy = OBIT_IO_byPlane;   dim[0] = 1;
  ObitInfoListAlwaysPut (outImage->info, "IOBy", OBIT_long, dim, &IOBy);
  outImage->extBuffer = TRUE;
  retCode = ObitImageOpen (outImage, OBIT_IO_WriteOnly, err);

  /* Update descriptor */
  outImage->myDesc = ObitImageDescCopy (outDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

  /* Close to update*/
  retCode = ObitImageClose (outImage, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, outImage->name);
    goto cleanup;
  }

  Obit_log_error(err, OBIT_InfoErr, "Loop reading planes, nx=%d, ny=%d",nx,ny);
  ObitSystemUsage (mySystem, err); /* Timing messages */

  specCor = g_malloc0(nlamb2*sizeof(ofloat));  /* Channel Weights */
  /* Loop reading planes */
  for (iplane=1; iplane<=nlamb2; iplane++) {
    lambplane= iplane-1; /* Plane in lambda^2 cubes */
    /* Lambda^2 Check for MFImage outputs */
    if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
      jplane = iplane+blc[2]-noffset-1; /* reduce compiler bitching */
      sprintf (keyword, "FREQ%4.4d",jplane);
      /* In case keyword missing */
      freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	(inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
      ObitInfoListGetTest (inQImage->myDesc->info, keyword, &type, dim, &freq);
    } else {   /* Normal spectral cube */
      freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	(inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
    }
    lamb2[lambplane] = (VELIGHT/freq)*(VELIGHT/freq);

    plane[0] = blc[2]+lambplane;  /* Select correct plane */
    retCode = ObitImageGetPlane (inQImage, inQFArrays[lambplane]->array, plane, err);
    retCode = ObitImageGetPlane (inUImage, inUFArrays[lambplane]->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) {
      Obit_traceback_msg (err, routine, inQImage->name);
      goto cleanup;
    }

    /* Deblank */
    ObitFArrayDeblank(inQFArrays[lambplane], 0.0);
    ObitFArrayDeblank(inUFArrays[lambplane], 0.0);

    /* Plane RMSes */
    QRMS[lambplane] = ObitFArrayRMS(inQFArrays[lambplane]);
    URMS[lambplane] = ObitFArrayRMS(inUFArrays[lambplane]);

    /* Plane weight - the average of 1/QRMS and 1/URMS */
    if ((QRMS[lambplane]!=0.0) && (QRMS[lambplane]!=fblank) &&
	(URMS[lambplane]!=0.0) && (URMS[lambplane]!=fblank)) {
      /* plnWt[lambplane] = 0.5*((1./QRMS[lambplane])+(1./URMS[lambplane]));*/
      plnWt[lambplane] = 1.0;  /* no weighting */
    } else plnWt[lambplane] = 0.0;
    fact = plnWt[lambplane]; /* factor for scaling & weighting */

    /* spectral correction (alpha!=0) */
    specCor[lambplane] = (ofloat)exp(-alpha*log(freq/refFreq));
    if (alpha!=0.0) fact *= specCor[lambplane];

    /* Scale for spectral index and weighting */
    ObitFArraySMul(inQFArrays[lambplane], fact);
    ObitFArraySMul(inUFArrays[lambplane], fact);
    plnWt[lambplane] = fact;  /* Full plane weighting */
  } /* end loop reading planes */

  /* Statistics to globals */
  work    = g_malloc0(nlamb2*sizeof(ofloat));  /* work array */
  for (i=0; i<nlamb2; i++) work[i] = QRMS[i];
  qMedian = medianValue(work, 1, nlamb2);
  qSigma  = MedianSigma(nlamb2, work, qMedian);
  for (i=0; i<nlamb2; i++) work[i] = URMS[i];
  uMedian = medianValue(work, 1, nlamb2);
  uSigma  = MedianSigma(nlamb2, work, uMedian);
  if (work) g_free(work);

  /* Close inputs */
  retCode = ObitImageClose (inQImage, err);
  inQImage->extBuffer = FALSE;   /* May need I/O buffer later */
  retCode = ObitImageClose (inUImage, err);
  inUImage->extBuffer = FALSE;   /* May need I/O buffer later */
  /* if it didn't work bail out */
  if ((retCode!=OBIT_IO_OK) || (err->error)) {
    Obit_traceback_msg (err, routine, inQImage->name);
    goto cleanup;
  }

  /* Drop bad lamb2 planes */
  for (i=0; i<nlamb2; i++) {
    if ((plnWt[i]<=0.0) || (QRMS[i]>qMedian+3*qSigma) || (URMS[i]>uMedian+3*uSigma)) {
      inQFArrays[i] = ObitFArrayUnref(inQFArrays[i]);
      inUFArrays[i] = ObitFArrayUnref(inUFArrays[i]);
      plnWt[i] = 0.0;
    }
  } /* end drop bad planes */
  
  //Obit_log_error(err, OBIT_InfoErr, "Begin");
  //ObitSystemUsage (mySystem, err); /* Timing messages */
  myFaraSyn = ObitFaraSynCreate 
    ("myFaraSyn",
     nlamb2, lamb2, refLamb2, inQFArrays, inUFArrays, plnWt, 
     QRMS, URMS, nOut, minRMSyn, maxRMSyn, delRMSyn, NULL, FALSE,
     NULL, NULL, outImage, minQUSNR, minFrac, doError, doRMSyn, maxChi2);

  /* Do analysis */
  ObitFaraSynRMAna(myFaraSyn, err);

  //Obit_log_error(err, OBIT_InfoErr, "End");
  //ObitSystemUsage (mySystem, err); /* Timing messages */

 cleanup:
  if (plnWt) g_free(plnWt);
  ObitFArrayUnref(workAmp);
  myFaraSyn = ObitFaraSynUnref(myFaraSyn);
  for (i=0; i<nlamb2; i++) {
    ObitFArrayUnref(inQFArrays[i]);
    ObitFArrayUnref(inUFArrays[i]);
  }

} /* end doFarad */

/*----------------------------------------------------------------------- */
/*  Write History for Farad                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void FaradHistory (ObitInfoList* myInput, ObitImage* inImage, 
		   ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inQFile",  "inQDisk", "inQName", "inQClass", "inQSeq",
    "inUFile",  "iUnDisk", "inUName", "inUClass", "inUSeq",
    "BLC", "TRC", "minRMSyn", "maxRMSyn", "delRMSyn", "nThreads",
    "doRMSyn", "doError", "minFrac", "minQUSNR", "refLambda2", 
    NULL};
  gchar *routine = "FaradHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));
  g_assert (ObitImageIsA(outImage));

  /* Do history  */
  inHistory  = newObitDataHistory ((ObitData*)inImage,  OBIT_IO_ReadOnly, err);
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

  /* If FITS copy header */
  if (inHistory->FileType==OBIT_IO_FITS) {
    ObitHistoryCopyHeader (inHistory, outHistory, err);
    ObitHistoryCopy (inHistory, outHistory, err);
 } else { /* simply copy history */
     ObitHistoryCopy (inHistory, outHistory, err);
  }
  if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  
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
 
} /* end FaradHistory  */

