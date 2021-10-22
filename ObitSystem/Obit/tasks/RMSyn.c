/* $Id$  */
/* Rotation measure synthesis */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2021                                               */
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

#include "ObitImage.h"
#include "ObitRMFit.h"
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
#include "ObitCArray.h"
#include "ObitUtil.h"
#include "ObitPlot.h"
void sincosf(float x, float *sin, float *cos); /* grumble */
#if HAVE_GSL==1  /* GSL stuff */
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_multifit_nlin.h"
#endif /* GSL stuff */

/* internal prototypes */
/* Get inputs */
ObitInfoList* RMSynIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void RMSynOut (ObitInfoList* outList, ObitErr *err);
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
/* Determine Faraday function */
void doRMSyn (ObitInfoList *myInput, ObitImage* inQImage, ObitImage* inUImage,
	      ObitImage* outImage, ObitErr *err);
/* Write history */
void RMSynHistory (ObitInfoList* myInput, ObitImage* inQImage, 
		   ObitImage* outImage, ObitErr* err);
/* Deconvolve cube */
void Decon(ObitInfoList *myInput, ObitImage* outImage, ObitErr *err);
/* Make Beam */
void MakeBeam (ObitInfoList *myInput, ObitImage* inQImage, ObitErr *err);
/* Get pixel RM Function */
void GetRM (olong pixel[2], olong niter, ofloat minFlux);
/* Deconvolve RM Function */
void DeconRM (ObitCArray *fn, ObitCArray *beam, ofloat minRM, olong niter, 
	      olong *nout, olong *pos, ocomplex *flux);

/* Program globals */
gchar *pgmName = "RMSyn";       /* Program name */
gchar *infile  = "RMSyn.inp";   /* File with program inputs */
gchar *outfile = "RMSyn.out";   /* File to contain program outputs */
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
odouble* lamb2=NULL;           /* Lambda^2 per channel  */
olong  nOutChan=0;             /* Number of output channels (nOutChan) */
ObitCArray **RMPlanes=NULL;    /* Planes in RM Image */
ObitCArray *Beam=NULL;         /* RM beam */
ObitCArray *Restor=NULL;       /* restoring function */

/*---------------Private structures----------------*/
/* Threaded function argument -  mostly uses globals */
typedef struct {
  /* ObitThread to use */
  ObitThread *thread;
  /* Size of image to loop over */
  olong        nx,ny;
  /* First row (1-rel) number */
  olong        first;
  /* Highest row (1-rel) number */
  olong        last;
  /* Function dependent arguments */
  olong        niter;
  ofloat       minFlux;
  /* thread number  */
  olong        ithread;
} RMFuncArg;

/** Private: Make Threaded args */
static olong MakeRMFuncArgs (ObitThread *thread, RMFuncArg ***ThreadArgs);

/** Private: Delete Threaded args */
static void KillRMFuncArgs (olong nargs, RMFuncArg **ThreadArgs);

/** Private: Threaded Deconvolve/Restore */
static gpointer ThreadDecon (gpointer arg);


int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Plot images                                                          */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inQImage=NULL, *inUImage=NULL, *outImage= NULL;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean doDecon=FALSE;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = RMSynIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) {ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Get parameters */
  ObitInfoListGetTest(myInput, "doDecon", &type, dim, &doDecon);

  /* Get input images */
  inQImage = getInputImage (myInput, 'Q', err);
  inUImage = getInputImage (myInput, 'U', err);
  outImage = getOutputImage (myInput, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Evaluate Faraday function */
  doRMSyn (myInput, inQImage, inUImage, outImage, err);
  if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}

  /* Deconvolving ? */
  if (doDecon) {
    /* Make complex Beam */
    MakeBeam (myInput, inQImage, err);
    if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}
    
    /* Get RM function, deconvolve, restore */
    Decon (myInput, outImage, err);
    if (err->error) {ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;}
  } /* end deconvolving */

   /* Do history */
  RMSynHistory ( myInput, inQImage, outImage, err);
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

ObitInfoList* RMSynIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "RMSynIn";

  /* Make default inputs InfoList */
  list = defaultInputs(err);

  /* Initialize output */
  myOutput = defaultOutputs(err);
  ObitReturnDumpRetCode (-999, outfile, myOutput, err);
  if (err->error) Obit_traceback_val (err, routine, "GetInput", list);

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
} /* end RMSynIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: RMSyn -input file -output ofile [args]\n");
    fprintf(stderr, "RMSyn Obit task to plot images \n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def RMSyn.in\n");
    fprintf(stderr, "  -output output result file, def RMSyn.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
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

  /* input FITS file name */
  strTemp = "Image.fits";
  dim[0] = strlen (strTemp); dim[1] = 1;
  ObitInfoListPut (out, "inFile", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file name */
  strTemp = "RMSynName";
  dim[0] = 12; dim[1] = 1;
  ObitInfoListPut (out, "inName", OBIT_string, dim, strTemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* input AIPS file class */
  strTemp = "Class ";
  dim[0] = 6; dim[1] = 1;
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
  olong        i;
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

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);


  return inImage;
} /* end getInputImage */

/*----------------------------------------------------------------------- */
/*  Define output image                                                   */
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
/*      outImage  Output cube                                             */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void doRMSyn (ObitInfoList *myInput, ObitImage* inQImage, ObitImage* inUImage, 
	     ObitImage* outImage,  ObitErr *err)
{
  olong i, iplane, nOut, iRM, plane[5]={1,1,1,1,1}, noffset=0;
  olong naxis[2];
  ObitIOSize IOBy;
  ObitInfoType type;
  ObitIOCode retCode;
  union ObitInfoListEquiv InfoReal; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitFArray **inQFArrays=NULL, **inUFArrays=NULL, *workAmp=NULL;
  ObitCArray *workPol=NULL, *workRot=NULL; 
  ObitImageDesc *outDesc=NULL;
  odouble freq, refLamb2;
  olong nlamb2, nx, ny;
  ofloat RM, cmplx[2], *work=NULL, minQ, minU;
  gchar *today=NULL, *RMSYN = "RMSYN   ", keyword[9];
  gboolean doDecon=FALSE;
  gchar *routine = "doRMSyn";

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

  /* Doing deconvolution? no need to write output here */
  ObitInfoListGetTest(myInput, "doDecon", &type, dim, &doDecon);

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

  refLamb2= 1.0e-6;  /* Reference lambda^2 - avoid zero divide */
  
  /* What plane does spectral data start on */
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
 
  /* How many output planes? */
  nOut = (olong)(0.999+(maxRMSyn-minRMSyn)/delRMSyn);
  nOutChan = nOut;  /* to global */

  /* Planes of output image */
  RMPlanes = g_malloc0(nOut*sizeof(ObitCArray*));

  /* Image size */
  nx = inQImage->myDesc->inaxes[0];
  ny = inQImage->myDesc->inaxes[1];
  /* Storage arrays */
  naxis[0] = (olong)nx;  naxis[1] = (olong)ny; 
  workPol  = ObitCArrayCreate (NULL, 2, naxis);
  workRot  = ObitCArrayCreate (NULL, 2, naxis);
  workAmp  = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nlamb2; i++) inQFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nlamb2; i++) inUFArrays[i] = ObitFArrayCreate (NULL, 2, naxis);
  for (i=0; i<nOut; i++)   RMPlanes[i]   = ObitCArrayCreate (NULL, 2, naxis);

  /* Output Image descriptor */
  outImage->myDesc = ObitImageDescCopy (inQImage->myDesc, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inQImage->name);

  /* Change third axis to type "RMSYN  " and set RM range */
  outImage->myDesc->inaxes[outImage->myDesc->jlocf] =  nOut;
  outImage->myDesc->crval[outImage->myDesc->jlocf]  =  VELIGHT/sqrt(refLamb2);
  outImage->myDesc->crpix[outImage->myDesc->jlocf]  =  1.0;
  outImage->myDesc->cdelt[outImage->myDesc->jlocf]  =  delRMSyn;
  outImage->myDesc->crval[outImage->myDesc->jlocf]  =  minRMSyn;
  strncpy (outImage->myDesc->ctype[outImage->myDesc->jlocf], RMSYN,  IMLEN_KEYWORD);
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
  /* Loop reading planes */
  for (iplane=0; iplane<nlamb2; iplane++) {
    /* Lambda^2 Check for MFImage outputs */
    if (!strncmp(inQImage->myDesc->ctype[inQImage->myDesc->jlocf], "SPECLNMF", 8)) {
	sprintf (keyword, "FREQ%4.4d",iplane+1);
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
	ObitInfoListGetTest (inQImage->myDesc->info, keyword, &type, dim, &freq);
       } else {   /* Normal spectral cube */
	freq = inQImage->myDesc->crval[inQImage->myDesc->jlocf] + 
	  inQImage->myDesc->cdelt[inQImage->myDesc->jlocf] * 
	  (inQImage->myDesc->plane - inQImage->myDesc->crpix[inQImage->myDesc->jlocf]);
      }
    lamb2[iplane] = (VELIGHT/freq)*(VELIGHT/freq);

    plane[0] = iplane+noffset+1;  /* Select correct plane */
    retCode = ObitImageGetPlane (inQImage, inQFArrays[iplane]->array, plane, err);
    retCode = ObitImageGetPlane (inUImage, inUFArrays[iplane]->array, plane, err);
    /* if it didn't work bail out */
    if ((retCode!=OBIT_IO_OK) || (err->error)) {
      Obit_traceback_msg (err, routine, inQImage->name);
      goto cleanup;
    }

    /* Deblank */
    ObitFArrayDeblank(inQFArrays[iplane], 0.0);
    ObitFArrayDeblank(inUFArrays[iplane], 0.0);

    /* Plane RMSes */
    QRMS[iplane] = ObitFArrayRMS(inQFArrays[iplane]);
    URMS[iplane] = ObitFArrayRMS(inUFArrays[iplane]);
  } /* end loop reading planes */

  /* Statistics for filtering only bedian within 2 sigma */
  work    = g_malloc0(nlamb2*sizeof(ofloat));  /* work array */
  for (i=0; i<nlamb2; i++) work[i] = QRMS[i];
  qMedian = medianValue(work, 1, nlamb2);
  qSigma  = MedianSigma(nlamb2, work, qMedian);
  for (i=0; i<nlamb2; i++) work[i] = URMS[i];
  uMedian = medianValue(work, 1, nlamb2);
  uSigma  = MedianSigma(nlamb2, work, uMedian);

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

  /* Lower limits */
  minQ = MAX(0.1*qMedian, (qMedian-5*qSigma));
  minU = MAX(0.1*uMedian, (uMedian-5*uSigma));

  /* Loop over RM planes */
  for (iRM=1; iRM<=nOut; iRM++) {
    RM = minRMSyn + (iRM-1)*delRMSyn;
    /* Zero acccumulator */
    cmplx[0] = 0.0; cmplx[1] = 0.0;
    ObitCArrayFill(RMPlanes[iRM-1], cmplx);
    /* loop over input planes */
    for (i=0; i<nlamb2; i++) {
      /* Want this one ? */
      if ((QRMS[i]>qMedian+2*qSigma) || (QRMS[i]<minQ) ||
	  (URMS[i]>uMedian+2*uSigma) || (URMS[i]<minU)) continue;
      sincosf((ofloat)(-2.0*RM*lamb2[i]), &cmplx[1], &cmplx[0]); /* sin/cos factors */
      ObitCArrayFill(workRot, cmplx);
      ObitCArrayComplex(inQFArrays[i], inUFArrays[i], workPol);
      ObitCArrayMul(workPol, workRot, workPol);
      /* Accumulate */
      ObitCArrayAdd(RMPlanes[iRM-1], workPol, RMPlanes[iRM-1]);
    } /* Loop over planes */
    if (!doDecon) {
      /* Get ampl - write to output if needed */
      ObitCArrayAmp(RMPlanes[iRM-1], workAmp);
      plane[0] = iRM;  /* Select correct plane */
      retCode = ObitImagePutPlane (outImage, workAmp->array, plane, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) {
	Obit_traceback_msg (err, routine, outImage->name);
	goto cleanup;
      } /* end write */
    }
   } /* end loop over RM */

 cleanup:
  ObitFArrayUnref(workAmp);
  ObitCArrayUnref(workPol);
  ObitCArrayUnref(workRot);
  for (i=0; i<nlamb2; i++) {
    ObitFArrayUnref(inQFArrays[i]);
    ObitFArrayUnref(inUFArrays[i]);
  }

} /* end doRMSyn */

/*----------------------------------------------------------------------- */
/*  Deconvolve/restore cube                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outImage  Output Image to update                                  */
/*   Output:                                                              */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void Decon (ObitInfoList *myInput, ObitImage* outImage, ObitErr *err)
{
  olong i, iRM;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong niter=20, naxis[2],plane[5]={1,1,1,1,1};
  ofloat minFlux=0.01;
  olong nTh, nRow, loRow, hiRow, nRowPerThread, nThreads, nx, ny;
  ObitFArray *workAmp=NULL;
  RMFuncArg **threadArgs;
  gchar *routine = "Decon";

  /* get parameters */
  ObitInfoListGetTest(myInput, "minFlux", &type, dim, &minFlux);
  ObitInfoListGetTest(myInput, "niter",   &type, dim, &niter);

  /* Initialize Threading */
  nThreads = MakeRMFuncArgs (outImage->thread, &threadArgs);
  /* Divide up work */
  nx = outImage->myDesc->inaxes[0];
  ny = outImage->myDesc->inaxes[1];
  nRow = ny;
  nRowPerThread = nRow/nThreads;
  nTh = nThreads;
  if (nRow<100) {nRowPerThread = nRow; nTh = 1;}
  loRow = 1;
  hiRow = nRowPerThread;
  hiRow = MIN (hiRow, nRow);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hiRow = nRow;  /* Make sure do all */
    threadArgs[i]->nx      = nx;
    threadArgs[i]->ny      = ny;
    threadArgs[i]->first   = loRow;
    threadArgs[i]->last    = hiRow;
    threadArgs[i]->minFlux = minFlux;
    threadArgs[i]->niter   = niter;
    if (nTh>1) threadArgs[i]->ithread = i;
    else threadArgs[i]->ithread = -1;
    /* Update which Row */
    loRow += nRowPerThread;
    hiRow += nRowPerThread;
    hiRow = MIN (hiRow, nRow);
  }
  /* DEBUG 
  fprintf (stderr,"nx %d, ny %d\n",nx,ny);
  for (i=0; i<nTh; i++) {
    fprintf (stderr,"%d, lo %d hi %d\n",i,threadArgs[i]->first,threadArgs[i]->last);
  }*/


  /* Do operation */
   ObitThreadIterator (outImage->thread, nTh, 
		       (ObitThreadFunc)ThreadDecon,
		       (gpointer**)threadArgs);
   
  /* Free local objects */
  KillRMFuncArgs(nThreads, threadArgs);
  
  /* Work FArray */
  naxis[0] = (olong)nx;  naxis[1] = (olong)ny; 
  workAmp  = ObitFArrayCreate (NULL, 2, naxis);
 /* Rewrite */
  for (iRM=1; iRM<=nOutChan; iRM++) {
    ObitCArrayAmp(RMPlanes[iRM-1], workAmp); /* Extract amplitude */
    plane[0] = iRM;  /* Select correct plane */
    ObitImagePutPlane (outImage, workAmp->array, plane, err);
  }
  ObitFArrayUnref(workAmp);
 if (err->error) Obit_traceback_msg (err, routine, outImage->name);
} /* end Decon */

/*----------------------------------------------------------------------- */
/*  Write History for RMSyn                                               */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void RMSynHistory (ObitInfoList* myInput, ObitImage* inImage, 
		    ObitImage* outImage, ObitErr* err)
{
  ObitHistory *inHistory=NULL, *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inQFile",  "inQDisk", "inQName", "inQClass", "inQSeq",
    "inUFile",  "inUDisk", "inUName", "inUClass", "inUSeq", "BeamSig",
    "BLC",  "TRC", "minRMSyn", "maxRMSyn", "delRMSyn", "nThreads",
    NULL};
  gchar *routine = "RMSynHistory";

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
 
} /* end RMSynHistory  */

/*----------------------------------------------------------------------- */
/*  Crude Gaussian approximation to center of beam                        */
/*  Distance to 1st point below 0.333 or first NULL                       */
/*   Input:                                                               */
/*      Beam   Beam to fit                                                */
/*  Return crude Sigma in cells.                                          */
/*----------------------------------------------------------------------- */
ofloat FitGauss(ObitFArray *Beam)
{
  ofloat sigma=0.0;
  ofloat sum, minSum, val, tsigma, newSigma, off;
  olong i, j, k, iCen;
  /* gchar *routine = "FitGauss"; */

  /* error checks */
  if (!Beam) return sigma;

  iCen = Beam->naxis[0]/2; j = 0;
  for (i=iCen+1; i<Beam->naxis[0]; i++, j++) {
    if (Beam->array[i]>Beam->array[i-1]) break;
    sigma = (ofloat)j;
    if (Beam->array[i]<0.33333) break;
  }
  /* tweak fit, brute force paramater search +/- 2 cells */
  minSum=1.0e20;
  for (k=0; k<121; k++) {
    sum = 0.0; 
    tsigma = MAX(0.0, (sigma-2)) + k*0.033;  /* test value */
    for (i=iCen+1; i<iCen+j; i++) {
      off = (iCen-i);
      val = exp (-(off*off)/(tsigma*tsigma));
      sum += (val-Beam->array[i]) * (val-Beam->array[i]);
    }
    if (sum<minSum) {minSum = sum; newSigma = tsigma;}
  }
  return newSigma;
} /* End FitGauss */

/*----------------------------------------------------------------------- */
/*  Generate complex RM "beam", leave as global Beam                      */
/*  Beam is made twice the size of RM syn function to allow deconvolution */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      inImage   Image to copy history from                              */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void MakeBeam(ObitInfoList* myInput, ObitImage* inImage, ObitErr* err)
{
  ofloat *x=NULL, RM, minQ, minU, amp, amp2, maxBeam=-1.0e6;
  olong i, iRM, naxis[2];
  ocomplex accum, SinCos, one, chVal;
  ofloat sigma, off, temp, BeamSig=0.0;
  ObitPlot *plot=NULL;
  ObitInfoType type;
  ObitFArray *workAmp=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar BeamPlot[204], strtmp[204];
  gboolean doPlot=FALSE;
  gchar *routine = "MakeBeam";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(inImage));

  /* Restoring sigma (rad/m^2) specified? */
  ObitInfoListGetTest(myInput, "BeamSig", &type, dim, &BeamSig);

  /* Want Beam plot? */
  ObitInfoListGetTest(myInput, "BeamPlot", &type, dim, BeamPlot);
  BeamPlot[dim[0]] = 0; ObitTrimTrail(BeamPlot);
  doPlot = (BeamPlot!=NULL) && (strlen(BeamPlot)>5) && (BeamPlot[0]!=' ');
  doPlot = FALSE;  /* PLPlot doesn't work */

  /* Allocate global beam array - twice nOutChan */
  naxis[0] = nOutChan*2;
  if (Beam==NULL) Beam = ObitCArrayCreate (NULL, 1, naxis);
  x = g_malloc0(naxis[0]*sizeof(ofloat));

  /* Set real 1 */
  COMPLEX_SET(one, 1., 0.);

  minQ = MAX(0.1*qMedian, (qMedian-5*qSigma));
  minU = MAX(0.1*uMedian, (uMedian-5*uSigma));
  /* Loop over RM planes */
  for (iRM=1; iRM<=nOutChan*2; iRM++) {
    RM = (iRM-nOutChan)*delRMSyn;
    x[iRM-1] = RM;
    /* Zero acccumulator */
    COMPLEX_SET(accum, 0., 0.);

    /* loop over input planes */
    for (i=0; i<nInChan; i++) {
      /* Want this one ? */
      if ((QRMS[i]>qMedian+2*qSigma) || (QRMS[i]<minQ) ||
	  (URMS[i]>uMedian+2*uSigma) || (URMS[i]<minU)) continue;
      sincosf((ofloat)(-2.0*RM*lamb2[i]), &SinCos.imag, &SinCos.real); /* sin/cos factors */
      COMPLEX_MUL2(chVal, one, SinCos);
      COMPLEX_ADD2(accum, accum, chVal);
    } /* end loop over input */
    Beam->array[2*iRM] = accum.real; Beam->array[2*iRM+1] = accum.imag; 
    amp = sqrtf(accum.real*accum.real+accum.imag*accum.imag);
    maxBeam = MAX (amp, maxBeam);
  } /* end loop over RM Planes */
  /* Normalize beam */
  for (iRM=1; iRM<=nOutChan*4; iRM++) 
    {
      Beam->array[iRM-1] /= maxBeam;
      /*fprintf (stderr, "%2.2d %6.2f %6.4f\n",iRM, minRMSyn + (iRM-1)*delRMSyn, Beam->array[iRM-1]);*/
    }

  /* Approximate beam size (pixels) */
  workAmp = ObitCArrayMakeF(Beam);
  ObitCArrayReal(Beam, workAmp);  /* to real part of complex beam */
  if (BeamSig>0.0) sigma = BeamSig/abs(delRMSyn);
  else             sigma = FitGauss(workAmp);
  Obit_log_error(err, OBIT_InfoErr, "Gaussian sigma %5.2f cells, %7.3f Rad/m^2 ", 
		 sigma, sigma*abs(delRMSyn));
  /* Save for history */
  dim[0] = 1;
  temp = sigma*abs(delRMSyn);
  ObitInfoListAlwaysPut (myInput, "BeamSig", OBIT_float, dim, &temp);


  /* Complex restoring function */
  naxis[0] = nOutChan*2;
  if (Restor==NULL) Restor = ObitCArrayCreate (NULL, 1, naxis);
  for (i=0; i<nOutChan*2; i++) {
    off = nOutChan-i;
    amp2 = exp (-(off*off)/(sigma*sigma));
    amp = sqrtf(Beam->array[2*i]*Beam->array[2*i] + Beam->array[2*i+1]*Beam->array[2*i+1]);
    Restor->array[2*i]   = amp2*Beam->array[2*i]/amp;
    Restor->array[2*i+1] = amp2*Beam->array[2*i+1]/amp;
  }

  for (iRM=1; iRM<=nOutChan*2; iRM++) {
    /*fprintf (stderr, "%2.2d %6.2f %6.4f %6.4f\n",iRM, minRMSyn + (iRM-1)*delRMSyn, 
      Beam->array[iRM-1], Restor->array[iRM-1]);*/
  }
  /* PLot? PLPlot doesn't seem to work */
  if (doPlot) {
    plot = newObitPlot ("Plot");
    strncpy(strtmp, BeamPlot,200);
    strncat(strtmp,".ps/psc",200);
    ObitPlotInitPlot (plot, strtmp, 5, 1, 1, err);
    strncpy(strtmp, "Beam/Restore Fn.",200);
    dim[0] = strlen(strtmp);
    ObitInfoListAlwaysPut (plot->info, "TITLE",OBIT_string, dim, strtmp);
    strncpy(strtmp, "Rad/m^2",200);
    dim[0] = strlen(strtmp);
    ObitInfoListAlwaysPut (plot->info, "XLABEL",OBIT_string, dim, strtmp);
    strncpy(strtmp, "Faraday Fn.",200);
    dim[0] = strlen(strtmp);
    ObitInfoListAlwaysPut (plot->info, "YLABEL",OBIT_string, dim, strtmp);
    ObitPlotXYPlot (plot, -2, workAmp->naxis[0], x, workAmp->array, err);
    /*ObitPlotXYOver (plot, -3, workAmp->naxis[0], x, Restor->array, err);*/
    ObitPlotFinishPlot (plot, err);
    ObitPlotUnref(plot);
    if (err->error) Obit_traceback_msg (err, routine, inImage->name);
  } /* end do Plot */
  /* Cleanup */
  if (x) g_free(x);
  ObitFArrayUnref(workAmp);
} /* end MakeBeam  */

/*----------------------------------------------------------------------- */
/*  get RM function for a pixel, deconvolve, restore                      */
/*   Input:                                                               */
/*     pixel   Zero relative pixel coordinate                             */
/*     niter   Max. CLEAN components                                      */
/*     minFlux Min. CLEAN residual                                        */
/*----------------------------------------------------------------------- */
void GetRM(olong pixel[2], olong niter, ofloat minFlux)
{
  ObitCArray *RM=NULL;
  ocomplex flux[1000], *crestor, *cRM, cwork;
  olong i, j, k, cenb, nbeam, iRM, ipx, lniter, nout, pos[1000], naxis[1];

  lniter = MIN(niter,1000); /* Number of iterations */

  /* extract function */
  nbeam = Beam->naxis[0];
  cenb = nbeam/2;  /* center of beam */
  naxis[0] = nOutChan;
  RM = ObitCArrayCreate (NULL, 1, naxis);
  ipx = pixel[0] + pixel[1]*RMPlanes[0]->naxis[0];
  for (i=0; i<nOutChan; i++) {
    RM->array[2*i]   = RMPlanes[i]->array[2*ipx];
    RM->array[2*i+1] = RMPlanes[i]->array[2*ipx+1];
  }

  /*for (iRM=1; iRM<=nOutChan; iRM++) {
      fprintf (stderr, "%2.2d %6.2f %6.4f\n",iRM, minRMSyn + (iRM-1)*delRMSyn, RM[iRM-1]);
      }*/

  /* Deconvolve */
  DeconRM (RM, Beam, minFlux, lniter, &nout, pos, flux);

  /* complex pointers to restor/output array */
  crestor  = (ocomplex*)Restor->array;
  cRM      = (ocomplex*)RM->array;

  /* Complex Restore using restoring beam Restor */
  for (k=0; k<nout; k++) {
    j = cenb - pos[k];
    for (i=0; i<nOutChan; i++,j++) {
      if ((j>0) && (j<nbeam)) {
	COMPLEX_MUL2(cwork, crestor[j], flux[k]); /* gain*flux*restoring_beam */
	COMPLEX_SUB (cRM[i], cRM[i], cwork); /* subtract cwork */
      }
    } /* end inner restore */
  } /* end loop over components */
  
  /* Debug listing */
  if ((pixel[0]==-100) && (pixel[1]==49)) {
    for (iRM=1; iRM<=nOutChan; iRM++) {
      fprintf (stderr, "%2.2d %6.2f %6.4f\n",iRM, minRMSyn + (iRM-1)*delRMSyn, RM->array[iRM-1]);
    }
  }
   /* Put deconvolved/restored pixel back */
  ipx = pixel[0] + pixel[1]*RMPlanes[0]->naxis[0];
  for (i=0; i<nOutChan; i++) {
    RMPlanes[i]->array[2*ipx]   = RM->array[2*i];
    RMPlanes[i]->array[2*ipx+1] = RM->array[2*i+1];
  }
  /* Cleanup */
  RM = ObitFArrayUnref(RM);
} /* end GetRM  */

void DeconRM (ObitCArray *fn, ObitCArray *beam, ofloat minRM, olong niter, 
	      olong *nout, olong *pos, ocomplex *flux)
{
/*----------------------------------------------------------------------- */
/*  Hogbom complex CLEAN deconvolution, gain 0.1                          */
/*   Input:                                                               */
/*      fn    1-D complex RM function to deconvolve                       */
/*      beam  RM complex "Beam", center in nbeam/2 (0 rel)                */
/*      minRM Minimum value to CLEAN                                      */
/*      niter Maximum number of components                                */
/*   Output:                                                              */
/*      fn      RM function residual                                      */
/*      nout    Number of actual components                               */
/*      pos     Pixel numbers (0 rel) of components                       */
/*      flux    complex "Flux" of component                               */
/*----------------------------------------------------------------------- */
  olong n, iter, i, j, cenb,  maxpos, nbeam;
  ofloat maxflux, amp, gain = 0.1;
  ocomplex cwork, val, *cfn, *cbeam;

  /* complex pointers to input/out arrays */
  cfn   = (ocomplex*)fn->array;
  cbeam = (ocomplex*)beam->array;

  *nout = 0;
  n = fn->naxis[0];
  if ((n<0) || (!fn) || (!beam)) return;  /* anything to do? */
  nbeam = Beam->naxis[0];
  cenb = nbeam/2;  /* center of beam */
  /* Outer loop */
  for (iter=0; iter<niter; iter++) {
    /* Next max */
    maxpos = -1; maxflux=-1.0e10; val.real=0.0; val.imag=0.0;
    for (i=0; i<n; i++) {
      amp = sqrtf(cfn[i].real*cfn[i].real + cfn[i].imag*cfn[i].imag);
      if (amp>maxflux) 
	{maxpos=i; maxflux=amp; val.real=cfn[i].real; val.imag=cfn[i].imag;}
    } /* end search for max */
    if (maxflux<minRM) return; /* Done? */
    /* save component */
    pos[iter] = maxpos; 
    COMPLEX_SET(flux[iter], val.real*gain, val.imag*gain);
    *nout = iter+1;
    /* Make residual */
    j = cenb - maxpos;
    for (i=0; i<n; i++,j++) {
      if ((j>0) && (j<nbeam)) {
	COMPLEX_MUL2(cwork, cbeam[j], flux[iter]); /*gain*flux*beam */
	COMPLEX_SUB (cfn[i], cfn[i], cwork); /* cfn[i] - cwork */
      }
    } /* end subtract */
  } /* end outer loop */

} /* end DeconRM */

/**
 * Make arguments for a Threaded ThreadRMFunc
 * \param thread     ObitThread object to be used
 * \param ThreadArgs[out] Created array of RMFuncArg, 
 *                   delete with KillRMFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeRMFuncArgs (ObitThread *thread, RMFuncArg ***ThreadArgs)

{
  olong i, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(RMFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(RMFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->thread= ObitThreadRef(thread);
    (*ThreadArgs)[i]->ithread  = i;
  }

  return nThreads;
} /*  end MakeRMFuncArgs */

/**
 * Delete arguments for ThreadRMFunc
 * \param nargs      number of elements in ThreadArgs.
 * \param ThreadArgs Array of RMFuncArg
 */
static void KillRMFuncArgs (olong nargs, RMFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      if (ThreadArgs[i]->thread) ObitThreadUnref(ThreadArgs[i]->thread);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillRMFuncArgs */

/**
 * Loop over a subset of rows deconvolving/restoring
 * Mostly used program globals
 * Callable as thread
 * \param arg Pointer to RMFuncArg argument with elements:
 * \li first    First Row (1-rel) number
 * \li last     Highest Row (1-rel) number
 * \li ithread  thread number, <0 -> no threading
 * \return NULL
 */
static gpointer ThreadDecon (gpointer arg)
{
  /* Get arguments from structure */
  RMFuncArg *largs = (RMFuncArg*)arg;
  olong      loRow      = largs->first-1;
  olong      hiRow      = largs->last;
  olong      niter      = largs->niter;
  ofloat     minFlux    = largs->minFlux;
  olong      nx         = largs->nx;

  /* local */
  olong      i, j, pixel[2];

  if ((hiRow<loRow) || (niter<=0)) goto finish;

  /* Loop over rows */
  for (i=loRow; i<hiRow; i++) {
     /* if (i==1000) fprintf (stderr,"Row %d\n",i);DEBUG */
    for (j=0; j<nx; j++) {
      pixel[0] = j; pixel[1] = i;
      GetRM (pixel, niter, minFlux);
    } /* end column loop */
  } /* end row loop */

  /* Indicate completion */
  finish: 
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
  
} /*  end ThreadDecon */

