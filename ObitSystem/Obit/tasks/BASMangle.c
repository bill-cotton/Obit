/* $Id:  $  */
/* BASMangle Mangle OTF/BAS images                                   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2016                                               */
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
#include "ObitFArray.h"
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* BASMANGLEIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void BASMANGLEOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Write history */
void BASMANGLEHistory (ObitInfoList* myInput, ObitImage* inImage, 
		       ObitErr* err);

/* Program globals */
gchar *pgmName = "BASMangle";       /* Program name */
gchar *infile  = "BASMangle.inp";   /* File with program inputs */
gchar *outfile = "BASMangle.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL;  /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL;  /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   BASMangle Obit program - compute mean and RMS of an image             */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *inImage= NULL;
  ObitInfoType type;
  ObitErr      *err= NULL;
  ObitImageDesc *Desc=NULL, *IODesc=NULL; 
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong        trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong        i, Aseq, disk, cno, iStokes, iPlane, plane[5]={1,1,1,1,1};
  olong        nStokes=1, nSigma=20, pos[2];
  ofloat       RMS1, RMS, val, maxv=-1.0e5, minv=1.0e5, fblank = ObitMagicF();
  gchar        *strTemp, inFile[128];
  gchar        Aname[13], Aclass[7], *Atype = "MA", *Stok="IQU";

  /* Startup - parse command line */
  err = newObitErr();
  myInput = BASMANGLEIn (argc, argv, err);
  if (err->error) {ierr = 1;  ObitErrLog(err);  goto exit;}

  /* Initialize logging */
  ObitErrInit (err, (gpointer)myInput);

   /* Initialize Threading */
  ObitThreadInit (myInput);

  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Create basic input Image Object */
  inImage = newObitImage("input Image");
  
  /* Get input parameters from myInput */
  ObitInfoListGetTest(myInput, "nStokes", &type, dim, &nStokes); 
  nStokes = MAX (1, MIN(3, nStokes));
  ObitInfoListGetTest(myInput, "nSigma",  &type, dim, &nSigma); 

  /* Loop over Stokes (1-rel) */
  for (iStokes=1; iStokes<=nStokes; iStokes++) {
    
    Obit_log_error(err, OBIT_InfoErr, "Stokes %c", Stok[iStokes-1]);
    /* File type - could be either AIPS or FITS */
    ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&strTemp);
    if (!strncmp (strTemp, "AIPS", 4)) { /* AIPS input */
      /* input AIPS disk */
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
      /* input AIPS name */
      ObitInfoListGet(myInput, "inName", &type, dim, Aname, err);
      /* input AIPS class */
      ObitInfoListGet(myInput, "inClass", &type, dim, Aclass, err);
      if (nStokes>1) Aclass[0] = Stok[iStokes-1];   /* Set Stokes */
      /* input AIPS sequence */
      ObitInfoListGet(myInput, "inSeq", &type, dim, &Aseq, err);
      
      /* if ASeq==0 want highest existing sequence */
      if (Aseq<=0) {
	Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
	/* Save on myInput*/
	dim[0] = dim[1] = 1;
	ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
      } 
      
      /* Find catalog number */
      cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
      
      /* define image */
      ObitImageSetAIPS (inImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
      
    } else if (!strncmp (strTemp, "FITS", 4)) {  /* FITS input */
      /* input FITS file name */
      for (i=0; i<128; i++) inFile[i] = 0;
      ObitInfoListGet(myInput, "inFile", &type, dim, inFile, err);
      
      /* input FITS disk */
      ObitInfoListGet(myInput, "inDisk", &type, dim, &disk, err);
      
      /* define image */
      ObitImageSetFITS (inImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
      
    } else { /* Unknown type - barf and bail */
      Obit_log_error(err, OBIT_Error, "%s: Unknown Image type %s", 
		     pgmName, strTemp);
      if (err->error) ierr = 1;
      ObitErrLog(err); /* show error messages */
      goto exit;
    }
    
    /* error check */
    if (err->error) ierr = 1;
    ObitErrLog(err); /* show any error messages on err */
    if (ierr!=0) goto exit;

    /* Get RMS in plane 1 */
    maxv = -1.0e5; minv = 1.0e5;
    plane[0] = 1;
    ObitImageGetPlane (inImage, NULL, plane, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
    RMS1 = ObitFArrayRMS(inImage->image);
    /* Max, min */
    val = ObitFArrayMax (inImage->image, pos);
    maxv = MAX (maxv, val);
    val = ObitFArrayMin (inImage->image, pos);
    minv = MIN (minv, val);

    /* Blank plane 2 */
    ObitFArrayFill(inImage->image, fblank);
    plane[0] = 2;
    ObitImagePutPlane (inImage, NULL, plane, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

    /* Loop over other planes */
    for (iPlane=3; iPlane<inImage->myDesc->inaxes[2]; iPlane++) {
      plane[0] = iPlane;
      ObitImageGetPlane (inImage, NULL, plane, err);
      if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
      RMS = ObitFArrayRMS(inImage->image);
      if (RMS>nSigma*RMS1) {  /* Blank */
	ObitFArrayFill(inImage->image, fblank);
	ObitImagePutPlane (inImage, NULL, plane, err);
	if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;
      } else { /* Max/min */
	/* Max, min */
	val = ObitFArrayMax (inImage->image, pos);
	maxv = MAX (maxv, val);
	val = ObitFArrayMin (inImage->image, pos);
	minv = MIN (minv, val);
      }
   }

    
    /* Open */
    ObitImageOpen (inImage, OBIT_IO_ReadWrite, err);
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

    /* Update both descriptors */
    Desc=inImage->myDesc; IODesc=(ObitImageDesc*)inImage->myIO->myDesc;
    
    if (maxv!=fblank) {Desc->maxval = maxv;   IODesc->maxval = maxv;}
    else              {Desc->maxval = fblank; IODesc->maxval = fblank;}
    if (minv!=fblank) {Desc->minval = minv;   IODesc->minval = minv;}
    else              {Desc->minval = fblank; IODesc->minval = fblank;}

    /* Mark as modified */
    inImage->myStatus = OBIT_Modified;

    ObitImageClose (inImage, err); /* Close */
    if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

    /* Tell results */
    Obit_log_error(err, OBIT_InfoErr, 
		   "New max %f min %f", maxv, minv);

    /* Do history */
    BASMANGLEHistory (myInput, inImage, err);
    if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  } /* end Stokes loop */
  /* show any messages and errors */
  if (err->error) ierr = 1;
  ObitErrLog(err);
  if (ierr!=0) goto exit;
  
  /* cleanup */
  myInput   = ObitInfoListUnref(myInput);    /* delete input list */
  inImage   = ObitUnref(inImage);
  
  /* Shutdown  */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* BASMANGLEIn (int argc, char **argv, ObitErr *err)
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
  ObitInfoList* list;
  gchar *routine = "BASMANGLEIn";

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
      
    } else if (strcmp(arg, "-nStokes") == 0) { /* Number of Stokes */
      dim[0] = 1;
      itemp = strtol(argv[++ax], NULL, 0);
      ObitInfoListPut (list, "nStokes", OBIT_oint, dim, &itemp, err);
      
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
} /* end BASMANGLEIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: BASMangle -input file -output ofile [args]\n");
    fprintf(stderr, "BASMangle  mangle OTF images from Big Ass Survey\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def BASMANGLE.in\n");
    fprintf(stderr, "  -output output result file, def BASMANGLE.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -nStokes Number of Stokes (IQU) \n");
    
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
  strTemp = "BASMANGLEName";
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

  /* BLC, TRC */
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
/*     Mean     Int        Image pixel mean  [0.0]                        */
/*     RMS      Int        Image pixel rms   [0.0]                        */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat ftemp;
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultOutputs";

  /* add parser items */
  /* Image mean */
  dim[0] = 1; dim[1] = 1;
  ftemp = 0.0;
  ObitInfoListPut (out, "Mean", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Image mean */
  dim[0] = 1; dim[1] = 1;
  ftemp = 0.0;
  ObitInfoListPut (out, "RMS", OBIT_float, dim, &ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Write History for BASMangle                                            */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      outImage  Image to write history to                               */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void BASMANGLEHistory (ObitInfoList* myInput, ObitImage* outImage, 
		      ObitErr* err)
{
  ObitHistory *outHistory=NULL;
  gchar        hicard[81];
  gchar        *hiEntries[] = {
    "DataType", "inFile",  "inDisk", "inName", "inClass", "inSeq",
    "nStokes",  "nSigma", "nThreads", 
    NULL};
  gchar *routine = "BASMANGLEHistory";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitInfoListIsA(myInput));
  g_assert (ObitImageIsA(outImage));

  /* Do history  */
  outHistory = newObitDataHistory ((ObitData*)outImage, OBIT_IO_WriteOnly, err);

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

  outHistory = ObitHistoryUnref(outHistory);
 
} /* end BASMANGLEHistory  */
