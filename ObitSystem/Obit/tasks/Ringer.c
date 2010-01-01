/* $Id$  */
/* Fit Rings to SiO maser images                                      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2007,2009                                          */
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
#include "ObitSystem.h"
#include "ObitParser.h"
#include "ObitReturn.h"
#include "ObitAIPSDir.h"

/* internal prototypes */
/* Get inputs */
ObitInfoList* RingerIn (int argc, char **argv, ObitErr *err);
/* Set outputs */
void RingerOut (ObitInfoList* outList, ObitErr *err);
/* Give basic usage on error */
void Usage(void);
/* Set default inputs */
ObitInfoList* defaultInputs(ObitErr *err);
/* Set default outputs */
ObitInfoList* defaultOutputs(ObitErr *err);
/* Get first input image */
ObitImage* getInputImage1 (ObitInfoList *myInput, ObitErr *err);
/* Get second input image */
ObitImage* getInputImage2 (ObitInfoList *myInput, ObitErr *err);
/* Check that inputs are compatable */
void checkInput (ObitImage* in1Image, ObitImage* in2Image, ObitErr *err);
/* do fitting of circular ring */
void doFit (ObitInfoList *myInput, ObitImage* in1Image, ObitImage* in2Image, 
	    ObitInfoList *myOutput, ObitErr *err);
/* parameter search */
void search (ofloat space, olong nsearch, ofloat spacing, olong prtLv, 
	     ofloat *Center, ofloat *Radius, ofloat *Width, ofloat Frac[2], 
	     ObitErr *err);
/* Histograms about new center */
void DoHisto(ofloat cen[2], ofloat spacing);
/* get moments about center */
void DoMom(gboolean final, ofloat Amom[2], ofloat Bmom[2]);
/* get Fraction of flux in annulus */
void GetFrac (ofloat Amom[2], ofloat Bmom[2], ofloat Frac[2]);
/* do fitting of Elliptical ring */
void doFitE (ObitInfoList *myInput, ObitImage* in1Image, ObitImage* in2Image, 
	    ObitInfoList *myOutput, ObitErr *err);
/* Elliptical ring parameter search */
void searchE (olong nsearch, ofloat spacing, olong prtLv, 
	      ofloat *Center, ofloat *Radius, ofloat *Width, ofloat Frac[2], 
	      ofloat *Ratio, ofloat dRat, ofloat *posAng, ofloat dPA,
	      ObitErr *err);
/* Histograms of new elliptical value  */
void DoHistoE(ofloat cen[2], ofloat ar[2], ofloat pa, ofloat spacing);

/* Program globals */
gchar *pgmName = "Ringer";       /* Program name */
gchar *infile  = "Ringer.inp";   /* File with program inputs */
gchar *outfile = "Ringer.out";   /* File to contain program outputs */
olong  pgmNumber;       /* Program number (like POPS no.) */
olong  AIPSuser;        /* AIPS user number number (like POPS no.) */
olong  nAIPS=0;         /* Number of AIPS directories */
gchar **AIPSdirs=NULL; /* List of AIPS data directories */
olong  nFITS=0;         /* Number of FITS directories */
gchar **FITSdirs=NULL; /* List of FITS data directories */
ObitInfoList *myInput  = NULL; /* Input parameter list */
ObitInfoList *myOutput = NULL; /* Output parameter list */
/* Data arrays */
ofloat *AXpixel=NULL, *AYpixel=NULL, *AData=NULL;
ofloat *BXpixel=NULL, *BYpixel=NULL, *BData=NULL;
olong Acount=0, Bcount=0;
/* histogram storage */
olong nhisto=0;
ofloat *Ahisto=NULL, *Bhisto=NULL;

int main ( int argc, char **argv )
/*----------------------------------------------------------------------- */
/*   Fit center and ring widths to SiO maser images                       */
/*----------------------------------------------------------------------- */
{
  oint ierr = 0;
  ObitSystem   *mySystem= NULL;
  ObitImage    *in1Image= NULL, *in2Image= NULL;
  gboolean     doEllip=FALSE;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitErr      *err= NULL;

   /* Startup - parse command line */
  err = newObitErr();
  myInput = RingerIn (argc, argv, err);
  if (err->error) ierr = 1;
  ObitErrLog(err); /* show any error messages on err */
  if (ierr!=0) goto exit;

  /* Initialize Obit */
  mySystem = ObitSystemStartup (pgmName, pgmNumber, AIPSuser, nAIPS, AIPSdirs, 
				nFITS, FITSdirs, (oint)TRUE, (oint)FALSE, err);
  if (err->error) ierr = 1; ObitErrLog(err); if (ierr!=0) goto exit;

  /* Get input images */
  in1Image = getInputImage1 (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  in2Image = getInputImage2 (myInput, err);
  if (err->error) ierr = 1;  ObitErrLog(err);  if (ierr!=0) goto exit;
  
  /* Check inputs */
  checkInput (in1Image, in2Image, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;

  /* do fitting */
  doFit (myInput, in1Image, in2Image, myOutput, err);
  if (err->error) ierr = 1;  ObitErrLog(err); if (ierr!=0) goto exit;

  /* want elliptical fit as well */
  ObitInfoListGetTest(myInput, "doEllip", &type, dim, &doEllip);
  if (doEllip) {
    doFitE (myInput, in1Image, in2Image, myOutput, err);
    if (err->error) ierr = 1; ObitErrLog(err);  if (ierr!=0) goto exit;
  }
  
  /* cleanup */
  myInput    = ObitInfoListUnref(myInput);    /* delete input list */
  in1Image   = ObitUnref(in1Image);
  in2Image   = ObitUnref(in2Image);
  
  /* Shutdown  */
 exit:
  ObitReturnDumpRetCode (ierr, outfile, myOutput, err);  /* Final output */
  myOutput = ObitInfoListUnref(myOutput);   /* delete output list */
  mySystem = ObitSystemShutdown (mySystem);
  
  return ierr;
} /* end of main */

ObitInfoList* RingerIn (int argc, char **argv, ObitErr *err)
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
  gchar *routine = "RingerIn";

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

  return list;
} /* end RingerIn */

void Usage(void)
/*----------------------------------------------------------------------- */
/*   Tells about usage of program and bails out                           */
/*----------------------------------------------------------------------- */
{
    fprintf(stderr, "Usage: Ringer -input file -output ofile [args]\n");
    fprintf(stderr, "Ringer Obit task = get image mean and RMS\n");
    fprintf(stderr, "Arguments:\n");
    fprintf(stderr, "  -input input parameter file, def Ringer.in\n");
    fprintf(stderr, "  -output output result file, def Ringer.out\n");
    fprintf(stderr, "  -pgmNumber Program (POPS) number, def 1 \n");
    fprintf(stderr, "  -DataType 'AIPS' or 'FITS' type for input image\n");
    fprintf(stderr, "  -inFile input FITS Image file\n");
    fprintf(stderr, "  -AIPSuser User AIPS number, def 2 \n");
    fprintf(stderr, "  -inName input AIPS file name\n");
    fprintf(stderr, "  -inClass input AIPS file class\n");
    fprintf(stderr, "  -inSeq input AIPS file sequence\n");
    fprintf(stderr, "  -inDisk input (AIPS or FITS) disk number (1-rel) \n");
    fprintf(stderr, "  -BLC bottom-left (1-rel) pixel def. [1,1,..] \n");
    fprintf(stderr, "  -TRC top-right (1-rel) pixel def. [nx,ny,..] \n");
    
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
  strTemp = "RingerName";
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
/*     Center    Center Pixel  [0., 0.]                                   */
/*     Width     Width in pixels  [0., 0.]                                */
/*----------------------------------------------------------------------- */
ObitInfoList* defaultOutputs(ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat ftemp[2] = {0.0, 0.0};
  ObitInfoList *out = newObitInfoList();
  gchar *routine = "defaultOutputs";

  /* add parser items */
  /* Center */
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Center", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Radius */
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Radius", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Width */
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Width", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Frac */
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Frac", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* Ratio */
  dim[0] = 2; dim[1] = 1;
  ObitInfoListPut (out, "Ratio", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  /* PosAng  */
  dim[0] = 1; dim[1] = 1;
  ObitInfoListPut (out, "PosAng", OBIT_float, dim, ftemp, err);
  if (err->error) Obit_traceback_val (err, routine, "DefInput", out);

  return out;
} /* end defaultOutputs */

/*----------------------------------------------------------------------- */
/*  Get first input image                                                 */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input image, read and on ObitFArray member        */
/*----------------------------------------------------------------------- */
ObitImage* getInputImage1 (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *inImage = NULL;
  ObitInfoType type;
  olong         Aseq, disk, cno, i;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "getInputImage1";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return inImage;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  inImage = newObitImage("input Image");
  
  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

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
      if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "inSeq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
    /* define object */
    ObitImageSetAIPS (inImage, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
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
    ObitImageSetFITS (inImage, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return inImage;
  }

  /* Ensure inImage fully instantiated and OK */
  ObitImageFullInstantiate (inImage, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);

  /* Set defaults BLC, TRC - use size on myIO as blc, trc incorporated into myDesc */
  for (i=0; i<IM_MAXDIM; i++) {
    if (blc[i]<=0) blc[i] = 1;
    blc[i] = MAX (1,  blc[i]);
    if (trc[i]<=0) trc[i] = ((ObitImageDesc*)inImage->myIO->myDesc)->inaxes[i];
    trc[i] = MIN (trc[i], ((ObitImageDesc*)inImage->myIO->myDesc)->inaxes[i]);
  }

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (inImage->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (inImage->info, "TRC", OBIT_long, dim, trc);

  /* Open and read Image, image on member image, an ObitFArray */
  ObitImageOpen (inImage, OBIT_IO_ReadOnly, err);
  ObitImageRead (inImage, NULL, err);
  ObitImageClose (inImage, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", inImage);

  return inImage;
} /* end getInputImage1 */

/*----------------------------------------------------------------------- */
/*  Get second input image                                                */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*   Return                                                               */
/*       ObitImage with input image, read and on ObitFArray member        */
/*----------------------------------------------------------------------- */
ObitImage* getInputImage2 (ObitInfoList *myInput, ObitErr *err)
{
  ObitImage    *in2Image = NULL;
  ObitInfoType type;
  olong         Aseq, disk, cno, i;
  gchar        *Type, *strTemp, inFile[129];
  gchar        Aname[13], Aclass[7], *Atype = "MA";
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong         blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong         trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "getInputImage2";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return in2Image;
  g_assert (ObitInfoListIsA(myInput));

  /* Create basic input Image data Object */
  in2Image = newObitImage("input Image");
  
  /* Get region from myInput */
  ObitInfoListGetTest(myInput, "BLC", &type, dim, blc); /* BLC */
  ObitInfoListGetTest(myInput, "TRC", &type, dim, trc); /* TRC */

  /* File type - could be either AIPS or FITS */
  ObitInfoListGetP (myInput, "DataType", &type, dim, (gpointer)&Type);
  if (!strncmp (Type, "AIPS", 4)) { /* AIPS input */
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
      strncpy (Aclass, strTemp, 7);
    } else { /* Didn't find */
      strncpy (Aclass, "NoClas", 7);
    }
    Aclass[6] = 0;
    /* input AIPS sequence */
    ObitInfoListGet(myInput, "in2Seq", &type, dim, &Aseq, err);

    /* if ASeq==0 want highest existing sequence */
    if (Aseq<=0) {
      Aseq = ObitAIPSDirHiSeq(disk, AIPSuser, Aname, Aclass, Atype, TRUE, err);
      if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);
      /* Save on myInput*/
      dim[0] = dim[1] = 1;
      ObitInfoListAlwaysPut(myInput, "in2Seq", OBIT_oint, dim, &Aseq);
    } 

    /* Find catalog number */
    cno = ObitAIPSDirFindCNO(disk, AIPSuser, Aname, Aclass, Atype, Aseq, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);
    
    /* define object */
    ObitImageSetAIPS (in2Image, OBIT_IO_byPlane, disk, cno, AIPSuser, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);
    
  } else if (!strncmp (Type, "FITS", 4)) {  /* FITS input */
    /* input FITS file name */
    if (ObitInfoListGetP(myInput, "in2File", &type, dim, (gpointer)&strTemp)) {
      strncpy (inFile, strTemp, 128);
    } else { 
      strncpy (inFile, "No_Filename_Given", 128);
    }
    
    /* input FITS disk */
    ObitInfoListGet(myInput, "in2Disk", &type, dim, &disk, err);

    /* define object */
    ObitImageSetFITS (in2Image, OBIT_IO_byPlane, disk, inFile, blc, trc, err);
    if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);
    
  } else { /* Unknown type - barf and bail */
    Obit_log_error(err, OBIT_Error, "%s: Unknown Data type %s", 
                   pgmName, Type);
    return in2Image;
  }

  /* Ensure in2Image fully instantiated and OK */
  ObitImageFullInstantiate (in2Image, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);

  /* Set defaults BLC, TRC - use size on myIO as blc, trc incorporated into myDesc */
  for (i=0; i<IM_MAXDIM; i++) {
    if (blc[i]<=0) blc[i] = 1;
    blc[i] = MAX (1,  blc[i]);
    if (trc[i]<=0) trc[i] = ((ObitImageDesc*)in2Image->myIO->myDesc)->inaxes[i];
    trc[i] = MIN (trc[i], ((ObitImageDesc*)in2Image->myIO->myDesc)->inaxes[i]);
  }

  /* Save blc, trc */
  dim[0] = IM_MAXDIM;
  ObitInfoListAlwaysPut (myInput, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (in2Image->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut (myInput, "TRC", OBIT_long, dim, trc);
  ObitInfoListAlwaysPut (in2Image->info, "TRC", OBIT_long, dim, trc);

  /* Open and read Image, image on member image, an ObitFArray */
  ObitImageOpen (in2Image, OBIT_IO_ReadOnly, err);
  ObitImageRead (in2Image, NULL, err);
  ObitImageClose (in2Image, err);
  if (err->error) Obit_traceback_val (err, routine, "myInput", in2Image);

  return in2Image;
} /* end getInputImage2 */

/*----------------------------------------------------------------------- */
/*  Check images for compatibility                                        */
/*   Input:                                                               */
/*      in1Image  First image                                             */
/*      in2Image  Second image                                            */
/*   Output:                                                              */
/*      err    Obit Error stack                                           */
/*----------------------------------------------------------------------- */
void checkInput (ObitImage* in1Image, ObitImage* in2Image, ObitErr *err)
{
  ObitImageDesc *in1Desc=NULL, *in2Desc=NULL;
  gchar *routine = "checkInput";

  in1Desc = in1Image->myDesc;
  in2Desc = in2Image->myDesc;
  /* Check size of planes */
  Obit_return_if_fail(((in1Desc->inaxes[0]==in2Desc->inaxes[0]) && 
		       (in1Desc->inaxes[1]==in2Desc->inaxes[1])), err,
		      "%s: Image planes incompatible  %d!= %d or  %d!= %d", 
		      routine, in1Desc->inaxes[0], in2Desc->inaxes[0], 
		      in1Desc->inaxes[1], in2Desc->inaxes[1]);
 
} /* end  checkInput */

/*----------------------------------------------------------------------- */
/*  Do fitting, fitted values saved as "Center" and "Width" on myOutput   */
/*  Also save fraction of flux in anulus as Frac                          */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      in1Image  First image                                             */
/*      in2Image  Second image                                            */
/*   Output:                                                              */
/*      myOutput  Output parameters on InfoList , fitted values added     */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void doFit (ObitInfoList *myInput, ObitImage* in1Image, ObitImage* in2Image, 
	    ObitInfoList *myOutput, ObitErr *err)
{
  ofloat Spacing=1.0, Cutoff=0.0;
  olong NSearch = 51, srch;
  ofloat Radius[2]={-1.0, -1.0}, Center[2]={-1.0, -1.0}, Width[2]={-1.0, -1.0};
  ofloat Frac[2]={-1.0, -1.0};
  ofloat fblank = ObitMagicF(), spc, sum1, sum2;
  ObitFArray *Apix=in1Image->image, *Bpix=in2Image->image; 
  olong prtLv=0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat scale;
  olong i, nx, ix, iy;
  gchar *routine = "doFit";

  if (err->error) return; /* previous error? */

  /* Get parameters from myInput */
  ObitInfoListGetTest(myInput, "Spacing", &type, dim, &Spacing);
  if (Spacing<=0.01) Spacing = 2.0;
  ObitInfoListGetTest(myInput, "NSearch", &type, dim, &NSearch);
  if (NSearch<=0) NSearch = 51;
  ObitInfoListGetTest(myInput, "Cutoff",  &type, dim, &Cutoff);
  Cutoff = MAX (0.01, Cutoff);
  ObitInfoListGetTest(myInput, "Center",  &type, dim, Center);
  /* Center 0-rel */
  Center[0] -= 1.0;
  Center[1] -= 1.0;
  ObitInfoListGetTest(myInput, "Radius",  &type, dim, Radius);
  ObitInfoListGetTest(myInput, "Width",   &type, dim, Width);
  ObitInfoListGetTest(myInput, "prtLv",   &type, dim, &prtLv);

  /* Convert to lists of pixels above Cutoff - first count */
  Acount = 0;
  for (i=0; i<Apix->arraySize; i++) 
    if ((Apix->array[i]!=fblank)&&(Apix->array[i]>=Cutoff)) Acount++;
  Bcount = 0;
  for (i=0; i<Bpix->arraySize; i++) 
    if ((Bpix->array[i]!=fblank)&&(Bpix->array[i]>=Cutoff)) Bcount++;

  /* Better have found more than 10 each */
  Obit_return_if_fail(((Acount+Bcount)>20), err,
		      "%s: Insufficient counts above Cutoff %d %d", 
		      routine, Acount, Bcount);

  AXpixel = g_malloc0(Acount*sizeof(ofloat));
  AYpixel = g_malloc0(Acount*sizeof(ofloat));
  AData   = g_malloc0(Acount*sizeof(ofloat));
  BXpixel = g_malloc0(Bcount*sizeof(ofloat));
  BYpixel = g_malloc0(Bcount*sizeof(ofloat));
  BData   = g_malloc0(Bcount*sizeof(ofloat));

  /* Histogram */ 
  nx = in1Image->myDesc->inaxes[0];
  nhisto = nx;
  Ahisto = g_malloc0 (nx*sizeof(float));
  Bhisto = g_malloc0 (nx*sizeof(float));
  
  /* Generate lists */
  Acount = 0;
  for (i=0; i<Apix->arraySize; i++) {
    if ((Apix->array[i]!=fblank)&&(Apix->array[i]>=Cutoff)) {
      AData[Acount] = Apix->array[i];
      iy = i / nx;
      ix = i - iy * nx;
      AXpixel[Acount] = (ofloat)ix;
      AYpixel[Acount] = (ofloat)iy;
      Acount++;
    }
  }
  
  Bcount = 0;
  for (i=0; i<Bpix->arraySize; i++) {
    if ((Bpix->array[i]!=fblank)&&(Bpix->array[i]>=Cutoff)) {
      BData[Bcount] = Bpix->array[i];
      iy = i / nx;
      ix = i - iy * nx;
      BXpixel[Bcount] = (ofloat)ix;
      BYpixel[Bcount] = (ofloat)iy;
      Bcount++;
    }
  }
  
  /* Default center : Use centroid */
  if (Center[0]<=0.01) {
    sum1 = sum2 = 0.0;
    for (i=0; i<Acount; i++) {sum1+=AData[i]; sum2+=AData[i]*AXpixel[i];}
    for (i=0; i<Bcount; i++) {sum1+=BData[i]; sum2+=BData[i]*BXpixel[i];}
    Center[0] = MAX (0.0, MIN((sum2/sum1), in1Image->myDesc->inaxes[0]-1.0));
  }
  if (Center[1]<=0.01) {
    sum1 = sum2 = 0.0;
    for (i=0; i<Acount; i++) {sum1+=AData[i]; sum2+=AData[i]*AXpixel[i];}
    for (i=0; i<Bcount; i++) {sum1+=BData[i]; sum2+=BData[i]*BXpixel[i];}
    Center[1] = MAX (0.0, MIN((sum2/sum1), in1Image->myDesc->inaxes[1]-1.0));
  }
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Initial Center %5.2f %5.2f ",
		   Center[0]+1.0, Center[1]+1.0);
    ObitErrLog(err); 
  }

  /* do coarse search  */
  spc = MAX (1.0, Spacing*2);
  srch = 1+2*(NSearch/2);   /* Make odd */
  search (spc, srch, Spacing, prtLv, Center, Radius, Width, Frac, err);
  if (err->error) Obit_traceback_msg (err, routine, "search");

  /* do fine search */
  spc = MAX (1.0, Spacing);
  srch = 1+2*((NSearch/2)/2);
  search (spc, srch, Spacing,  prtLv, Center, Radius, Width, Frac, err);
  if (err->error) Obit_traceback_msg (err, routine,"search");

  /* do fine2 search */
  spc = MAX (1.0, Spacing/2);
  srch = 1+2*((NSearch/4)/2);
  search (spc, srch, Spacing, prtLv, Center, Radius, Width, Frac, err);
  if (err->error) Obit_traceback_msg (err, routine, "search");

  /* Center 1-rel */
  Center[0] += 1.0;
  Center[1] += 1.0;

  /* report results */
  scale = fabs(in1Image->myDesc->cdelt[1])*3600000.0;
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Center %5.2f %5.2f Radius %5.2f %5.2f Width %5.2f %5.2f pixels",
		   Center[0], Center[1], Radius[0], Radius[1], Width[0], Width[1]);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Radius %5.2f %5.2f Width %5.2f %5.2f mas",
		   Radius[0]*scale, Radius[1]*scale, 
		   Width[0]*scale, Width[1]*scale);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Fraction of flux in annulus %5.3f %5.3f ",
		   Frac[0], Frac[1]);
 }

  if (prtLv>=2) {
    /* show histograms */
    scale = fabs(in1Image->myDesc->cdelt[1])*3600000.0*Spacing;
    Obit_log_error(err, OBIT_InfoErr, "  Radius   A Hist   B Hist");
    for (i = 0; i<200; i++) {
     Obit_log_error(err, OBIT_InfoErr, 
		    "%8.3f %8.3f %8.3f", 
		    i*scale, 100.0*Ahisto[i], 100.0*Bhisto[i]);
    } 
  }
  
  ObitErrLog(err); /* show any error messages on err */
  /* Cleanup */
  if (AXpixel) g_free(AXpixel);
  if (AYpixel) g_free(AYpixel);
  if (AData)   g_free(AData);
  if (BXpixel) g_free(BXpixel);
  if (BYpixel) g_free(BYpixel);
  if (BData)   g_free(BData);
  if (Ahisto)  g_free(Ahisto);
  if (Bhisto)  g_free(Bhisto);

  /* Save results on output */
  dim[0] = 2;
  ObitInfoListAlwaysPut (myOutput, "Center", OBIT_float, dim, Center);
  ObitInfoListAlwaysPut (myOutput, "Radius", OBIT_float, dim, Radius);
  ObitInfoListAlwaysPut (myOutput, "Width",  OBIT_float, dim, Width);
  ObitInfoListAlwaysPut (myOutput, "Frac",   OBIT_float, dim, Frac);
} /* end doFit */

  void search (ofloat space, olong nsearch, ofloat spacing, olong prtLv, 
	       ofloat *Center, ofloat *Radius, ofloat *Width, ofloat Frac[2], 
	       ObitErr *err)
/*----------------------------------------------------------------------- */
/*   Fit ring                                                             */
/*   Inputs:                                                              */
/*      space    Spacing (pixels) for search                              */
/*      nsearch  Number of searches per parameter                         */
/*      spacing  Spacing (pixels) for histogram                           */
/*      prtLv    print level, debug if >=3                                */
/*   Output:                                                              */
/*      Center   Center pixel                                             */
/*      Radius   Ring radii in pixels                                     */
/*      Width    Ring widths in pixel                                     */
/*      Frac     Fraction of flux in annulus defined by Amom, Bmom        */
/*      err      Obit Error stack                                         */
/*----------------------------------------------------------------------- */
{
  ofloat Amom[2], Bmom[2], cen[2], minxy[2], minval, v, dx, dy;
  ofloat r1, r2;
  olong i, j, nx, ny;

  minxy[0] = minxy[1] = 0;
  /* search grid */
  minval = 1.0e20;
  nx = nsearch; ny = nsearch;
  dx = space; dy = space;
  for (j=0; j<nx; j++) {
    cen[1] = Center[1]+(j-ny/2)*dy;
    for (i=0; i<ny; i++) {
      cen[0] = Center[0]+(i-nx/2)*dx;
      DoHisto(cen, spacing);
      DoMom (FALSE, Amom, Bmom);
      v = Amom[1]+Bmom[1];
      /*      map[i][j] = v;*/
      if (minval>v) {
	minval = v;
	minxy[0] = cen[0];
	minxy[1] = cen[1];
	r1 = Amom[0];
	r2 = Bmom[0];
      }
    }
  }

  /* save results */
  Center[0] = minxy[0];
  Center[1] = minxy[1];
  DoHisto(Center, spacing);
  DoMom (TRUE, Amom, Bmom);
  Radius[0] = Amom[0]*spacing;
  Radius[1] = Bmom[0]*spacing;
  Width[0]  = Amom[1]*spacing;
  Width[1]  = Bmom[1]*spacing;

  /* Get fraction of flux in annulus */
  GetFrac(Amom, Bmom, Frac);

  if (prtLv>=3) {
    Obit_log_error(err, OBIT_InfoErr, 
		   " min %5.3f @ %5.2f %5.2f radii %5.2f %5.2f width %5.2f %5.2f fract %5.3f %5.3f \n",
		   minval, minxy[0]+1.0, minxy[1]+1.0, Radius[0], Radius[1], 
		   Width[0], Width[1], Frac[0], Frac[1]);
    ObitErrLog(err); /* show any error messages on err */
  }
} /* end  search */

void DoHisto(ofloat cen[2], ofloat spacing)
/*----------------------------------------------------------------------- */
/*   Calculate histograms of summed flux vs distance from center          */
/*   Normalize by area                                                    */
/*   Input :                                                              */
/*      cen  R[*] center pixel about which to make the histogram          */
/*      spacing  Spacing (pixels) for histogram                           */
/*----------------------------------------------------------------------- */
{
  olong  i, icell;
  ofloat dx, dy, dist, maxd2, ispace=1.0/spacing;

  /* initialize */
  for (i=0; i<nhisto; i++) {Ahisto[i] = 0.0; Bhisto[i] = 0.0;}

  /* maximum distance squared */
  maxd2 = (nhisto*spacing) * (nhisto*spacing);

  /* A data */
  for (i=0; i<Acount; i++) {
    dx = AXpixel[i]-cen[0];
    dy = AYpixel[i]-cen[1];
    dist = dx*dx + dy*dy;
    if (dist<maxd2) {
      dist = sqrt (dist);
      icell = dist * ispace + 0.5;
      if (icell>=nhisto) icell = nhisto-1;
      Ahisto[icell] += AData[i];
    }
  }
  /* normalize by area */
  for (i=1; i<nhisto; i++) Ahisto[i] /= (float)i;

  /* B data */
  for (i=0; i<Bcount; i++) {
    dx = BXpixel[i]-cen[0];
    dy = BYpixel[i]-cen[1];
    dist = dx*dx + dy*dy;
    if (dist<maxd2) {
      dist = sqrt (dist);
      icell = dist * ispace + 0.5;
      if (icell>=nhisto) icell = nhisto-1;
      Bhisto[icell] += BData[i];
    }
  }
 /* normalize by area */
  for (i=1; i<nhisto; i++) Bhisto[i] /= (float)i;
} /* end  DoHisto */

void DoMom(gboolean final, ofloat Amom[2], ofloat Bmom[2])
/*----------------------------------------------------------------------- */
/*   Get Moments of histograms                                            */
/*   Input:                                                               */
/*      final L    if true, restrict to mom1 -> mom0+2mom1                */
/*   Output :                                                             */
/*      Amom  R[*] moments of histogram, second and higher about the first*/
/*      Bmom  R[*] moments of histogram, second and higher about the first*/
/*----------------------------------------------------------------------- */
{
  olong  i, imaxx, iminn;
  odouble dx, sumx, sumxy;

  /* A image */
  /* first moment */
  sumx = 0.0; sumxy = 0.0;
  for (i=0; i<nhisto; i++) {
    sumx  += Ahisto[i];
    sumxy += (odouble)i * Ahisto[i];
  }
  Amom[0] = sumxy / sumx;

  /* second moment */
  sumxy = 0.0;
  for (i=0; i<nhisto; i++) {
    dx = (odouble)i - Amom[0];
    sumxy += dx * dx * Ahisto[i];
  }
  Amom[1] = sqrt (sumxy / sumx);

  if (final) {
    /* Do it again excluding points closer than mom0-mom1 or
       further than mom0+2Mom1 */
    iminn = 0;
    imaxx = MIN( nhisto-1,(olong)((Amom[0] + 2.0 * Amom[1]) + 0.5));
    /* first moment */
    sumx = 0.0; sumxy = 0.0;
    for (i=iminn; i<=imaxx; i++) {
      sumx  += Ahisto[i];
      sumxy += (odouble)i * Ahisto[i];
    }
    Amom[0] = sumxy / sumx;
    
    /* second moment */
    sumxy = 0.0;
    for (i=0; i<nhisto; i++) {
      if (i>imaxx) break;
      dx = (odouble)i - Amom[0];
      sumxy += dx * dx * Ahisto[i];
    }
    Amom[1] = sqrt (sumxy / sumx);
  }

  /* B image */
  /* first moment */
  sumx = 0.0; sumxy = 0.0;
  for (i=0; i<nhisto; i++) {
    sumx  += Bhisto[i];
    sumxy += (odouble)i * Bhisto[i];
  }
  Bmom[0] = sumxy / sumx;

  /* second moment */
  sumxy = 0.0;
  for (i=0; i<nhisto; i++) {
    dx = (odouble)i - Bmom[0];
    sumxy += dx * dx * Bhisto[i];
  }
  Bmom[1] = sqrt (sumxy / sumx);

  /* Do it again excluding points further than mom0+Mom1 */
  if (final) {
    iminn = 0;
    imaxx = MIN( nhisto-1,(olong)((Bmom[0] + 2.0 * Bmom[1]) + 0.5));
    /* first moment */
    sumx = 0.0; sumxy = 0.0;
    for (i=iminn; i<=imaxx; i++) {
      sumx  += Bhisto[i];
      sumxy += (odouble)i * Bhisto[i];
    }
    Bmom[0] = sumxy / sumx;
    
    /* second moment */
    sumxy = 0.0;
    for (i=0; i<nhisto; i++) {
      if (i>imaxx) break;
      dx = (odouble)i - Bmom[0];
      sumxy += dx * dx * Bhisto[i];
    }
    Bmom[1] = sqrt (sumxy / sumx);
  }
} /* end  DoMom */

void GetFrac(ofloat Amom[2], ofloat Bmom[2], ofloat Frac[2])
/*----------------------------------------------------------------------- */
/*   Calculate histograms of summed flux vs distance from center          */
/*   Normalize by area                                                    */
/*   Input :                                                              */
/*      Amom     Moments of histogram, second and higher about the first  */
/*      Bmom     Moments of histogram, second and higher about the first  */
/*   Output:                                                              */
/*      Frac     Fraction of flux in annulus defined by Amom, Bmom        */
/*----------------------------------------------------------------------- */
{
  olong  i, i1, i2, icen, iwid;
  ofloat sum, total;

  /* A histogram */
  total = 0.0;
  for (i=0; i<nhisto; i++) total +=Ahisto[i];
  /* Boundaries of annulus */
  icen = (olong)(Amom[0]+0.5);
  iwid = (olong)(Amom[1]+0.5);
  i1 = MAX (0, icen-iwid);
  i2 = MIN (nhisto-1, icen+iwid);
  sum = 0.0;
  for (i=i1; i<=i2; i++)   sum +=Ahisto[i];
  Frac[0] = sum/total;

  /* B histogram */
  total = 0.0;
  for (i=0; i<nhisto; i++) total +=Bhisto[i];
  /* Boundaries of annulus */
  icen = (olong)(Bmom[0]+0.5);
  iwid = (olong)(Bmom[1]+0.5);
  i1 = MAX (0, icen-iwid);
  i2 = MIN (nhisto-1, icen+iwid);
  sum = 0.0;
  for (i=i1; i<=i2; i++)   sum +=Bhisto[i];
  Frac[1] = sum/total;

} /* end  GetFrac */

/*----------------------------------------------------------------------- */
/*  Fit axial ratioes and position angle starting from Center and Radius  */
/*                          */
/*   Input:                                                               */
/*      myInput   Input parameters on InfoList                            */
/*      in1Image  First image                                             */
/*      in2Image  Second image                                            */
/*   Output:                                                              */
/*      myOutput  Output parameters on InfoList , fitted values added     */
/*      err       Obit Error stack                                        */
/*----------------------------------------------------------------------- */
void doFitE (ObitInfoList *myInput, ObitImage* in1Image, ObitImage* in2Image, 
	    ObitInfoList *myOutput, ObitErr *err)
{
  ofloat Spacing=1.0, Cutoff=0.0;
  olong NSearch = 51, srch;
  ofloat Radius[2]={-1.0, -1.0}, Center[2]={-1.0, -1.0}, Width[2]={-1.0, -1.0};
  ofloat Ratio[2]={1.0,1.0},  Frac[2], PosAng=0.0, dRat, dPA, PARange;
  ofloat fblank = ObitMagicF();
  ObitFArray *Apix=in1Image->image, *Bpix=in2Image->image; 
  olong prtLv=0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat scale;
  olong i, nx, ix, iy;
  gchar *routine = "doFitE";

  if (err->error) return; /* previous error? */

  /* Get parameters from myInput/myOutput */
  ObitInfoListGetTest(myInput, "Spacing", &type, dim, &Spacing);
  if (Spacing<=0.01) Spacing = 2.0;
  ObitInfoListGetTest(myInput, "NSearch", &type, dim, &NSearch);
  if (NSearch<=0) NSearch = 51;
  ObitInfoListGetTest(myInput, "Cutoff",  &type, dim, &Cutoff);
  Cutoff = MAX (0.01, Cutoff);
  ObitInfoListGetTest(myInput, "Ratio",   &type, dim, Ratio);
  ObitInfoListGetTest(myInput, "PosAng",  &type, dim, &PosAng);
  ObitInfoListGetTest(myOutput, "Center",  &type, dim, Center);
  /* Center 0-rel */
  Center[0] -= 1.0;
  Center[1] -= 1.0;
  ObitInfoListGetTest(myOutput, "Radius",  &type, dim, Radius);
  ObitInfoListGetTest(myOutput, "Width",   &type, dim, Width);
  ObitInfoListGetTest(myInput,  "prtLv",   &type, dim, &prtLv);

  /* Convert to lists of pixels above Cutoff - first count */
  Acount = 0;
  for (i=0; i<Apix->arraySize; i++) 
    if ((Apix->array[i]!=fblank)&&(Apix->array[i]>=Cutoff)) Acount++;
  Bcount = 0;
  for (i=0; i<Bpix->arraySize; i++) 
    if ((Bpix->array[i]!=fblank)&&(Bpix->array[i]>=Cutoff)) Bcount++;

  /* Better have found more than 10 each */
  Obit_return_if_fail(((Acount+Bcount)>20), err,
		      "%s: Insufficient counts above Cutoff %d %d", 
		      routine, Acount, Bcount);

  AXpixel = g_malloc0(Acount*sizeof(ofloat));
  AYpixel = g_malloc0(Acount*sizeof(ofloat));
  AData   = g_malloc0(Acount*sizeof(ofloat));
  BXpixel = g_malloc0(Bcount*sizeof(ofloat));
  BYpixel = g_malloc0(Bcount*sizeof(ofloat));
  BData   = g_malloc0(Bcount*sizeof(ofloat));

  /* Histogram */ 
  nx = in1Image->myDesc->inaxes[0];
  nhisto = nx;
  Ahisto = g_malloc0 (nx*sizeof(float));
  Bhisto = g_malloc0 (nx*sizeof(float));
  
  /* Generate lists */
  Acount = 0;
  for (i=0; i<Apix->arraySize; i++) {
    if ((Apix->array[i]!=fblank)&&(Apix->array[i]>=Cutoff)) {
      AData[Acount] = Apix->array[i];
      iy = i / nx;
      ix = i - iy * nx;
      AXpixel[Acount] = (ofloat)ix;
      AYpixel[Acount] = (ofloat)iy;
      Acount++;
    }
  }
  
  Bcount = 0;
  for (i=0; i<Bpix->arraySize; i++) {
    if ((Bpix->array[i]!=fblank)&&(Bpix->array[i]>=Cutoff)) {
      BData[Bcount] = Bpix->array[i];
      iy = i / nx;
      ix = i - iy * nx;
      BXpixel[Bcount] = (ofloat)ix;
      BYpixel[Bcount] = (ofloat)iy;
      Bcount++;
    }
  }
  
  /* Position angle to Radians */
  PosAng *= DG2RAD;
  /* If position angle given only search +/- 20 deg */
  if (fabs(PosAng)>1.0e-5) PARange = 20.0*2.0*DG2RAD;
  else PARange = 3.1415926;

  /* Default ratio - 1.5 */ 
  if (Ratio[0]<=0.0) Ratio[0] = 1.5;
  if (Ratio[1]<=0.0) Ratio[1] = 1.5;
  
  /* do coarse search  */
  srch = 1+2*(NSearch/2);   /* Make odd */
  srch = MAX (5, srch);
  dRat = 1.0/srch;
  dPA  = PARange/srch;
  searchE (srch, Spacing, prtLv, Center, Radius, Width, Frac, 
	   Ratio, dRat, &PosAng, dPA, err);
  if (err->error) Obit_traceback_msg (err, routine, "search");

  /* do fine search */
  srch = 1+2*((NSearch/2)/2);
  srch = MAX (5, srch);
  dRat = 0.05/srch;
  dPA  = 0.25*PARange/srch;
  searchE (srch, Spacing,  prtLv, Center, Radius, Width, Frac,  
	   Ratio, dRat, &PosAng, dPA, err);
  if (err->error) Obit_traceback_msg (err, routine,"search");

  /* do fine2 search */
  srch = 1+2*((NSearch/4)/2);
  srch = MAX (5, srch);
  dRat = 0.01/srch;
  dPA  = 0.05*PARange/srch;
  searchE (srch, Spacing, prtLv, Center, Radius, Width, Frac,  
	   Ratio, dRat, &PosAng, dPA, err);
  if (err->error) Obit_traceback_msg (err, routine, "search");

  /* Center 1-rel */
  Center[0] += 1.0;
  Center[1] += 1.0;

  /* Position angle to degrees */
  PosAng *= RAD2DG;

  /* report results */
  scale = fabs(in1Image->myDesc->cdelt[1])*3600000.0;
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "Axial Ratio %5.2f %5.2f Position angle %5.2f deg",
		   Ratio[0], Ratio[1], PosAng);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Center %5.2f %5.2f Radius %5.2f %5.2f Width %5.2f %5.2f pixels",
		   Center[0], Center[1], Radius[0], Radius[1], Width[0], Width[1]);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Radius %5.2f %5.2f Width %5.2f %5.2f mas",
		   Radius[0]*scale, Radius[1]*scale, 
		   Width[0]*scale, Width[1]*scale);
    Obit_log_error(err, OBIT_InfoErr, 
		   "Fraction of flux in annulus %5.3f %5.3f ",
		   Frac[0], Frac[1]);
 }

  if (prtLv>=2) {
    /* show histograms */
    scale = fabs(in1Image->myDesc->cdelt[1])*3600000.0*Spacing;
    Obit_log_error(err, OBIT_InfoErr, "  Radius   A Hist   B Hist");
    for (i = 0; i<200; i++) {
     Obit_log_error(err, OBIT_InfoErr, 
		    "%8.3f %8.3f %8.3f", 
		    i*scale, 100.0*Ahisto[i], 100.0*Bhisto[i]);
    } 
  }
  
  ObitErrLog(err); /* show any error messages on err */
  /* Cleanup */
  if (AXpixel) g_free(AXpixel);
  if (AYpixel) g_free(AYpixel);
  if (AData)   g_free(AData);
  if (BXpixel) g_free(BXpixel);
  if (BYpixel) g_free(BYpixel);
  if (BData)   g_free(BData);
  if (Ahisto)  g_free(Ahisto);
  if (Bhisto)  g_free(Bhisto);

  /* Save results on output */
  dim[0] = 2;
  ObitInfoListAlwaysPut (myOutput, "Ratio",  OBIT_float, dim, Ratio);
  ObitInfoListAlwaysPut (myOutput, "Radius", OBIT_float, dim, Radius);
  ObitInfoListAlwaysPut (myOutput, "Width",  OBIT_float, dim, Width);
  ObitInfoListAlwaysPut (myOutput, "Frac",   OBIT_float, dim, Frac);
  dim[0] = 1;
  ObitInfoListAlwaysPut (myOutput, "PosAng", OBIT_float, dim, &PosAng);
} /* end doFitE */

void searchE (olong nsearch, ofloat spacing, olong prtLv, 
	      ofloat *Center, ofloat *Radius, ofloat *Width, ofloat Frac[2],
	      ofloat *Ratio, ofloat dRat, ofloat *posAng, ofloat dPA,
	      ObitErr *err)
/*----------------------------------------------------------------------- */
/*   Fit Axial ratio and position angle                                   */
/*   Does 3 parameter direct search                                       */
/*   Inputs:                                                              */
/*      nsearch  Number of searches per parameter                         */
/*      spacing  Spacing (pixels) for histogram                           */
/*      prtLv    print level, debug if >=3                                */
/*      Center   Center pixel                                             */
/*      dRat     Spacing in ratio search                                  */
/*      dPA      Spacing in position angle search (rad)                   */
/*   Output:                                                              */
/*      Radius   Ring radii in pixels                                     */
/*      Width    Ring widths in pixel                                     */
/*      Ratio    Axial ratio per image, center value on input             */
/*      Frac     Fraction of flux in annulus defined by Amom, Bmom        */
/*      PosAng   Position angle (rad) center value on input               */
/*      err      Obit Error stack                                         */
/*----------------------------------------------------------------------- */
{
  ofloat Amom[2], Bmom[2],  minxy[3], minval, v;
  ofloat ar[2], pa;
  olong i, j, k, ns;

  minxy[0] = minxy[1] = minxy[2] = 0;
  /* search grid */
  minval = 1.0e20;
  ns = 1 + 2*(nsearch/4);
  for (k=0; k<nsearch; k++) {  /* Loop in PA */
    pa = (*posAng) + (k-nsearch/2)*dPA;
    for (j=0; j<ns; j++) { /* Loop in first ratio */
      ar[0] = Ratio[0] + (j-nsearch/2)*dRat;
      if (ar[0]<1.0) continue;
      for (i=0; i<ns; i++) { /* Loop in second ratio */
	ar[1] = Ratio[1] + (i-nsearch/2)*dRat;
	if (ar[1]<1.0) continue;
	DoHistoE(Center, ar, pa, spacing);
	DoMom (FALSE, Amom, Bmom);
	v = Amom[1]+Bmom[1];
	/*      map[i][j] = v;*/
	if (minval>v) {
	  minval = v;
	  minxy[0] = pa;     /* Save min pa */
	  minxy[1] = ar[0];  /* Save min first axial ratio */
	  minxy[2] = ar[1];  /* Save min second axial ratio */
	}
      } /* end Loop in second ratio */
    } /* end Loop in first ratio */
  } /* end PA Loop */

  /* best value statistics */
  DoHistoE(Center, &minxy[1], minxy[0], spacing);
  DoMom (TRUE, Amom, Bmom);
  Radius[0] = Amom[0]*spacing;
  Radius[1] = Bmom[0]*spacing;
  Width[0]  = Amom[1]*spacing;
  Width[1]  = Bmom[1]*spacing;
  
  /* Get fraction of flux in annulus */
  GetFrac(Amom, Bmom, Frac);
  
  /* save results */
  *posAng  = minxy[0];
  Ratio[0] = minxy[1];
  Ratio[1] = minxy[2];
  
  if (prtLv>=3) {
    Obit_log_error(err, OBIT_InfoErr, 
		   " min %5.3f Ratio %5.2f %5.2f PosAngle %5.2f ",
		   minval, Ratio[0], Ratio[1], (*posAng)*57.296);
    ObitErrLog(err); /* show any messages on err */
  }
} /* end  searchE */

void DoHistoE(ofloat cen[2], ofloat ar[2], ofloat pa, ofloat spacing)
/*----------------------------------------------------------------------- */
/*   Calculate histograms of summed flux vs distance from center          */
/*   Normalize by area                                                    */
/*   Input :                                                              */
/*      cen  R[2] center pixel about which to make the histogram          */
/*      ar   R[2] Axial ratio of first and second image                   */
/*      pa   Position angle of major axis                                 */
/*      spacing  Spacing (pixels) for histogram                           */
/*----------------------------------------------------------------------- */
{
  olong  i, icell;
  ofloat dx, dy, x2, y2, dist, maxd2, ispace=1.0/spacing;
  ofloat theta, phi, tanphi, maja1, mina1, maja2, mina2, rad;

  /* initialize */
  for (i=0; i<nhisto; i++) {Ahisto[i] = 0.0; Bhisto[i] = 0.0;}

  /* maximum distance squared */
  maxd2 = (nhisto*spacing) * (nhisto*spacing);

  /* Major and minor axes for unit ellipse */
  mina1 = 2.0 / (1.0 + ar[0]);
  maja1 = mina1 * ar[0];
  mina1 *= mina1;  /* Square */
  maja1 *= maja1;
  mina2 = 2.0 / (1.0 + ar[1]);
  maja2 = mina2 * ar[1];
  mina2 *= mina2;  /* Square */
  maja2 *= maja2;

  /* A data */
  for (i=0; i<Acount; i++) {
    dx = AXpixel[i]-cen[0];
    dy = AYpixel[i]-cen[1];
    theta = atan2(-dx, dy);
    /* counter rotation in ellipse frame */
    phi = theta - pa;
    /* distance to center from puncture point along angle phi */
    tanphi = tan(phi);
    y2 = (maja1*mina1/(maja1*tanphi*tanphi + mina1));
    x2 = (mina1 - (y2*mina1)/maja1);
    rad = sqrt(x2 + y2);

    /* Distort by inverse of distance on ellipse from center.*/
    dx /= rad;
    dy /= rad;
    dist = dx*dx + dy*dy;
    if (dist<maxd2) {
      dist = sqrt (dist);
      icell = dist * ispace + 0.5;
      if (icell>=nhisto) icell = nhisto-1;
      Ahisto[icell] += AData[i];
    }
  }
  /* normalize by area */
  for (i=1; i<nhisto; i++) Ahisto[i] /= (float)i;

  /* B data */
  for (i=0; i<Bcount; i++) {
    dx = BXpixel[i]-cen[0];
    dy = BYpixel[i]-cen[1];
    theta = atan2(-dx, dy);
    /* rotation in ellipse frame */
    phi = theta - pa;
    /* distance to center from puncture point along angle phi */
    tanphi = tan(phi);
    y2 = (maja2*mina2/(maja2*tanphi*tanphi + mina2));
    x2 = (mina2 - (y2*mina2)/maja2);
    rad = sqrt(x2 + y2);

    /* Distort by inverse of distance on ellipse from center.*/
    dx /= rad;
    dy /= rad;
    dist = dx*dx + dy*dy;
    if (dist<maxd2) {
      dist = sqrt (dist);
      icell = dist * ispace + 0.5;
      if (icell>=nhisto) icell = nhisto-1;
      Bhisto[icell] += BData[i];
    }
  }
 /* normalize by area */
  for (i=1; i<nhisto; i++) Bhisto[i] /= (float)i;
} /* end  DoHistoE */
